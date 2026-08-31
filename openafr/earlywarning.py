"""Snapshot store + baseline for the C. auris azole early-warning track (issue #21).

Pipeline stage: resistance early-warning / persistence layer.

Why this exists
---------------
Emergence detection (#22) needs memory: to say a mutation or isolate is *new* or
*rising* we must know what we had seen before. This module turns each NCBI
Pathogen Detection pull into a normalized, versioned, byte-stable snapshot and
gives the two operations the rest of the track builds on:

  * baseline_as_of(records, date) -- the set of isolates NCBI had as of a date
  * diff_snapshots(old, new)      -- what isolates are new between two pulls

What the spike (#20) forces on the schema
-----------------------------------------
NCBI runs *no* AMR pipeline on C. auris: AMR_genotypes / computed_types are empty
for every isolate, so there is no resistance call to read off the feed. The
resistance signal is ours to manufacture from the linked SRA reads (the ERG11
re-caller, #23). This store therefore carries the isolate metadata and a
genome pointer now, and an explicitly *pending* resistance field that #23 will
populate. We never fake a call: `erg11_call` stays empty and `resistance_source`
records why (reads available vs. no reads), so an unresolved isolate is honestly
unresolved rather than silently "wild type".

Design decisions
----------------
* The stable isolate identity is `target_acc` with its version suffix stripped
  (PDT000585001.1 -> PDT000585001). NCBI bumps the `.N` when it re-processes the
  same isolate; the base accession is the thing that persists across releases.
* A snapshot data file is *byte-stable per release*: rows are sorted by isolate
  key and the pull timestamp is kept OUT of the data. Re-pulling the same release
  yields an identical file, so `git diff` between releases shows only real change,
  and re-running is idempotent. Pull provenance (timestamp, tool version, source
  URL, content hash) lives in the separate INDEX manifest, not in the rows.
* Stdlib only (urllib + csv + json + hashlib), matching scripts/probe_ncbi_auris.py.
  A persistence layer should not drag pandas into a repo whose dependency surface
  is rdkit + the docking toolkit.
"""
import csv
import datetime as _dt
import hashlib
import io
import re
import urllib.request

BASE = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Candidozyma_auris/"
ORG_PREFIX = "PDG000000067"

# Values that mean "no real data" in BioSample-derived metadata (see probe).
NULLISH = {"", "null", "missing", "not provided", "not collected",
           "not applicable", "n/a", "na", "unknown", "not determined"}

# The snapshot schema. Order is fixed so snapshot files are byte-stable and
# git-diffable. Every isolate row carries exactly these columns.
#
# Schema evolution note: the FKS1/echinocandin early-warning track (v2, T3) appended
# `fks1_call`/`fks1_resistance_source` AFTER the azole-only columns, so a snapshot written
# before v2 lacks them. `read_snapshot` migrates such a snapshot in memory on load
# (`migrate_record`), seeding the two fields to their honest `pending` state from run_acc --
# never a silent "wild type". Appending (not inserting) keeps the pre-v2 columns in their
# original positions so a diff between a pre-v2 and a migrated snapshot shows only the two
# new trailing columns.
SNAPSHOT_COLUMNS = (
    "isolate_key",          # stable identity: target_acc without the .N version
    "target_acc",           # full NCBI target accession, e.g. PDT000585001.1
    "biosample_acc",        # BioSample, e.g. SAMN05379575
    "run_acc",              # SRA run -- the substrate for our ERG11 re-caller (#23)
    "asm_acc",              # assembly accession if NCBI assembled it (rare, ~3%)
    "collection_date",      # when the sample was taken (raw, may be a bare year)
    "target_creation_date", # when the isolate became visible in NCBI (ISO date)
    "country",              # geo_loc_name normalized to country
    "geo_loc_name",         # raw geography (may carry region after a colon)
    "lat_lon",              # coordinates when present (sparse, ~24%)
    "erg11_call",           # azole-resistance substitutions -- PENDING (#23), never faked
    "resistance_source",    # provenance/why-empty for erg11_call
    "fks1_call",            # echinocandin-resistance substitutions -- PENDING (v2/T3), never faked
    "fks1_resistance_source",  # provenance/why-empty for fks1_call
)

# The columns the v2 FKS1 track appended. A snapshot missing these was written before v2;
# migrate_record back-fills them on read (seeded from run_acc, honestly pending).
_V2_FKS1_COLUMNS = ("fks1_call", "fks1_resistance_source")


def _get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "openafr-earlywarning/1"})
    with urllib.request.urlopen(req, timeout=180) as resp:
        return resp.read()


def latest_release():
    """Highest PDG release number currently published for C. auris."""
    html = _get(BASE).decode("utf-8", "replace")
    versions = re.findall(rf'href="{ORG_PREFIX}\.(\d+)/"', html)
    if not versions:
        raise RuntimeError("could not list releases at " + BASE)
    return max(int(v) for v in versions)


def metadata_url(release):
    tag = f"{ORG_PREFIX}.{release}"
    return tag, f"{BASE}{tag}/Metadata/{tag}.metadata.tsv"


def fetch_release_rows(release=None):
    """Fetch one NCBI release's metadata TSV. Returns (rows, tag, url)."""
    if release is None:
        release = latest_release()
    tag, url = metadata_url(release)
    text = _get(url).decode("utf-8", "replace")
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    return rows, tag, url


def _clean(value):
    """Blank out nullish placeholders; strip whitespace on everything else."""
    v = (value or "").strip()
    return "" if v.lower() in NULLISH else v


def _isolate_key(target_acc):
    """PDT000585001.1 -> PDT000585001 (strip the re-processing version suffix)."""
    return target_acc.split(".")[0] if target_acc else ""


def normalize(rows):
    """Turn raw NCBI metadata dicts into snapshot records (list of dicts).

    Rows without a target accession cannot be keyed, so they are dropped and
    reported to the caller instead of being invented into the store.
    Returns (records, dropped) where `dropped` counts unkeyable rows.
    """
    records, dropped = [], 0
    for r in rows:
        target = _clean(r.get("target_acc"))
        key = _isolate_key(target)
        if not key:
            dropped += 1
            continue
        run = _clean(r.get("Run"))
        geo = _clean(r.get("geo_loc_name"))
        country = geo.split(":")[0].strip() if geo else ""
        records.append({
            "isolate_key": key,
            "target_acc": target,
            "biosample_acc": _clean(r.get("biosample_acc")),
            "run_acc": run,
            "asm_acc": _clean(r.get("asm_acc")),
            "collection_date": _clean(r.get("collection_date")),
            "target_creation_date": _clean(r.get("target_creation_date")),
            "country": country,
            "geo_loc_name": geo,
            "lat_lon": _clean(r.get("lat_lon")),
            # NCBI exposes no azole call for C. auris (spike #20). Leave the call
            # empty and record why, so #23 can fill it and #22 never mistakes an
            # uncalled isolate for a susceptible one.
            "erg11_call": "",
            "resistance_source": "pending:sra-recaller" if run else "pending:no-reads",
            # NCBI likewise exposes no echinocandin/FKS1 call. Same discipline as ERG11:
            # empty call, a pending source recording why, filled by the FKS1 re-caller (v2).
            "fks1_call": "",
            "fks1_resistance_source": _pending_fks1_source(run),
        })
    # Deduplicate on isolate_key, keeping the highest target_acc version so a
    # re-processed isolate does not appear twice. NCBI can carry both .1 and .2.
    best = {}
    for rec in records:
        cur = best.get(rec["isolate_key"])
        if cur is None or rec["target_acc"] > cur["target_acc"]:
            best[rec["isolate_key"]] = rec
    out = sorted(best.values(), key=lambda x: x["isolate_key"])
    return out, dropped


def _pending_fks1_source(run_acc):
    """The honest pending FKS1 source for an isolate, seeded from its run_acc exactly as
    the ERG11 source is: a linked SRA run means the re-caller can attempt it
    (`pending:sra-fks1-recaller`); no run means there is nothing to call (`pending:no-reads`)."""
    return "pending:sra-fks1-recaller" if (run_acc or "").strip() else "pending:no-reads"


def migrate_record(rec):
    """Back-fill any v2 schema column a pre-FKS1 snapshot record lacks, IN PLACE, and return it.

    A snapshot written before the FKS1/echinocandin track carries no `fks1_call` /
    `fks1_resistance_source`. Reading one back and writing it (e.g. an ERG11 `fill`) would,
    without this, drop those columns to empty via the DictWriter default -- turning an
    honestly-uncalled FKS1 status into a blank that could later read as "wild type". So on
    load we seed them to the SAME pending state `normalize` would have: an empty call and a
    run_acc-derived pending source. A record that already carries a real FKS1 source (an
    isolate the re-caller has resolved) is left UNTOUCHED -- we never overwrite a real call
    with a pending one."""
    if rec.get("fks1_call") is None:
        rec["fks1_call"] = ""
    if not (rec.get("fks1_resistance_source") or "").strip():
        rec["fks1_resistance_source"] = _pending_fks1_source(rec.get("run_acc"))
    return rec


def _serialize(records):
    """Render records to the canonical, byte-stable snapshot TSV text.

    Each record is migrated (on a copy) first, so a record that predates the FKS1 columns is
    written with its honest pending FKS1 state rather than an empty cell -- which keeps
    write-then-read an identity (an empty cell would migrate back to pending on read and the
    digests would disagree)."""
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=list(SNAPSHOT_COLUMNS),
                       delimiter="\t", lineterminator="\n",
                       extrasaction="ignore")
    w.writeheader()
    for rec in sorted(records, key=lambda x: x["isolate_key"]):
        w.writerow(migrate_record(dict(rec)))
    return buf.getvalue()


def snapshot_sha256(records):
    return hashlib.sha256(_serialize(records).encode("utf-8")).hexdigest()


def write_snapshot(path, records):
    """Write records to `path` as a canonical snapshot. Idempotent: the same
    records always produce the same bytes. Returns the content sha256."""
    text = _serialize(records)
    with open(path, "w", newline="") as fh:
        fh.write(text)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def read_snapshot(path):
    """Load a snapshot TSV back into records (list of dicts), migrating pre-v2 rows.

    A snapshot written before the FKS1 track lacks the two v2 columns; each record is
    passed through `migrate_record` so every returned record carries the full schema with
    an honest pending FKS1 status. For a current-schema snapshot the migration is a no-op,
    so `read_snapshot(write_snapshot(x)) == x` still holds."""
    with open(path, newline="") as fh:
        return [migrate_record(r) for r in csv.DictReader(fh, delimiter="\t")]


def migrate_snapshot_file(path):
    """Rewrite a snapshot file to the current schema in place (a pre-v2 -> v2 migration).

    Idempotent: reads (which migrates in memory), then writes the canonical serialization.
    Returns (sha256, changed) where `changed` is True iff the on-disk bytes moved -- so a
    caller can report "already current" without a spurious rewrite."""
    with open(path, newline="") as fh:
        before = fh.read()
    records = [migrate_record(r) for r in csv.DictReader(io.StringIO(before), delimiter="\t")]
    after = _serialize(records)
    if after != before:
        with open(path, "w", newline="") as fh:
            fh.write(after)
    return hashlib.sha256(after.encode("utf-8")).hexdigest(), after != before


def baseline_as_of(records, date):
    """The baseline: isolates NCBI had made visible on or before `date`.

    Keyed on target_creation_date (sample-visible-in-NCBI), not collection_date,
    because emergence is measured against what we *could have seen* in a prior
    pull. Isolates with no/partial creation date are excluded from the baseline
    (we cannot place them in time) and returned as a count so the gap is visible.
    Returns (in_baseline, undated) where both are lists of records.
    """
    date = str(date)[:10]
    in_baseline, undated = [], []
    for rec in records:
        created = (rec.get("target_creation_date") or "")[:10]
        if len(created) == 10 and created[:4].isdigit():
            if created <= date:
                in_baseline.append(rec)
        else:
            undated.append(rec)
    return in_baseline, undated


def diff_snapshots(old_records, new_records):
    """What changed between an earlier and a later snapshot.

    Returns a dict:
      added   -- records whose isolate_key is in new but not old (emergence input)
      removed -- isolate_keys in old but not new (isolates NCBI withdrew/renamed)
      revised -- isolate_keys in both but with a bumped target_acc version
    """
    old_by = {r["isolate_key"]: r for r in old_records}
    new_by = {r["isolate_key"]: r for r in new_records}
    added = [new_by[k] for k in sorted(new_by.keys() - old_by.keys())]
    removed = sorted(old_by.keys() - new_by.keys())
    revised = [k for k in sorted(new_by.keys() & old_by.keys())
               if new_by[k].get("target_acc") != old_by[k].get("target_acc")]
    return {"added": added, "removed": removed, "revised": revised}


# ----- INDEX manifest: git-tracked provenance without committing big blobs -----

INDEX_COLUMNS = ("release_tag", "pulled_at_utc", "n_isolates", "n_dropped",
                 "sha256", "source_url", "tool_version")
TOOL_VERSION = "1"


def update_index(index_path, entry):
    """Upsert one release's provenance row into the git-tracked INDEX.

    The INDEX is the small, committable memory of the store: which releases were
    pulled, when, how many isolates, and the content hash of each snapshot. The
    multi-MB snapshot TSVs themselves stay out of git (see data/earlywarning).
    Idempotent per release_tag: re-pulling replaces that release's row.
    """
    rows = {}
    try:
        with open(index_path, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rows[r["release_tag"]] = r
    except FileNotFoundError:
        pass
    rows[entry["release_tag"]] = entry
    with open(index_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(INDEX_COLUMNS), delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for tag in sorted(rows, key=lambda t: int(t.split(".")[-1])):
            w.writerow(rows[tag])


def utcnow_iso():
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
