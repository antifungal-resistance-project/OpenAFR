"""Build the issue-#82 TEMPORAL holdout — a genuinely prospective active/decoy set drawn
from ChEMBL documents that POST-DATE the entire tuning corpus, so the frozen geometric
criterion is tested on chemistry that did not exist when the criterion was built.

Pipeline stage: ligands/ (temporal holdout)

The stronger follow-up #82 asks for. The disjoint-population holdout already landed (#108,
scripts/make_external_holdout.py + validate_gate_external.py): a never-*scored* subset of the
SAME curated ChEMBL population. Its one honest residual — "a clean never-touched holdout, but
not a later time slice" — is exactly what this closes. Here the split is TEMPORAL, not just
never-scored:

  every molecule's earliest antifungal / CYP51 activity is reported in a ChEMBL document dated
  strictly AFTER the tuning corpus (document year > CUTOFF_YEAR), so it could not have informed
  the criterion, its 3.0 A cutoff, or any decoy-selection choice even in principle.

RELEASE-GATED (the honest current state). Today the current CYP51/azole universe is exhausted:
the tuning corpus was pulled from the ChEMBL release current in 2026-08, so no document dated
> CUTOFF_YEAR (= 2026, i.e. >= 2027) exists yet. This builder therefore writes NOTHING and
reports NOT-YET-AVAILABLE until a newer ChEMBL release accrues post-cutoff depositions — the
same "ship the reproducible scaffold, defer the numeric run to the gated host" discipline the
metal-aware baseline harness (#83, scripts/rescore_baselines.py) uses. The rule, the cutoff,
the disjointness and the domain checks are all fixed now and unit-tested offline; only the
fetch waits on the newer release.

Discipline (mirrors data/ligands/PROVENANCE.md — only the provenance of "held out" changes,
to "deposited after the tuning corpus"):

  actives  measured-active azole CYP51 antifungals whose earliest antifungal / CYP51 activity
           document year > CUTOFF_YEAR, azole-bearing, in the frozen domain (<= 45 heavy, >= 1
           aromatic N), novel (max Morgan r=2 Tanimoto < 0.70 vs every prior active + the 5TZ1
           co-crystal), and disjoint by ChEMBL id from EVERY set the project has used.

  decoys   property-matched presumed-inactives, built by the SAME frozen rule as every other
           decoy set (openafr.inactives.select_inactives), themselves drawn only from post-
           cutoff documents and disjoint by id from every used set.

Deterministic: same newer-release ChEMBL snapshot -> byte-identical .smi (compounds emitted in
ascending ChEMBL-id order).

Usage:
    python scripts/make_temporal_holdout.py status              # report the release gate
    python scripts/make_temporal_holdout.py actives             # fetches ChEMBL (release-gated)
    python scripts/make_temporal_holdout.py decoys              # fetches ChEMBL (release-gated)
"""
import argparse
import hashlib
import json
import pathlib
import sys
import time
import urllib.parse
import urllib.request

from rdkit import Chem, DataStructs, RDLogger

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import inactives
from scripts.make_external_holdout import USED_ID_FILES, _id_of, _read_smi, _sha
from scripts.make_verified_actives import NOVELTY_MAX, OTESECONAZOLE, is_azole
from scripts.make_verified_inactives import ALL_ACTIVES

RDLogger.DisableLog("rdApp.*")

# The temporal boundary (FROZEN — see work/PREREGISTRATION_temporal.md).
# The tuning corpus was pulled from the ChEMBL release current in 2026-08; a document can carry
# year 2026, so a molecule that provably post-dates the whole corpus needs year > 2026. Strict
# ">" (not ">=") guarantees no overlap with the corpus's latest possible year.
CUTOFF_YEAR = 2026

ACTIVES_OUT = "data/ligands/temporal_actives.smi"
DECOYS_OUT = "data/ligands/temporal_decoys.smi"
SEED = 82                             # reuse #82's provenance seed; NOT 42 (the tuning seed)
PER_ACTIVE = 8

# Every set the criterion has ever seen, INCLUDING the #108 disjoint-population holdout — the
# temporal set must be disjoint from all of them by id (nothing already used can leak back in).
USED_FILES = list(USED_ID_FILES) + [
    "data/ligands/external_actives.smi", "data/ligands/external_decoys.smi",
]
# Whole-cell C. albicans + the C. albicans CYP51 enzyme target — the same activity universe the
# verified sets were built from (scripts/make_verified_inactives.QUERIES).
TARGET_QUERIES = (
    {"target_organism": "Candida albicans"},
    {"target_chembl_id": "CHEMBL1780"},
)
ACTIVITY = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
STATUS = "https://www.ebi.ac.uk/chembl/api/data/status.json"
PAGE = 1000


# --- pure, offline-testable temporal logic -------------------------------------------------

def earliest_document_year(records):
    """The earliest document year across a compound's activity records, or None if none carry
    a usable year. This is the compound's first *measured* antifungal appearance — the quantity
    the temporal split is defined on."""
    years = [r["document_year"] for r in records
             if isinstance(r.get("document_year"), int)]
    return min(years) if years else None


def is_temporal(earliest_year, cutoff=CUTOFF_YEAR):
    """A compound is in the temporal holdout iff its first measured antifungal activity was
    reported strictly AFTER the tuning corpus's release year — so it could not have informed
    the criterion. An unknown year (None) is NOT temporal: absence of proof is not proof."""
    return earliest_year is not None and earliest_year > cutoff


def used_ids(files=USED_FILES):
    """The union of ChEMBL ids across every set the criterion has used."""
    used = set()
    for f in files:
        if pathlib.Path(f).exists():
            used |= {_id_of(n) for _s, n in _read_smi(f)}
    return used


# --- network shell (release-gated; not exercised offline) ----------------------------------

def _get(url, retries=4):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except Exception as e:
            if attempt == retries - 1:
                raise SystemExit(f"ChEMBL fetch failed after {retries} tries: {e}\n  {url}")
            time.sleep(2 ** attempt)


def _fetch_activities(sel):
    """All (molecule_id, smiles, document_year) activity rows for one target selector."""
    fields = ("molecule_chembl_id,canonical_smiles,document_year,standard_type,"
              "standard_relation,standard_value,standard_units,activity_comment,"
              "pchembl_value,target_chembl_id")
    params = dict(sel, only=fields, limit=PAGE, offset=0, format="json")
    total = _get(f"{ACTIVITY}?{urllib.parse.urlencode(dict(params, limit=1))}")["page_meta"]["total_count"]
    rows, n = [], 0
    while n < total:
        page = _get(f"{ACTIVITY}?{urllib.parse.urlencode(dict(params, offset=n))}")
        got = page.get("activities", [])
        if not got:
            break
        rows += got
        n += len(got)
        print(f"    {n}/{total}", end="\r", flush=True)
    print()
    return rows


def _grouped_records():
    """{compound_id: [records]}, {compound_id: smiles} across the antifungal targets."""
    by_cmp, smiles = {}, {}
    for sel in TARGET_QUERIES:
        print(f"  fetch {sel} ...")
        for rec in _fetch_activities(sel):
            cid = rec.get("molecule_chembl_id")
            if not cid:
                continue
            by_cmp.setdefault(cid, []).append(rec)
            if rec.get("canonical_smiles"):
                smiles.setdefault(cid, rec["canonical_smiles"])
    return by_cmp, smiles


def _release():
    """(version, date) of the live ChEMBL release, or (None, None) if unreachable."""
    try:
        s = _get(STATUS)
        return s.get("chembl_db_version"), s.get("chembl_release_date")
    except SystemExit:
        return None, None


def cmd_status(args):
    ver, date = _release()
    print(f"temporal cutoff          : document year > {CUTOFF_YEAR}  (i.e. >= {CUTOFF_YEAR + 1})")
    print(f"live ChEMBL release      : {ver or 'unreachable'}  ({date or '-'})")
    print("gate                     : the holdout materializes only once a release carries")
    print(f"                           antifungal / CYP51 documents dated > {CUTOFF_YEAR}.")
    print("Run `actives` then `decoys`; each writes NOTHING and reports NOT-YET-AVAILABLE")
    print("until such documents exist (the current CYP51/azole universe is exhausted).")
    return 0


def _write_actives(by_cmp, smiles):
    ref_fps = []
    for path in ALL_ACTIVES:
        ref_fps += [a["fp"] for a in inactives.load_actives(path)]
    ref_fps.append(inactives.fingerprint(Chem.MolFromSmiles(OTESECONAZOLE)))
    already = used_ids()

    chosen, seen = [], set()
    for cid in sorted(by_cmp, key=lambda c: int(c[6:]) if c[6:].isdigit() else 10**12):
        if cid in already or cid not in smiles:
            continue
        verdict, _ = inactives.compound_verdict(by_cmp[cid])
        if verdict != "has-active-evidence":
            continue
        if not is_temporal(earliest_document_year(by_cmp[cid])):
            continue
        smi = smiles[cid]
        mol = Chem.MolFromSmiles(smi)
        if mol is None or smi in seen or not is_azole(mol):
            continue
        if not inactives.has_aromatic_nitrogen(mol) or mol.GetNumHeavyAtoms() > inactives.MAX_HEAVY:
            continue
        fp = inactives.fingerprint(mol)
        if max(DataStructs.BulkTanimotoSimilarity(fp, ref_fps), default=0) >= NOVELTY_MAX:
            continue
        seen.add(smi)
        chosen.append((cid, smi))
    return chosen


def cmd_actives(args):
    ver, date = _release()
    print(f"live ChEMBL release      : {ver or 'unreachable'}  ({date or '-'})")
    by_cmp, smiles = _grouped_records()
    print(f"antifungal compounds     : {len(by_cmp)}")
    chosen = _write_actives(by_cmp, smiles)
    if not chosen:
        print(f"\nNOT-YET-AVAILABLE: 0 measured-active azoles with a document year > {CUTOFF_YEAR}.")
        print("The temporal holdout is release-gated — the current CYP51/azole universe is")
        print("exhausted. Re-run against a newer ChEMBL release; nothing was written.")
        return 2
    with open(ACTIVES_OUT, "w") as fh:
        for cid, smi in chosen:
            fh.write(f"{smi}\tactive_{cid}\n")
    heavy = [Chem.MolFromSmiles(s).GetNumHeavyAtoms() for _c, s in chosen]
    print(f"\nwrote {len(chosen)} temporal actives (year > {CUTOFF_YEAR}) to {ACTIVES_OUT}")
    print(f"heavy atoms: min {min(heavy)}  mean {sum(heavy)/len(heavy):.1f}  max {max(heavy)}")
    print(f"sha256 {_sha(ACTIVES_OUT)}")
    return 0


def cmd_decoys(args):
    if not pathlib.Path(ACTIVES_OUT).exists():
        sys.exit(f"no {ACTIVES_OUT} — run `actives` first (release-gated).")
    actives = inactives.load_actives(ACTIVES_OUT)
    active_fps = [a["fp"] for a in actives]
    already = used_ids() | {_id_of(a["name"]) for a in actives}

    by_cmp, smiles = _grouped_records()
    # decoys must also be post-cutoff AND untouched by any used set
    cand = {}
    for cid in by_cmp:
        if cid in already or cid not in smiles:
            continue
        if not is_temporal(earliest_document_year(by_cmp[cid])):
            continue
        cand[cid] = smiles[cid]
    candidates = sorted(cand.items(), key=lambda kv: int(kv[0][6:]) if kv[0][6:].isdigit() else 10**12)
    print(f"post-cutoff decoy pool   : {len(candidates)} distinct, disjoint by id")

    chosen, rejects = inactives.select_inactives(candidates, actives, active_fps, PER_ACTIVE)
    print("\nfilter funnel:")
    for k, v in sorted(rejects.items(), key=lambda kv: -kv[1]):
        print(f"  rejected {k:20} {v:6d}")
    if not chosen:
        print(f"\nNOT-YET-AVAILABLE: no post-cutoff property-matched decoys. Nothing written.")
        return 2
    chosen.sort(key=lambda c: int(c["chembl_id"][6:]) if c["chembl_id"][6:].isdigit() else 10**12)
    with open(DECOYS_OUT, "w") as fh:
        for c in chosen:
            fh.write(f"{c['smiles']}\tdecoy_{c['chembl_id']}_for_{c['paired']}\n")
    print(f"\nwrote {len(chosen)} temporal decoys to {DECOYS_OUT}")
    print(f"sha256 {_sha(DECOYS_OUT)}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("status", help="report the release gate (needs network for the live release)"
                   ).set_defaults(func=cmd_status)
    sub.add_parser("actives", help="fetch + build the post-cutoff actives (release-gated)"
                   ).set_defaults(func=cmd_actives)
    sub.add_parser("decoys", help="fetch + build post-cutoff property-matched decoys (release-gated)"
                   ).set_defaults(func=cmd_decoys)
    args = ap.parse_args(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    sys.exit(main())
