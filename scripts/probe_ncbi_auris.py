"""Spike probe: is NCBI Pathogen Detection a usable feed for a C. auris azole early-warning?

Pipeline stage: resistance early-warning / source validation (issue #20)

What this answers
-----------------
The early-warning track (issues #20-27) rests on one premise: that a public feed
delivers *fresh, dated, geolocated* C. auris isolates whose azole-resistance state
we can read. Before building anything downstream, this probe interrogates the feed
and reports the four go/no-go metrics from issue #20:

  1. collection dates present at a usable fill rate
  2. geography present at a usable fill rate
  3. azole-resistance calls (ERG11/CYP51 mutations, or AST phenotype) EXPOSED by
     NCBI, or at least DERIVABLE because raw reads / assemblies are linked
  4. new isolates arrive at a weeks-to-months cadence, not years

Where the data lives
--------------------
NCBI Pathogen Detection publishes per-organism "PDG" result sets over plain HTTP:

  https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Candidozyma_auris/

(C. auris now sits under its reassigned genus *Candidozyma*.) Each release folder
PDG000000067.<N>/ carries Metadata/<...>.metadata.tsv (one row per isolate) and
AMR/<...>.amr.metadata.tsv. This probe pulls the metadata TSV for the latest
release and computes the metrics. Stdlib only (urllib + csv) on purpose: the repo
keeps its dependency surface to rdkit, and a source probe should not need pandas.

Run
---
  python scripts/probe_ncbi_auris.py                 # latest release, live
  python scripts/probe_ncbi_auris.py --release 670   # pin a release
  python scripts/probe_ncbi_auris.py --tsv local.tsv # offline, from a saved TSV
"""
import argparse
import collections
import csv
import io
import re
import sys
import urllib.request

BASE = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Candidozyma_auris/"
ORG_PREFIX = "PDG000000067"

# Values that mean "no real data" in BioSample-derived metadata.
NULLISH = {"", "null", "missing", "not provided", "not collected",
           "not applicable", "n/a", "na", "unknown", "not determined"}

# Azole-resistance ERG11/CYP51 hot-spot substitutions we would need to see called.
AZOLE_HOTSPOTS = ("Y132F", "K143R", "F126L", "F126T", "VF125AL", "Y132")


def _get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "openafr-earlywarning-probe/1"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        return resp.read()


def latest_release():
    html = _get(BASE).decode("utf-8", "replace")
    versions = re.findall(rf'href="{ORG_PREFIX}\.(\d+)/"', html)
    if not versions:
        sys.exit("could not list releases at " + BASE)
    return max(int(v) for v in versions)


def load_rows(release=None, tsv_path=None):
    if tsv_path:
        with open(tsv_path, newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t")), "local:" + tsv_path
    if release is None:
        release = latest_release()
    tag = f"{ORG_PREFIX}.{release}"
    url = f"{BASE}{tag}/Metadata/{tag}.metadata.tsv"
    text = _get(url).decode("utf-8", "replace")
    return list(csv.DictReader(io.StringIO(text), delimiter="\t")), tag


def filled(rows, col):
    return sum(1 for r in rows
               if (r.get(col) or "").strip().lower() not in NULLISH)


def month_hist(rows, col):
    h = collections.Counter()
    for r in rows:
        v = (r.get(col) or "").strip()
        if len(v) >= 7 and v[:4].isdigit():
            h[v[:7]] += 1
    return h


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", type=int, help="pin a PDG release number")
    ap.add_argument("--tsv", help="read a locally saved metadata.tsv instead of fetching")
    ap.add_argument("--recent-months", type=int, default=12)
    args = ap.parse_args()

    rows, tag = load_rows(args.release, args.tsv)
    n = len(rows)
    if not n:
        sys.exit("no isolate rows found")

    def pct(col):
        c = filled(rows, col)
        return c, 100.0 * c / n

    print(f"# NCBI Pathogen Detection probe — Candidozyma (Candida) auris")
    print(f"release: {tag}   isolates: {n}\n")

    print("[1] dates")
    for col in ("collection_date", "target_creation_date", "sra_release_date"):
        c, p = pct(col)
        print(f"    {col:22s} {c:6d}/{n} = {p:5.1f}%")

    print("[2] geography")
    for col in ("geo_loc_name", "lat_lon"):
        c, p = pct(col)
        print(f"    {col:22s} {c:6d}/{n} = {p:5.1f}%")

    print("[3] azole-resistance calls")
    for col in ("AMR_genotypes", "AST_phenotypes", "computed_types"):
        c, p = pct(col)
        print(f"    {col:22s} {c:6d}/{n} = {p:5.1f}%")
    applied = filled(rows, "amrfinder_applied") and sum(
        1 for r in rows if (r.get("amrfinder_applied") or "0") not in ("0", "", "NULL"))
    hotspot_hits = sum(1 for r in rows
                       if any(h in (r.get("AMR_genotypes") or "") for h in AZOLE_HOTSPOTS))
    print(f"    amrfinder actually applied to: {applied} isolates")
    print(f"    rows mentioning an ERG11 azole hot-spot: {hotspot_hits}")
    #   derivability: can WE call the mutations even if NCBI does not?
    for col in ("Run", "asm_acc"):
        c, p = pct(col)
        label = "SRA raw reads" if col == "Run" else "assembled genome"
        print(f"    derivable via {label:16s} ({col}): {c:6d}/{n} = {p:5.1f}%")

    print("[4] cadence (new DB targets per month, target_creation_date)")
    hist = month_hist(rows, "target_creation_date")
    for k in sorted(hist)[-args.recent_months:]:
        print(f"    {k}: {hist[k]}")

    print("\ntop geographies:")
    geo = collections.Counter((r.get("geo_loc_name") or "NULL").split(":")[0].strip()
                              for r in rows)
    for name, c in geo.most_common(8):
        print(f"    {name:20s} {c}")


if __name__ == "__main__":
    main()
