#!/usr/bin/env python
"""Annotate a .smi library with synthesizability / availability (#86).

Per molecule, reports the Ertl SAscore (local, deterministic) and its band
(easy/moderate/hard), plus an availability annotation. The SAscore half needs no network.

Availability is deliberately conservative:
  --catalogued       mark every row 'catalogued-drug' (buyable) — correct for the
                     repurposing shortlist, drawn from the Broad Drug Repurposing Hub
                     (approved / clinical / catalogued tool compounds, all vendor-listed).
  --online           best-effort PubChem presence probe (a WEAK proxy for buyable; a found
                     record is not a vendor quote). Off by default so the run stays offline
                     and deterministic.
  (neither)          availability = 'not-checked'.

Informational, never a gate.

Usage:
    python scripts/synth_filter.py work/candidates/shortlist_focus.smi --approved-drugs
    python scripts/synth_filter.py lib.smi --online --json > out.json
"""
import argparse
import json
import pathlib
import sys
import urllib.request

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import synth
from scripts.admet_filter import read_smi

TSV_COLS = ["name", "sa_score", "sa_band", "availability", "pubchem_cid"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("smi", help="input .smi file (<smiles>\\t<name> per line)")
    ap.add_argument("--catalogued", action="store_true",
                    help="mark all rows buyable ('catalogued-drug') — true for the repurposing shortlist")
    ap.add_argument("--online", action="store_true",
                    help="best-effort PubChem presence probe (weak availability proxy)")
    ap.add_argument("--json", action="store_true", help="emit full profiles as JSON")
    ap.add_argument("--tsv", action="store_true", help="emit a flat TSV table")
    args = ap.parse_args(argv)

    opener = urllib.request.urlopen if args.online else None
    results = []
    for name, smiles in read_smi(args.smi):
        try:
            p = synth.profile(smiles, name=name,
                              catalogued=args.catalogued, opener=opener)
        except ValueError as e:
            sys.stderr.write("skip %s: %s\n" % (name, e))
            continue
        results.append(p)

    if args.json:
        json.dump(results, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 0

    if args.tsv:
        print("\t".join(TSV_COLS))
        for p in results:
            print("\t".join(str(p[c]) for c in TSV_COLS))
        return 0

    hdr = "%-16s %8s %-9s %-16s" % ("name", "SAscore", "band", "availability")
    print(hdr)
    print("-" * len(hdr))
    for p in results:
        print("%-16s %8.2f %-9s %-16s" % (
            (p["name"] or "")[:16], p["sa_score"], p["sa_band"], p["availability"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
