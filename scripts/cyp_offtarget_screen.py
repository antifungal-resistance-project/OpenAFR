#!/usr/bin/env python
"""Counter-screen a .smi library for human-CYP off-target liability (#74).

Reads a whitespace/tab-separated SMILES file ("<smiles>\\t<name>" per line — the format
build_shortlist.py and the shortlist_*.smi files use) and, per molecule, reports the
per-isoform inhibition-liability profile from openafr.cyp_offtarget for the five human
CYPs that dominate drug metabolism (3A4, 2C9, 2D6, 2C19, 1A2). Everything is a LOCAL
RDKit computation — no network, no docking.

The headline signal is mechanistic: an azole warhead (imidazole/triazole/tetrazole) is a
strong type-II heme ligand and inhibits every human CYP the same way it inhibits fungal
CYP51 — the very feature the geometry rank selected for. These are informational
annotations on a rank, never a gate. Read the flags; for an approved drug the label is
the authority.

Usage:
    python scripts/cyp_offtarget_screen.py work/candidates/shortlist_focus.smi
    python scripts/cyp_offtarget_screen.py work/candidates/shortlist_focus.smi --json > out.json
    python scripts/cyp_offtarget_screen.py work/candidates/shortlist_focus.smi --tsv  > out.tsv
"""
import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import cyp_offtarget
from scripts.admet_filter import read_smi

PANEL = cyp_offtarget.PANEL
TSV_COLS = ["name", "mw", "clogp", "azole_warhead", "type_ii_ligator", "n_high",
            "n_moderate"] + list(PANEL) + ["flags"]


def _record(name, smiles):
    p = cyp_offtarget.profile(smiles)
    p["name"] = name
    p["flags"] = cyp_offtarget.flag_summary(p)
    return p


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("smi", help="input .smi file (<smiles>\\t<name> per line)")
    ap.add_argument("--json", action="store_true", help="emit full profiles as JSON")
    ap.add_argument("--tsv", action="store_true", help="emit a flat TSV table")
    args = ap.parse_args(argv)

    results = []
    for name, smiles in read_smi(args.smi):
        try:
            results.append(_record(name, smiles))
        except ValueError as e:
            sys.stderr.write("skip %s: %s\n" % (name, e))

    if args.json:
        json.dump({"panel": list(PANEL), "candidates": results}, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 0

    if args.tsv:
        print("\t".join(TSV_COLS))
        for p in results:
            print("\t".join(str(x) for x in (
                [p["name"], p["mw"], p["clogp"], p["azole_warhead"] or "-",
                 p["type_ii_ligator"], p["n_high"], p["n_moderate"]]
                + [p["panel"][c]["liability"] for c in PANEL]
                + [p["flags"]])))
        return 0

    # default: human-readable table
    hdr = "%-14s %6s %6s  %-8s %-8s %-8s %-8s %-8s  %s" % (
        ("name", "MW", "cLogP") + PANEL + ("flags",))
    print(hdr)
    print("-" * len(hdr))
    for p in results:
        cells = tuple(p["panel"][c]["liability"] for c in PANEL)
        print("%-14s %6.1f %6.2f  %-8s %-8s %-8s %-8s %-8s  %s" % (
            (p["name"][:14], p["mw"], p["clogp"]) + cells + (p["flags"],)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
