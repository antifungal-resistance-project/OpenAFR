#!/usr/bin/env python
"""Profile a .smi library for ADMET / drug-likeness + structural alerts (#75).

Reads a whitespace-separated SMILES file ("<smiles>\\t<name>" per line, the format
build_shortlist.py and the shortlist_*.smi files already use) and, per molecule, prints
the physicochemical profile plus any PAINS/BRENK/NIH structural-alert and hERG-risk flags
from openafr.admet. Everything is a LOCAL RDKit computation — no network.

These are informational annotations on a geometry rank, never a gate: a marketed azole
can (and ketoconazole does) violate a rule and still be a drug. Read the flags, don't
auto-reject on them.

Usage:
    python scripts/admet_filter.py work/candidates/shortlist_focus.smi
    python scripts/admet_filter.py work/candidates/shortlist_focus.smi --json > out.json
    python scripts/admet_filter.py work/candidates/shortlist_focus.smi --tsv  > out.tsv
"""
import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import admet


def read_smi(path):
    """Yield (name, smiles) from a .smi file. A trailing '|...|' CXSMILES block is kept."""
    rows = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split()
            if len(parts) < 2:
                continue
            smiles, name = parts[0], parts[-1]
            rows.append((name, smiles))
    return rows


def _flag_summary(p):
    """One-glance flag string: what a reader should go check."""
    flags = []
    if p["ro5_violations"] >= 2:
        flags.append("Ro5x%d" % p["ro5_violations"])
    if not p["veber_ok"]:
        flags.append("Veber")
    if p["esol_logs"] <= -6.0:
        flags.append("insoluble")
    if p["pains"]:
        flags.append("PAINS")
    if p["structural_alerts"]:
        flags.append("alert")
    if p["herg_risk_factors"]:
        flags.append("hERG?")
    return ",".join(flags) if flags else "clean"


TSV_COLS = [
    "name", "mw", "clogp", "hbd", "hba", "tpsa", "rotatable_bonds",
    "aromatic_rings", "esol_logs", "ro5_violations", "veber_ok",
    "n_pains", "n_alerts", "herg_risk_factors", "flags",
]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("smi", help="input .smi file (<smiles>\\t<name> per line)")
    ap.add_argument("--json", action="store_true", help="emit full profiles as JSON")
    ap.add_argument("--tsv", action="store_true", help="emit a flat TSV table")
    args = ap.parse_args(argv)

    rows = read_smi(args.smi)
    results = []
    for name, smiles in rows:
        try:
            p = admet.profile(smiles)
        except ValueError as e:
            sys.stderr.write("skip %s: %s\n" % (name, e))
            continue
        p["name"] = name
        p["flags"] = _flag_summary(p)
        results.append(p)

    if args.json:
        json.dump(results, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 0

    if args.tsv:
        print("\t".join(TSV_COLS))
        for p in results:
            print("\t".join(str(x) for x in (
                p["name"], p["mw"], p["clogp"], p["hbd"], p["hba"], p["tpsa"],
                p["rotatable_bonds"], p["aromatic_rings"], p["esol_logs"],
                p["ro5_violations"], p["veber_ok"], len(p["pains"]),
                len(p["structural_alerts"]), p["herg_risk_factors"], p["flags"],
            )))
        return 0

    # default: human-readable table
    hdr = "%-16s %6s %6s %5s %5s %6s %5s %6s %5s %-24s" % (
        "name", "MW", "cLogP", "HBD", "HBA", "TPSA", "RotB", "logS", "Ro5v", "flags")
    print(hdr)
    print("-" * len(hdr))
    for p in results:
        print("%-16s %6.1f %6.2f %5d %5d %6.1f %5d %6.2f %5d %-24s" % (
            p["name"][:16], p["mw"], p["clogp"], p["hbd"], p["hba"], p["tpsa"],
            p["rotatable_bonds"], p["esol_logs"], p["ro5_violations"], p["flags"]))
        for d in p["pains"]:
            print("    PAINS: %s" % d)
        for fam, d in p["structural_alerts"]:
            print("    %s: %s" % (fam, d))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
