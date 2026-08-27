#!/usr/bin/env python
"""Sensitivity of the applicability-domain thresholds (#87).

The shortlist is only valid "within the declared applicability domain (<=45 heavy atoms,
>=1 aromatic N)". This asks whether those two lines are principled — set by where the known
actives live and by docking reliability — or arbitrary cutoffs that happen to include our
hits. It is a fully LOCAL, deterministic RDKit reanalysis (no docking, no network) over
molecules whose SMILES persist in the repo.

Three questions, three answers:

  A. Envelope vs. the labelled set.  Sweep each threshold and report how many of the 15
     validated actives and the 399 property-matched decoys stay in-domain. Shows whether
     the chosen line tracks the actives' own envelope rather than the shortlist.

  B. Shortlist membership.  Over the same shortlist, sweep the thresholds and report which
     candidates enter/leave — i.e. how fragile the shortlist is to moving a line.

  C. Per-candidate fragility flags.  For each shortlist candidate, its heavy-atom count and
     aromatic-nitrogen count, flagged:
       size-fragile         heavy atoms within MARGIN of the cap (a small tightening drops it)
       coordinator-fragile  exactly one aromatic N (a >=2 rule would drop it; its whole
                            in-domain status rests on a single potential iron coordinator)

The two frozen thresholds are imported from scripts.filter_library (single source of truth),
never re-hardcoded here.

Usage:
    conda run -n openafr python scripts/domain_sensitivity.py
    conda run -n openafr python scripts/domain_sensitivity.py --json
"""
import argparse
import json
import os
import sys

from rdkit import Chem, RDLogger

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from scripts.filter_library import AROMATIC_N, MAX_HEAVY_ATOMS

RDLogger.DisableLog("rdApp.*")
LIG = os.path.join(ROOT, "data", "ligands")
ACTIVE_FILES = [os.path.join(LIG, "actives.smi"),
                os.path.join(LIG, "actives_holdout_final.smi")]
DECOY_FILE = os.path.join(LIG, "decoys.smi")
SHORTLIST = os.path.join(ROOT, "work", "candidates", "shortlist_focus.smi")

SIZE_MARGIN = 5   # heavy atoms within this many of the cap -> size-fragile


def descriptors(path):
    """[(name, heavy_atoms, n_aromatic_N)] for a SMILES<TAB>name file."""
    out = []
    for line in open(path):
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        smi, _, name = line.partition("\t")
        mol = Chem.MolFromSmiles(smi.strip())
        if mol is None:
            continue
        aromatic_n = sum(1 for a in mol.GetAtoms()
                         if a.GetIsAromatic() and a.GetSymbol() == "N")
        out.append((name.strip() or "(unnamed)", mol.GetNumHeavyAtoms(), aromatic_n))
    return out


def in_domain(heavy, aromatic_n, cap, min_arom_n):
    return heavy <= cap and aromatic_n >= min_arom_n


def sweep(mols, caps, min_arom_ns):
    """For each (cap, min_arom_n), fraction of `mols` retained in-domain."""
    grid = {}
    n = len(mols)
    for cap in caps:
        for m in min_arom_ns:
            k = "cap=%d,aromN>=%d" % (cap, m)
            grid[k] = sum(1 for _, h, a in mols if in_domain(h, a, cap, m))
    return grid, n


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", action="store_true", help="emit the full analysis as JSON")
    args = ap.parse_args(argv)

    actives = descriptors(ACTIVE_FILES[0]) + descriptors(ACTIVE_FILES[1])
    decoys = descriptors(DECOY_FILE)
    shortlist = descriptors(SHORTLIST)

    caps = [30, 35, 36, 40, 45, 48, 50, 55]
    aromNs = [0, 1, 2]

    act_grid, n_act = sweep(actives, caps, aromNs)
    dec_grid, n_dec = sweep(decoys, caps, aromNs)
    short_grid, n_short = sweep(shortlist, caps, aromNs)

    # per-candidate fragility
    frag = []
    for name, heavy, aromatic_n in shortlist:
        flags = []
        if heavy > (MAX_HEAVY_ATOMS - SIZE_MARGIN):
            flags.append("size-fragile")
        if aromatic_n == 1:
            flags.append("coordinator-fragile")
        frag.append({"name": name, "heavy": heavy, "aromatic_n": aromatic_n,
                     "flags": flags})

    # the gap around the cap in the validated actives (why 45 is robust)
    act_heavy = sorted(h for _, h, _ in actives)
    below_cap = [h for h in act_heavy if h <= MAX_HEAVY_ATOMS]
    above_cap = [h for h in act_heavy if h > MAX_HEAVY_ATOMS]
    gap = {"largest_in_domain_active_heavy": max(below_cap) if below_cap else None,
           "smallest_out_of_domain_active_heavy": min(above_cap) if above_cap else None,
           "largest_shortlist_heavy": max(h for _, h, _ in shortlist),
           "actives_min_aromatic_n": min(a for _, _, a in actives),
           "decoys_with_zero_aromatic_n": sum(1 for _, _, a in decoys if a == 0)}

    result = {
        "frozen_thresholds": {"max_heavy_atoms": MAX_HEAVY_ATOMS, "min_aromatic_n": 1},
        "n": {"actives": n_act, "decoys": n_dec, "shortlist": n_short},
        "sweep_actives_retained": act_grid,
        "sweep_decoys_retained": dec_grid,
        "sweep_shortlist_retained": short_grid,
        "candidate_fragility": frag,
        "envelope_gap": gap,
    }

    if args.json:
        json.dump(result, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 0

    print("Applicability-domain sensitivity (frozen: <=%d heavy atoms, >=1 aromatic N)"
          % MAX_HEAVY_ATOMS)
    print("n: %d validated actives, %d decoys, %d shortlist\n"
          % (n_act, n_dec, n_short))

    print("A. Heavy-atom cap sweep (aromatic-N held at >=1) — retained / total")
    print("   %-8s %-14s %-12s %-12s" % ("cap", "actives", "decoys", "shortlist"))
    for cap in caps:
        k = "cap=%d,aromN>=1" % cap
        print("   %-8d %-14s %-12s %-12s" % (
            cap, "%d/%d" % (act_grid[k], n_act),
            "%d/%d" % (dec_grid[k], n_dec),
            "%d/%d" % (short_grid[k], n_short)))

    print("\nB. Aromatic-N floor sweep (heavy cap held at %d) — retained / total"
          % MAX_HEAVY_ATOMS)
    print("   %-10s %-14s %-12s %-12s" % ("aromN>=", "actives", "decoys", "shortlist"))
    for m in aromNs:
        k = "cap=%d,aromN>=%d" % (MAX_HEAVY_ATOMS, m)
        print("   %-10d %-14s %-12s %-12s" % (
            m, "%d/%d" % (act_grid[k], n_act),
            "%d/%d" % (dec_grid[k], n_dec),
            "%d/%d" % (short_grid[k], n_short)))

    print("\nC. Per-candidate fragility (heavy, aromatic-N)")
    for c in sorted(frag, key=lambda d: (-len(d["flags"]), d["name"])):
        tag = ", ".join(c["flags"]) if c["flags"] else "robust"
        print("   %-16s heavy=%2d aromN=%d   %s"
              % (c["name"][:16], c["heavy"], c["aromatic_n"], tag))

    g = gap
    print("\nWhy the cap is robust: largest in-domain active = %s heavy, smallest "
          "out-of-domain active = %s heavy (a gap with no active in it); largest "
          "shortlist candidate = %s heavy."
          % (g["largest_in_domain_active_heavy"],
             g["smallest_out_of_domain_active_heavy"], g["largest_shortlist_heavy"]))
    print("Why aromatic-N >=1 is a scope gate, not a discriminator: %d/%d decoys already "
          "have >=1 aromatic N, and every validated active has >=%d."
          % (n_dec - g["decoys_with_zero_aromatic_n"], n_dec, g["actives_min_aromatic_n"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
