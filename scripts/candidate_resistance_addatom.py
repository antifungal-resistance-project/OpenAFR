"""Per-candidate resistance-pocket retention for the ADD-ATOM C. auris mutants (issue #77).

The sibling of scripts/candidate_resistance.py (the Y132F run). Same criterion, same bar,
same interpretation — this grades the two clade-specific substitutions that ADD side-chain
atoms and therefore need the pinned repacker (scripts/repack_mutant.py):

    K143R   clade I     Lys -> Arg      work/receptor_auris_K143R.pdbqt
    F126L   clade III   Phe -> Leu      work/receptor_auris_F126L.pdbqt

Both receptors carry a MODELED (minimum-strain) rotamer, not an observed geometry — a strictly
weaker structural claim than the Y132F deletion (see
work/PREREGISTRATION_candidate_resistance_addatom.md). The verdict rule itself is unchanged:
the frozen `classify` from candidate_resistance.py, pinned by
work/PREREGISTRATION_candidate_resistance.md and its tests, is reused verbatim — nothing about
the pass bar is re-tuned for these mutants.

Method (frozen in the add-atom pre-registration, inheriting the Y132F run):
  - mutant consensus = the SAME consensus() logic as scripts/shortlist_confidence.py: MIN
    nitrogen-to-iron distance over seeds that coordinated (< FE_CUTOFF 3.0 A), with
    n_coord/n_seeds and spread, over the per-seed dock tree work/candidates/dock_auris_<M>/s*.
  - WT baseline = the already-frozen `fungal_consensus` in
    work/candidates/shortlist_confidence.json (NOT re-docked).
  - retention_shift = mutant_consensus - wt_consensus  (+ = moved AWAY from the iron).

Usage:
    python scripts/candidate_resistance_addatom.py K143R \
        --dock work/candidates/dock_auris_K143R \
        --receptor work/receptor_auris_K143R.pdbqt \
        --out work/candidates/shortlist_resistance_K143R.json
"""
import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr.pdbqt import iron_position
# Reuse the frozen verdict rule + consensus logic unchanged — do not redefine them here.
from scripts.candidate_resistance import MIN_COORD_SEEDS, RETENTION_TOL, classify
from scripts.shortlist_confidence import FE_CUTOFF, consensus, per_seed_distances


def grade(dock_dir, receptor_pdb, wt_baseline, mutant):
    """Per-candidate retention records for one add-atom mutant (pure w.r.t. the dock tree).

    `receptor_pdb` is the heavy-atom .pdb (HETATM heme), the same convention the Y132F and
    ensemble graders use to read the iron — not the obabel .pdbqt (which relabels the FE
    record). The iron position is identical either way (the mutation never moves it)."""
    fe = iron_position(receptor_pdb)
    byname = per_seed_distances(dock_dir, fe)
    wt = {r["name"]: r["fungal_consensus"] for r in json.load(open(wt_baseline))["candidates"]}

    recs = []
    for name in sorted(set(byname) | set(wt)):
        m_cons, m_spread, m_ncoord, m_nseed = consensus(byname.get(name, {}))
        wt_cons = wt.get(name)
        verdict, shift = classify(m_cons, m_ncoord, wt_cons)
        recs.append({
            "name": name,
            "wt_consensus": wt_cons,
            "mutant_consensus": None if m_cons is None else round(m_cons, 3),
            "mutant_seeds_coordinating": f"{m_ncoord}/{m_nseed}",
            "mutant_spread": round(m_spread, 3),
            "retention_shift": None if shift is None else round(shift, 3),
            "verdict": verdict,
            "per_seed": {k: (None if v is None else round(v, 3))
                        for k, v in sorted(byname.get(name, {}).items())},
        })
    # Sort by WT consensus (the shortlist order); unrankable last.
    recs.sort(key=lambda r: (r["wt_consensus"] is None, r["wt_consensus"] or 9e9))
    return {
        "mutant": mutant,
        "clade": {"K143R": "I", "F126L": "III"}.get(mutant),
        "modeled_rotamer": True,
        "fe_cutoff": FE_CUTOFF,
        "retention_tol": RETENTION_TOL,
        "min_coord_seeds": MIN_COORD_SEEDS,
        "candidates": recs,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mutant", choices=["K143R", "F126L"])
    ap.add_argument("--dock")
    ap.add_argument("--receptor", help="heavy-atom .pdb to read the iron from "
                    "(default: work/receptor_auris_<mutant>.pdb)")
    ap.add_argument("--wt-baseline", default="work/candidates/shortlist_confidence.json")
    ap.add_argument("--out", help="write JSON here (default: stdout)")
    args = ap.parse_args(argv)

    dock = args.dock or f"work/candidates/dock_auris_{args.mutant}"
    receptor = args.receptor or f"work/receptor_auris_{args.mutant}.pdb"
    out = grade(dock, receptor, args.wt_baseline, args.mutant)

    text = json.dumps(out, indent=2)
    if args.out:
        pathlib.Path(args.out).write_text(text + "\n")
        n_lost = sum(1 for r in out["candidates"] if r["verdict"] == "LOST")
        print(f"{args.mutant}: wrote {args.out}  "
              f"({len(out['candidates'])} candidates, {n_lost} LOST)")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
