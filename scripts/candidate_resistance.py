"""Per-candidate resistance-pocket retention (issues #69, #77).

Grades whether the shortlist candidates keep their iron-coordinating geometry in the
*C. auris* Y132F resistance pocket, relative to the wild-type *C. albicans* pocket where
they were ranked. This is a per-candidate CONFIDENCE AXIS, not a gate — it demotes/annotates,
never promotes. Criterion + bar are frozen in
work/PREREGISTRATION_candidate_resistance.md (sha work/PREREG_candidate_resistance.sha256)
BEFORE any Y132F distance was read.

Method (all frozen in the pre-registration):
  - Y132F consensus = the SAME consensus() logic as scripts/shortlist_confidence.py:
    MIN nitrogen-to-iron distance over seeds that coordinated (< FE_CUTOFF 3.0 A), with
    n_coord/n_seeds and spread, over the per-seed dock tree work/candidates/dock_auris_Y132F/s*.
  - WT baseline = the already-frozen `fungal_consensus` in
    work/candidates/shortlist_confidence.json (same 16 ligands, same 5 seeds, same frozen
    protocol, receptor work/receptor_A.pdb). NOT re-docked here.
  - retention_shift = y132f_consensus - wt_consensus  (+ = moved AWAY from the iron).

Pre-committed per-candidate classification:
  RETAINED  coordinates Y132F iron in >= MIN_COORD_SEEDS AND |retention_shift| <= RETENTION_TOL
  SHIFTED   coordinates in >= MIN_COORD_SEEDS but |retention_shift| > RETENTION_TOL (direction reported)
  LOST      coordinates in <  MIN_COORD_SEEDS on the Y132F pocket (reliably fails to reach the iron)

Retention != resistance-breaking: the held-out Y132F run showed the criterion is nearly
insensitive to Y132F, so RETAINED is a robustness check, not evidence a molecule defeats
resistance. LOST is the informative demotion this axis exists to surface.

Usage:
    python scripts/candidate_resistance.py \
        [--y132f work/candidates/dock_auris_Y132F] \
        [--y132f-receptor work/receptor_auris_Y132F.pdb] \
        [--wt-baseline work/candidates/shortlist_confidence.json] [--json]
"""
import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr.pdbqt import iron_position
# Reuse the frozen consensus logic unchanged — do not redefine it here.
from scripts.shortlist_confidence import FE_CUTOFF, consensus, per_seed_distances

RETENTION_TOL = 0.5     # A; |shift| beyond this = SHIFTED (frozen in the pre-registration)
MIN_COORD_SEEDS = 3     # coordinating seeds required to count as reaching the iron (frozen)


def classify(y132f_cons, y132f_ncoord, wt_cons):
    """Pre-committed RETAINED / SHIFTED / LOST verdict + the signed shift."""
    if y132f_ncoord < MIN_COORD_SEEDS:
        return "LOST", None
    if wt_cons is None or y132f_cons is None:
        return "UNRANKABLE", None
    shift = y132f_cons - wt_cons
    if abs(shift) <= RETENTION_TOL:
        return "RETAINED", shift
    return ("SHIFTED-away" if shift > 0 else "SHIFTED-closer"), shift


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--y132f", default="work/candidates/dock_auris_Y132F")
    ap.add_argument("--y132f-receptor", default="work/receptor_auris_Y132F.pdb")
    ap.add_argument("--wt-baseline", default="work/candidates/shortlist_confidence.json")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    fe = iron_position(args.y132f_receptor)
    y132f = per_seed_distances(args.y132f, fe)

    wt_conf = json.load(open(args.wt_baseline))
    wt = {r["name"]: r["fungal_consensus"] for r in wt_conf["candidates"]}

    recs = []
    for name in sorted(set(y132f) | set(wt)):
        y_cons, y_spread, y_ncoord, y_nseed = consensus(y132f.get(name, {}))
        wt_cons = wt.get(name)
        verdict, shift = classify(y_cons, y_ncoord, wt_cons)
        recs.append({
            "name": name,
            "wt_consensus": wt_cons,
            "y132f_consensus": None if y_cons is None else round(y_cons, 3),
            "y132f_seeds_coordinating": f"{y_ncoord}/{y_nseed}",
            "y132f_spread": round(y_spread, 3),
            "retention_shift": None if shift is None else round(shift, 3),
            "verdict": verdict,
            "per_seed_y132f": {k: (None if v is None else round(v, 3))
                               for k, v in sorted(y132f.get(name, {}).items())},
        })
    # Sort by WT consensus (the shortlist order); unrankable last.
    recs.sort(key=lambda r: (r["wt_consensus"] is None, r["wt_consensus"] or 9e9))

    if args.json:
        print(json.dumps({"fe_cutoff": FE_CUTOFF, "retention_tol": RETENTION_TOL,
                          "min_coord_seeds": MIN_COORD_SEEDS, "candidates": recs}, indent=2))
        return 0

    print("=" * 92)
    print("CANDIDATE RESISTANCE RETENTION — C. auris Y132F pocket (#69, #77)")
    print(f"receptor: {args.y132f_receptor}   coordination cutoff {FE_CUTOFF} A   "
          f"tol {RETENTION_TOL} A   min {MIN_COORD_SEEDS}/5 seeds")
    print("=" * 92)
    print("%-16s %8s %10s %8s %8s   %s" %
          ("candidate", "WT_min", "Y132F_min", "coord", "shift", "verdict"))
    print("-" * 92)
    for r in recs:
        wm = "-" if r["wt_consensus"] is None else "%.3f" % r["wt_consensus"]
        ym = "-" if r["y132f_consensus"] is None else "%.3f" % r["y132f_consensus"]
        sh = "-" if r["retention_shift"] is None else "%+.3f" % r["retention_shift"]
        print("%-16s %8s %10s %8s %8s   %s" %
              (r["name"], wm, ym, r["y132f_seeds_coordinating"], sh, r["verdict"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
