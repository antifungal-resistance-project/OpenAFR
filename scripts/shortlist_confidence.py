"""Per-candidate confidence: seed-stability (#68) + human-CYP51 selectivity (#73).

Reads two consensus dock trees (fungal + human), each a parent dir of per-seed subdirs
s<seed>/ produced by scripts/screen_consensus.sh / screen_human.sh, and for every ligand
computes:

  #68  Seed stability (fungal receptor)
    - per-seed best nitrogen-to-iron distance (min over that seed's poses)
    - consensus  = MIN over seeds that coordinated  (the pre-registered consensus criterion,
                   work/PREREGISTRATION_consensus.md — a less-biased estimate of the achievable
                   coordination)
    - spread     = max - min of the per-seed best over coordinating seeds (run-to-run noise)
    - n_coord    = how many of the seeds produced a coordinating pose (< FE_CUTOFF)
    - STABLE iff it coordinated in ALL seeds AND spread <= SPREAD_MAX. An unstable candidate's
      tight single-seed distance was not reproducible — exactly the failure the reliability
      result warned about (a single dock is a noisy estimator; clotrimazole once ranked
      #2774/2776). Unstable candidates are demoted, not trusted.

  #73  Human-CYP51 selectivity (fungal vs human receptor)
    - human_consensus = MIN over seeds on the HUMAN receptor (3LD6), same criterion
    - margin  = human_consensus - fungal_consensus  (Angstrom). POSITIVE => the nitrogen
                reaches the FUNGAL iron closer than the human iron = fungal-selective = good.
    - reaches_human = human_consensus < FE_CUTOFF  (does it coordinate the human iron at all?)
      A candidate that coordinates the fungal iron but NOT the human iron is cleanly selective;
      one that coordinates both is a selectivity concern to weigh against the real-azole
      reference drugs (which DO reach the human iron — that is their known toxicity).

FE_CUTOFF (3.0 A) is the established coordination cutoff (scripts/analyze_poses.py). The two
thresholds this script adds (SPREAD_MAX, and the sign convention on margin) are frozen in the
pre-registrations before any distance is read:
    work/PREREGISTRATION_shortlist_consensus.md   (#68)
    work/PREREGISTRATION_human_selectivity.md      (#73)

Usage:
    python scripts/shortlist_confidence.py \
        --fungal work/candidates/dock_fungal --human work/candidates/dock_human \
        [--fungal-receptor work/receptor_A.pdb] [--human-receptor work/receptor_human_A.pdb] \
        [--json]
"""
import argparse
import glob
import json
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses

FE_CUTOFF = 3.0    # coordination cutoff, from scripts/analyze_poses.py (unchanged)
SPREAD_MAX = 1.0   # A; max tolerated run-to-run spread for a "stable" candidate (pre-registered)


def best_distance(pose_path, fe):
    """Minimum nitrogen-to-iron distance over all poses in one docked file, or None."""
    best = None
    for _num, atoms in read_poses(pose_path):
        d = min_nitrogen_iron_distance(atoms, fe)
        if d is not None and (best is None or d < best):
            best = d
    return best


def per_seed_distances(parent, fe):
    """{ligand: {seed: best_N_Fe or None}} across all s<seed>/ subdirs under parent."""
    seed_dirs = sorted(d for d in glob.glob(os.path.join(parent, "s*"))
                       if os.path.isdir(d) and os.path.basename(d)[1:].isdigit())
    out = {}
    for d in seed_dirs:
        seed = os.path.basename(d)[1:]
        for f in sorted(os.listdir(d)):
            if f.endswith(".pdbqt"):
                out.setdefault(f[:-6], {})[seed] = best_distance(os.path.join(d, f), fe)
    return out


def consensus(byseed):
    """(consensus_min, spread, n_coord, n_seeds) from a {seed: best_N_Fe} dict.

    consensus = min over seeds that produced a coordinating pose (< FE_CUTOFF); spread =
    max-min over those coordinating seeds (0 if <2); n_coord = how many seeds coordinated."""
    coord = [v for v in byseed.values() if v is not None and v < FE_CUTOFF]
    n_seeds = len(byseed)
    if not coord:
        # No seed coordinated: fall back to the closest approach seen (still informative),
        # but n_coord = 0 marks it as not reliably reaching the iron.
        seen = [v for v in byseed.values() if v is not None]
        return (min(seen) if seen else None, 0.0, 0, n_seeds)
    return (min(coord), (max(coord) - min(coord) if len(coord) > 1 else 0.0), len(coord), n_seeds)


def assess(name, fungal_byseed, human_byseed):
    """One ligand's full confidence record (pure)."""
    f_cons, f_spread, f_ncoord, f_nseed = consensus(fungal_byseed)
    h_cons, h_spread, h_ncoord, h_nseed = consensus(human_byseed)
    stable = (f_ncoord == f_nseed and f_nseed > 0 and f_spread <= SPREAD_MAX)
    reaches_fungal = f_ncoord > 0
    reaches_human = h_ncoord > 0
    margin = (h_cons - f_cons) if (f_cons is not None and h_cons is not None) else None
    return {
        "name": name,
        "fungal_consensus": None if f_cons is None else round(f_cons, 3),
        "fungal_spread": round(f_spread, 3),
        "fungal_seeds_coordinating": f"{f_ncoord}/{f_nseed}",
        "reaches_fungal_iron": reaches_fungal,
        "stable": stable,
        "human_consensus": None if h_cons is None else round(h_cons, 3),
        "human_seeds_coordinating": f"{h_ncoord}/{h_nseed}",
        "reaches_human_iron": reaches_human,
        "selectivity_margin": None if margin is None else round(margin, 3),
        "per_seed_fungal": {k: (None if v is None else round(v, 3))
                            for k, v in sorted(fungal_byseed.items())},
        "per_seed_human": {k: (None if v is None else round(v, 3))
                           for k, v in sorted(human_byseed.items())},
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fungal", default="work/candidates/dock_fungal")
    ap.add_argument("--human", default="work/candidates/dock_human")
    ap.add_argument("--fungal-receptor", default="work/receptor_A.pdb")
    ap.add_argument("--human-receptor", default="work/receptor_human_A.pdb")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    f_fe = iron_position(args.fungal_receptor)
    h_fe = iron_position(args.human_receptor)
    fungal = per_seed_distances(args.fungal, f_fe)
    human = per_seed_distances(args.human, h_fe)
    names = sorted(set(fungal) | set(human))

    recs = [assess(n, fungal.get(n, {}), human.get(n, {})) for n in names]
    # Sort by fungal consensus (tightest first); unrankable (None) last.
    recs.sort(key=lambda r: (r["fungal_consensus"] is None, r["fungal_consensus"] or 9e9))

    if args.json:
        print(json.dumps({"fe_cutoff": FE_CUTOFF, "spread_max": SPREAD_MAX,
                          "candidates": recs}, indent=2))
        return 0

    print("=" * 100)
    print("SHORTLIST CONFIDENCE — seed stability (#68) + human-CYP51 selectivity (#73)")
    print(f"fungal receptor: {args.fungal_receptor}   human receptor: {args.human_receptor}")
    print(f"coordination cutoff {FE_CUTOFF} A   stable = coordinates every seed AND "
          f"spread <= {SPREAD_MAX} A")
    print("=" * 100)
    hdr = ("candidate", "fung_min", "spread", "f_coord", "stable",
           "hum_min", "h_coord", "margin", "sel")
    print("%-16s %8s %7s %8s %7s %8s %8s %8s  %s" % hdr)
    print("-" * 100)
    for r in recs:
        fm = "-" if r["fungal_consensus"] is None else "%.3f" % r["fungal_consensus"]
        hm = "-" if r["human_consensus"] is None else "%.3f" % r["human_consensus"]
        mg = "-" if r["selectivity_margin"] is None else "%+.3f" % r["selectivity_margin"]
        # selectivity verdict: fungal-selective if reaches fungal iron and (not human, or margin>0)
        if not r["reaches_fungal_iron"]:
            sel = "no-fungal"
        elif not r["reaches_human_iron"]:
            sel = "SELECTIVE (fungal only)"
        elif r["selectivity_margin"] is not None and r["selectivity_margin"] > 0:
            sel = "fungal-closer"
        else:
            sel = "hits-human>=fungal"
        print("%-16s %8s %7s %8s %7s %8s %8s %8s  %s" % (
            r["name"], fm, "%.3f" % r["fungal_spread"], r["fungal_seeds_coordinating"],
            "yes" if r["stable"] else "NO", hm, r["human_seeds_coordinating"], mg, sel))
    return 0


if __name__ == "__main__":
    sys.exit(main())
