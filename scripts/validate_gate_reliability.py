"""Reliability validation gate — implements work/PREREGISTRATION_reliability.md.

Tests a pool-size-neutral aggregator on the same already-docked wild-type 5-seed poses the
consensus run used. Consensus ranked by the MINIMUM N-Fe over seeds and failed EF@1% because
`min` is an extreme-value statistic biased toward the larger class (348 decoys vs 7 actives).
This ranks by the MEDIAN of each molecule's per-seed best distances — how *reliably* it
reaches the iron — which a fixed number of draws makes pool-size-neutral.

Pass bar and grading math are identical to run 2 (from protocol.yaml): AUC >= 0.70 AND
EF@1% >= 5.0x on mode C. Unrankable actives (no coordinating pose in any seed) are ranked
LAST, never dropped. Consensus (min) and the seed-42 single search are reported alongside for
context but do NOT gate.

This is the SECOND and FINAL aggregator tested on these poses (see the pre-registration).

Usage:
    python scripts/validate_gate_reliability.py [CONSENSUS_DIR] [RECEPTOR]
        CONSENSUS_DIR  parent dir with per-seed subdirs s<seed>/  (default work/screen_consensus_wt)
        RECEPTOR       structure to read the heme iron from        (default work/receptor_A.pdb)
"""
import glob
import os
import pathlib
import statistics
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import protocol
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses
from openafr.scoring import report
from scripts.validate_gate2 import ACTIVES, check_prereg

PREREG = "work/PREREGISTRATION_reliability.md"
PREREG_SHA = "f8aa1ff93a8d77645e36af9a1b13f3f3d2a782001c0a6afa4b12a1bdd76ca4a7"
BASELINE_SEED = "42"  # the frozen protocol seed


def reliability_stats(byseed):
    """(median, minimum) of the per-seed best N-Fe distances.

    median = the pool-size-neutral reliability criterion (typical achievable coordination
    across independent searches); minimum = the consensus value, reported for context. Seeds
    with no coordinating pose (None) are ignored; if no seed coordinated, both are None so
    the gate ranks the molecule LAST rather than dropping it."""
    vals = [v for v in byseed.values() if v is not None]
    if not vals:
        return None, None
    return statistics.median(vals), min(vals)


def best_distance(pose_path, fe):
    """Minimum nitrogen-to-iron distance over all poses in one docked file, or None."""
    best = None
    for _num, atoms in read_poses(pose_path):
        d = min_nitrogen_iron_distance(atoms, fe)
        if d is not None and (best is None or d < best):
            best = d
    return best


def main():
    parent = sys.argv[1] if len(sys.argv) > 1 else "work/screen_consensus_wt"
    receptor = sys.argv[2] if len(sys.argv) > 2 else "work/receptor_A.pdb"
    fe = iron_position(receptor)
    actives = {l.split("\t")[1].strip() for l in open(ACTIVES) if "\t" in l}

    seed_dirs = sorted(d for d in glob.glob(os.path.join(parent, "s*"))
                       if os.path.isdir(d) and os.path.basename(d)[1:].isdigit())
    if not seed_dirs:
        sys.exit(f"no per-seed subdirs (s<seed>/) under {parent} — run scripts/screen_consensus.sh")
    seeds = [os.path.basename(d)[1:] for d in seed_dirs]

    gate = protocol.load()["gate"]
    min_auc, min_ef1 = gate["min_auc"], gate["min_ef1"]
    rep = lambda label, rows, gating=False: report(
        label, rows, gating, gate["ef_fractions"], gate["bedroc_alpha"])

    print("=" * 70)
    print("RELIABILITY GATE — PRE-REGISTERED WILD-TYPE HELD-OUT RE-VALIDATION")
    print("=" * 70)
    print(f"receptor         : {receptor}")
    print(f"seeds            : {seeds}  (median over per-seed minima; baseline = seed {BASELINE_SEED})")
    check_prereg(PREREG, PREREG_SHA)
    protocol.verify()

    per_mol = {}
    for d in seed_dirs:
        seed = os.path.basename(d)[1:]
        for f in sorted(os.listdir(d)):
            if f.endswith(".pdbqt"):
                per_mol.setdefault(f[:-6], {})[seed] = best_distance(os.path.join(d, f), fe)

    rel_rows, min_rows, base_rows, med = [], [], [], {}
    for name, byseed in per_mol.items():
        median, minimum = reliability_stats(byseed)
        med[name] = median
        is_act = name in actives
        rel_rows.append((name, median, is_act))
        min_rows.append((name, minimum, is_act))
        base_rows.append((name, byseed.get(BASELINE_SEED), is_act))

    seen = {r[0] for r in rel_rows if r[2]}
    if seen != actives:
        print(f"\n!! MISSING ACTIVES: {actives - seen} — results not valid")

    r = rep("MODE C (RELIABILITY) — median N-Fe over seeds  [PRIMARY — gates]", rel_rows, gating=True)
    m = rep("MODE C (CONSENSUS) — min N-Fe over seeds", min_rows)
    b = rep(f"MODE C (BASELINE, seed {BASELINE_SEED}) — single search", base_rows)

    print("\nheld-out actives by reliability (median), best first:")
    print("  %-16s %8s" % ("active", "median"))
    for a in sorted(actives, key=lambda a: med.get(a) or 9e9):
        print("  %-16s %8s" % (a, "-" if med.get(a) is None else "%.3f" % med[a]))

    print("\n" + "=" * 70)
    passed = r[0] >= min_auc and r[1][0] >= min_ef1
    print(f"RELIABILITY mode C: AUC {r[0]:.3f} (need >= {min_auc})   EF@1% {r[1][0]:.2f}x (need >= {min_ef1})")
    print(f"  (context) CONSENSUS AUC {m[0]:.3f} EF@1% {m[1][0]:.2f}x | BASELINE AUC {b[0]:.3f} EF@1% {b[1][0]:.2f}x")
    if passed:
        print("GATE: PASS — reliability criterion supported on the wild-type held-out set.")
        return 0
    print("GATE: FAIL — reliability criterion not supported. Per pre-registration, this is the")
    print("final aggregator on these poses; the honest conclusion is a geometry ceiling at the")
    print("top, not a fixable noise problem. No third aggregator will be tried.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
