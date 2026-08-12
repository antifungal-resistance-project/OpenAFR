"""Axial-geometry validation gate — implements work/PREREGISTRATION_axial.md.

The FIRST feature tested on the run-2 wild-type poses that is *orthogonal* to N-Fe
distance. Every prior look (run-2 baseline, consensus, reliability) ranked on the same one
number — how CLOSE the coordinating nitrogen gets to the heme iron — and the reliability
pre-registration closed the "aggregate that distance differently" line for good. Azole
inhibition is *axial*: the drug nitrogen approaches the iron along the porphyrin normal, not
from the side. Distance is blind to direction; a property-matched decoy can reach a short
N-Fe distance with an off-axis nitrogen doing no real axial coordination. This gate adds the
angle.

Per pose, using the SAME nitrogen the distance criterion uses (the one nearest the iron):
    v = N - Fe                      |v| = the N-Fe distance
    theta = acute angle between v and the heme normal n, folded to [0, 90] deg
    S = |v| * (1 + sin theta)       mode E — rewards close AND on-axis (lambda = 1, untuned)

Heme normal n is fixed from the receptor once: from HEM A 601's four pyrrole nitrogens,
n = normalize((NC - NA) x (ND - NB)). Because theta is folded to the acute angle and the
feature uses sin theta, the orientation (sign) of n is irrelevant — nothing to tune.

Aggregation carries the reliability lesson forward: per molecule take the per-seed MINIMUM S
(R = 5 seeds), then rank by the MEDIAN of those 5 values (smaller = better). Median, not min,
because consensus proved min-over-R is an extreme-value statistic biased toward the larger
class (348 decoys vs 7 actives).

Pass bar and grading math are identical to run 2 (from protocol.yaml): AUC >= 0.70 AND
EF@1% >= 5.0x on mode E. Unrankable actives (no nitrogen in any seed) are ranked LAST, never
dropped. Reported alongside but NOT gating: pure axial angle (median per-seed min theta),
reliability distance (median per-seed min N-Fe), and the seed-42 single-search S.

This is the 4th look at these poses and the 1st orthogonal to distance. Per the
pre-registration, it is the first AND ONLY orthogonal feature tried on this dataset:
whatever the outcome, no second feature is fished on these poses.

Usage:
    python scripts/validate_gate_axial.py [CONSENSUS_DIR] [RECEPTOR]
        CONSENSUS_DIR  parent dir with per-seed subdirs s<seed>/  (default work/screen_consensus_wt)
        RECEPTOR       structure to read heme iron + pyrrole N from (default work/receptor_A.pdb)
"""
import glob
import math
import os
import pathlib
import statistics
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import protocol
from openafr.pdbqt import iron_position, read_poses
from openafr.scoring import report
from scripts.validate_gate2 import ACTIVES, check_prereg

PREREG = "work/PREREGISTRATION_axial.md"
PREREG_SHA = "3ff59c29cd88c2dfcbca2ddd726494a565b39d5fa7fb647721a0a11a0d1c4e33"
BASELINE_SEED = "42"  # the frozen protocol seed


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _norm(a):
    n = math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
    return (a[0] / n, a[1] / n, a[2] / n)


def heme_normal(receptor_pdb):
    """Unit normal of the porphyrin plane, from HEM A 601's four pyrrole nitrogens.

    n = normalize((NC - NA) x (ND - NB)) — the cross product of the two porphyrin diagonals.
    Deterministic; computed once from the fixed receptor. Its sign is irrelevant downstream
    because the axial angle is folded to the acute value.
    """
    want = {"NA": None, "NB": None, "NC": None, "ND": None}
    for l in open(receptor_pdb):
        if l.startswith("HETATM") and l[17:20].strip() == "HEM":
            name = l[12:16].strip()
            if name in want and want[name] is None:
                want[name] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    missing = [k for k, v in want.items() if v is None]
    if missing:
        raise SystemExit(f"heme pyrrole nitrogens {missing} not found in {receptor_pdb}")
    return _norm(_cross(_sub(want["NC"], want["NA"]), _sub(want["ND"], want["NB"])))


def pose_axial(pose_path, fe, n_hat):
    """Per pose, measured on the nitrogen NEAREST the iron (the atom the distance criterion
    uses), yield (dist, theta_deg, S). Poses with no nitrogen are skipped."""
    out = []
    for _num, atoms in read_poses(pose_path):
        ns = [a for a in atoms if a[0].startswith("N")]
        if not ns:
            continue
        nb = min(ns, key=lambda a: math.dist(a[1:], fe))
        d = math.dist(nb[1:], fe)
        vhat = _norm(_sub(nb[1:], fe))
        cos = abs(vhat[0] * n_hat[0] + vhat[1] * n_hat[1] + vhat[2] * n_hat[2])
        theta = math.degrees(math.acos(min(1.0, max(-1.0, cos))))  # folded to [0, 90]
        out.append((d, theta, d * (1.0 + math.sin(math.radians(theta)))))
    return out


def per_seed_minima(pose_path, fe, n_hat):
    """(min S, min theta, min N-Fe) over all poses in one docked file; each None if the file
    has no nitrogen-bearing pose. The three minima are taken independently per the
    pre-registration (S gates; theta and distance are reported alongside)."""
    vals = pose_axial(pose_path, fe, n_hat)
    if not vals:
        return None, None, None
    return (min(v[2] for v in vals), min(v[1] for v in vals), min(v[0] for v in vals))


def _median(byseed):
    """Median of the non-None per-seed minima, or None if no seed produced a value."""
    vals = [v for v in byseed.values() if v is not None]
    return statistics.median(vals) if vals else None


def main():
    parent = sys.argv[1] if len(sys.argv) > 1 else "work/screen_consensus_wt"
    receptor = sys.argv[2] if len(sys.argv) > 2 else "work/receptor_A.pdb"
    fe = iron_position(receptor)
    n_hat = heme_normal(receptor)
    actives = {l.split("\t")[1].strip() for l in open(ACTIVES) if "\t" in l}

    # per-seed pose dirs are s<digits>/; exclude screen.sh's s<seed>_pdbqt/ input dirs.
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
    print("AXIAL GATE — PRE-REGISTERED ORTHOGONAL (DIRECTION) FEATURE, WILD-TYPE HELD-OUT")
    print("=" * 70)
    print(f"receptor         : {receptor}")
    print(f"heme normal n    : ({n_hat[0]:+.4f}, {n_hat[1]:+.4f}, {n_hat[2]:+.4f})")
    print(f"seeds            : {seeds}  (mode E gates; baseline slice = seed {BASELINE_SEED})")
    check_prereg(PREREG, PREREG_SHA)
    protocol.verify()

    # molecule -> {seed: (minS, minTheta, minDist)}
    per_mol = {}
    for d in seed_dirs:
        seed = os.path.basename(d)[1:]
        for f in sorted(os.listdir(d)):
            if f.endswith(".pdbqt"):
                per_mol.setdefault(f[:-6], {})[seed] = per_seed_minima(
                    os.path.join(d, f), fe, n_hat)

    s_rows, theta_rows, dist_rows, base_rows = [], [], [], []
    med_s = {}
    for name, byseed in per_mol.items():
        is_act = name in actives
        s_by = {s: v[0] for s, v in byseed.items()}
        th_by = {s: v[1] for s, v in byseed.items()}
        d_by = {s: v[2] for s, v in byseed.items()}
        med_s[name] = _median(s_by)
        s_rows.append((name, med_s[name], is_act))
        theta_rows.append((name, _median(th_by), is_act))
        dist_rows.append((name, _median(d_by), is_act))
        base_rows.append((name, s_by.get(BASELINE_SEED), is_act))

    seen = {r[0] for r in s_rows if r[2]}
    if seen != actives:
        print(f"\n!! MISSING ACTIVES: {actives - seen} — results not valid")

    e = rep("MODE E (AXIAL) — median per-seed min S = |v|*(1+sin theta)  [PRIMARY — gates]",
            s_rows, gating=True)
    rep("axial angle alone — median per-seed min theta (deg; direction only)", theta_rows)
    rep("reliability distance — median per-seed min N-Fe (for reference)", dist_rows)
    b = rep(f"seed-{BASELINE_SEED} single-search S — does the single-shot top-1% survive?",
            base_rows)

    print("\nheld-out actives by axial S (median), best first:")
    print("  %-16s %10s" % ("active", "median S"))
    for a in sorted(actives, key=lambda a: med_s.get(a) if med_s.get(a) is not None else 9e9):
        print("  %-16s %10s" % (a, "-" if med_s.get(a) is None else "%.3f" % med_s[a]))

    print("\n" + "=" * 70)
    passed = e[0] >= min_auc and e[1][0] >= min_ef1
    print(f"AXIAL mode E: AUC {e[0]:.3f} (need >= {min_auc})   EF@1% {e[1][0]:.2f}x (need >= {min_ef1})")
    print(f"  (context) seed-{BASELINE_SEED} single-search S: AUC {b[0]:.3f} EF@1% {b[1][0]:.2f}x")
    if passed:
        print("GATE: PASS — the orthogonal axial-direction feature recovers top-1% separation")
        print("that distance-aggregation alone could not. Rescuable by geometry of coordination,")
        print("not distance. Carry forward ONLY under its own separate pre-registration.")
        return 0
    print("GATE: FAIL — adding axial direction still does not separate real azoles from")
    print("property-matched decoys at the top 1%. Per pre-registration the conclusion is")
    print("stronger than reliability's: the decoys coordinate the iron with azole-like axial")
    print("geometry too, so pose geometry (distance + direction) cannot own the extreme top")
    print("here. Defensible product stays a trustworthy top-~5% shortlist. No second feature")
    print("is fished on these poses; a further attempt requires an independent dataset.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
