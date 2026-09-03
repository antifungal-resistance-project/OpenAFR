"""Criterion-fix validation — does restricting the iron-coordinator to an AROMATIC
ring nitrogen preserve the pre-registered run-2 gate?

Pipeline stage: validate/ (a criterion comparison, not a new gate).

Background. The frozen gate ranks by `openafr.pdbqt.min_nitrogen_iron_distance` — the
nearest nitrogen of ANY type to the heme iron. The coordinator-identity axis
(work/RESULTS_coordinator.md) showed that on a NOVEL library this silently rewards
nitrile/amine/azide nitrogens, which is out-of-mechanism (the method was validated on an
aromatic ring N donating to the iron). The principled fix is to rank on
`openafr.coordinator.min_aromatic_nitrogen_iron_distance` instead.

This script re-grades the run-2 held-out gate (7 azole actives vs 350 property-matched
decoys, work/PREREGISTRATION_run2.md) under BOTH criteria on the SAME poses, so the only
variable is which nitrogen the minimum is allowed to use. The validated actives are azoles
(aromatic coordinators) so their distance should be ~unchanged; nitrile/amine-bearing decoys
that scored via a non-aromatic N are pushed away. The pre-registered hypothesis
(work/PREREGISTRATION_coordinator.md) is that the aromatic-restricted AUC does NOT fall below
the any-N AUC and the gate still PASSES (AUC >= 0.70, EF@1% >= 5.0x).

A molecule with no aromatic-ring nitrogen reaching the pocket is `unrankable` under the
restricted criterion and — exactly as the frozen gate treats an active with no usable pose —
ranked LAST, never dropped (dropping unrankable rows is the easiest way to fake a pass).
"""
import argparse
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import protocol
from openafr.coordinator import min_aromatic_nitrogen_iron_distance
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses
from openafr.scoring import report

ACTIVES = "data/ligands/actives_holdout_final.smi"


def _min_over_poses(pose_file, fe, criterion):
    """Smallest distance the criterion returns over all poses, or None if never defined."""
    best = None
    for _num, atoms in read_poses(pose_file):
        d = criterion(atoms, fe)
        if d is not None and (best is None or d < best):
            best = d
    return best


def grade(screen, receptor="work/receptor_A.pdb", actives_file=ACTIVES):
    fe = iron_position(receptor)
    actives = {l.split("\t")[1].strip() for l in open(actives_file) if "\t" in l}
    gate = protocol.load()["gate"]
    min_auc, min_ef1 = gate["min_auc"], gate["min_ef1"]
    ef_fractions, bedroc_alpha = gate["ef_fractions"], gate["bedroc_alpha"]

    any_rows, arom_rows = [], []
    for f in sorted(os.listdir(screen)):
        if not f.endswith(".pdbqt"):
            continue
        name = f[:-6]
        is_act = name in actives
        p = os.path.join(screen, f)
        any_rows.append((name, _min_over_poses(p, fe, min_nitrogen_iron_distance), is_act))
        arom_rows.append((name, _min_over_poses(p, fe, min_aromatic_nitrogen_iron_distance),
                          is_act))

    seen = {r[0] for r in any_rows if r[2]}
    if seen != actives:
        print(f"!! MISSING ACTIVES: {actives - seen} — results not valid")

    print("=" * 70)
    print("CRITERION-FIX VALIDATION — any-N vs aromatic-N on the run-2 held-out gate")
    print("=" * 70)
    print(f"receptor         : {receptor}")
    protocol.verify()
    print(f"held-out actives : {len(actives)}   pool: {len(any_rows)}")
    print(f"pass bar         : AUC >= {min_auc} AND EF@1% >= {min_ef1}x\n")

    a = report("MODE C (any nitrogen) — the FROZEN criterion", any_rows, True,
               ef_fractions, bedroc_alpha)
    b = report("MODE C-arom (aromatic ring nitrogen only) — the FIX", arom_rows, True,
               ef_fractions, bedroc_alpha)

    n_unrankable = sum(1 for r in arom_rows if r[1] is None)
    n_unrankable_act = sum(1 for r in arom_rows if r[1] is None and r[2])
    print("\n" + "=" * 70)
    print(f"any-N       : AUC {a[0]:.3f}   EF@1% {a[1][0]:.2f}x")
    print(f"aromatic-N  : AUC {b[0]:.3f}   EF@1% {b[1][0]:.2f}x   "
          f"(unrankable {n_unrankable}, of which actives {n_unrankable_act})")
    preserved = b[0] >= a[0] - 1e-9
    passes = b[0] >= min_auc and b[1][0] >= min_ef1
    print(f"AUC preserved vs any-N : {'YES' if preserved else 'NO'} "
          f"({b[0]-a[0]:+.3f})")
    print(f"aromatic-N still PASSES : {'YES' if passes else 'NO'}")
    return 0 if (preserved and passes) else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("screen", nargs="?", default="work/screen2",
                    help="docked pose dir (run-2 gate set)")
    ap.add_argument("--receptor", default="work/receptor_A.pdb")
    ap.add_argument("--actives", default=ACTIVES)
    args = ap.parse_args()
    return grade(args.screen, args.receptor, args.actives)


if __name__ == "__main__":
    sys.exit(main())
