"""Multi-baseline comparison gate — implements work/PREREGISTRATION_baselines_multi.md (#83).

Quantifies the geometric iron-approach criterion against a TABLE of independent baselines a
skeptic would reach for, all on the IDENTICAL held-out poses (no re-docking, only the scorer
changes). Extends the single-opponent GNINA gate (scripts/rescore_gnina.py, #63) to the
registry in openafr.baselines.

    C (ours)   iron filter (N-Fe < 3.0 A), then best Vina affinity among survivors
               — the run-2 PRIMARY criterion, reproduced on these poses
    <each baseline>  best score over ALL poses — the "should I just use this tool instead of
               your whole pipeline?" question, one row per opponent

Endpoints AUC and EF@5% (EF@1% excluded — RESULTS_reliability), each with a bootstrap 95% CI,
from the same openafr.scoring functions every gate uses. Unrankable molecules (no score) are
ranked LAST, never dropped.

The scorer binaries (GNINA, PLANTS) are Linux/x86 native and run on a cloud host via
scripts/run_gnina_rescore.sh / run_plants_rescore.sh; this grader reads their saved outputs
and is reproducible anywhere. A baseline whose output dir is absent is reported as
`not-yet-run (cloud-gated)` — the current state until the cloud run fills it — so mode C and
any baselines already present still grade.

Usage:
    python scripts/rescore_baselines.py [POSE_DIR] [RECEPTOR_PDB]
        POSE_DIR      dir of Vina .pdbqt poses      (default work/screen2)
        RECEPTOR_PDB  structure to read the iron from (default work/receptor_A.pdb)
"""
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import baselines, protocol
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses, read_scores
from openafr.scoring import bootstrap_auc_ci, report
from scripts.validate_gate2 import ACTIVES, check_prereg

PREREG = "work/PREREGISTRATION_baselines_multi.md"
PREREG_SHA = "b7cbb8df4770fe84920e1d130fa7eedabfe896caee82f42f97240fc145766d46"
FE_CUTOFF = 3.0


def _ours_key(pose_path, fe):
    """Mode C: best (most negative) Vina affinity among the iron-bound poses, or None."""
    poses = list(read_poses(pose_path))
    vina = read_scores(pose_path)
    if len(vina) != len(poses):
        return None
    best = None
    for (_num, atoms), score in zip(poses, vina):
        d = min_nitrogen_iron_distance(atoms, fe)
        if d is not None and d < FE_CUTOFF and (best is None or score < best):
            best = score
    return best


def _baseline_key(spec, name, pose_path):
    """One opponent's molecule key for a ligand, or None if its output is absent/misaligned.

    Aligned 1:1 with the poses by input-file order; a pose/score count mismatch degrades to
    None (ranked last) rather than silently mis-scoring — the same honesty as the GNINA gate.
    """
    out = os.path.join(spec.subdir, name + spec.ext)
    if not os.path.exists(out):
        return None
    per_pose = spec.parse(open(out).read())
    n_poses = sum(1 for _ in read_poses(pose_path))
    if len(per_pose) != n_poses:
        print(f"  !! {spec.key}: pose/score count mismatch for {name} "
              f"({n_poses} poses vs {len(per_pose)} scores) — excluded for it")
        return None
    return baselines.molecule_key(per_pose)


def _ordered_flags(rows):
    ok = sorted((r for r in rows if r[1] is not None), key=lambda r: r[1])
    return [r[2] for r in ok] + [r[2] for r in rows if r[1] is None]


def main():
    pose_dir = sys.argv[1] if len(sys.argv) > 1 else "work/screen2"
    receptor = sys.argv[2] if len(sys.argv) > 2 else "work/receptor_A.pdb"
    if not os.path.isdir(pose_dir):
        sys.exit(f"no pose dir: {pose_dir} — dock the held-out set first (scripts/run_screen2.sh)")

    fe = iron_position(receptor)
    actives = {l.split("\t")[1].strip() for l in open(ACTIVES) if "\t" in l}
    gate = protocol.load()["gate"]
    ef_fractions, bedroc_alpha = gate["ef_fractions"], gate["bedroc_alpha"]
    ef5 = ef_fractions.index(0.05)
    rep = lambda label, rows, g=False: report(label, rows, g, ef_fractions, bedroc_alpha)

    print("=" * 70)
    print("MULTI-BASELINE COMPARISON — geometry vs the tools a skeptic reaches for (#83)")
    print("=" * 70)
    print(f"receptor         : {receptor}   heme Fe ({fe[0]:.2f}, {fe[1]:.2f}, {fe[2]:.2f})")
    print(f"poses            : {pose_dir}   metal cutoff N-Fe < {FE_CUTOFF} A")
    check_prereg(PREREG, PREREG_SHA)
    protocol.verify()

    names = [f[:-6] for f in sorted(os.listdir(pose_dir)) if f.endswith(".pdbqt")]

    # Mode C (ours) — always available; reproduced from the poses.
    c_rows = [(n, _ours_key(os.path.join(pose_dir, n + ".pdbqt"), fe), n in actives)
              for n in names]
    seen = {r[0] for r in c_rows if r[2]}
    if seen != actives:
        print(f"\n!! actives missing from the poses: {actives - seen}")

    table = []  # (label, auc, ef5, ci_lo, ci_hi, metal_aware)
    c = rep("MODE C (OURS) — iron filter, then Vina affinity  [run-2 primary]", c_rows, g=True)
    clo, chi = bootstrap_auc_ci(_ordered_flags(c_rows))
    table.append(("C (ours) — geometry", c[0], c[1][ef5], clo, chi, True))

    pending = []
    for spec in baselines.REGISTRY:
        if not os.path.isdir(spec.subdir):
            pending.append(spec)
            continue
        rows = [(n, _baseline_key(spec, n, os.path.join(pose_dir, n + ".pdbqt")), n in actives)
                for n in names]
        if not any(r[1] is not None for r in rows):
            pending.append(spec)
            continue
        b = rep(f"BASELINE {spec.label}  [{spec.citation}]", rows)
        lo, hi = bootstrap_auc_ci(_ordered_flags(rows))
        table.append((f"{spec.label}", b[0], b[1][ef5], lo, hi, spec.metal_aware))

    print("\n" + "=" * 70)
    print("COMPARISON TABLE  (identical poses; AUC + EF@5%, bootstrap 95% CI)")
    print(f"  {'method':28} {'metal?':6} {'AUC':>6} {'EF@5%':>7}   95% CI")
    for label, auc, e5, lo, hi, metal in table:
        print(f"  {label:28} {'yes' if metal else 'no':6} {auc:6.3f} {e5:6.2f}x   "
              f"[{lo:.3f}, {hi:.3f}]")
    for spec in pending:
        print(f"  {spec.label:28} {'yes' if spec.metal_aware else 'no':6} "
              f"{'not-yet-run (cloud-gated)':>28}")

    if pending:
        print(f"\n{len(pending)} baseline(s) not yet run — Linux/x86 binaries. Fill on a cloud host:")
        for spec in pending:
            runner = f"scripts/run_{spec.key}_rescore.sh"
            print(f"  {spec.key:8} -> {runner}  ->  outputs in {spec.subdir}/")
        print("Then re-run this grader; the comparison is incomplete until then (no headline).")
        return 2

    print("\nAll baselines present — comparison complete. Interpret per "
          "work/PREREGISTRATION_baselines_multi.md (overlap in CIs = tie, not a win).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
