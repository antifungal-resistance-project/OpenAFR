"""Per-candidate pose-convergence diagnostics (issue #71).

Grades whether each shortlist candidate's docking search CONVERGED on one binding geometry
or DISPERSED across the pocket — a per-candidate CONFIDENCE AXIS, not a gate. It
demotes/annotates, never promotes. Criterion + bar are frozen in
work/PREREGISTRATION_pose_convergence.md (sha work/PREREG_pose_convergence.sha256) BEFORE any
pose-to-pose RMSD was read.

Method (all frozen in the pre-registration):
  - Read the 20 poses of the ONE frozen seed-42 dock per candidate (num_modes 20,
    exhaustiveness 32, the frozen protocol) with their Vina affinities.
  - Cluster them by unaligned, index-matched heavy-atom RMSD at 2.0 A (AutoDock convention),
    energy-ordered leader clustering. The top cluster is led by the best-affinity pose.
  - CONVERGED iff >= 50% of the returned poses fall in that top cluster; else DISPERSED.
  - Also report where the coordinating (min N-Fe) pose sits: its affinity rank and whether it
    lies inside the dominant cluster.

This is the WITHIN-search complement to the ACROSS-seed stability axis (#68): a candidate can
be seed-stable yet dispersed, or seed-unstable yet locally converged. Both belong on the card.

Docking (when --from-screen is not given) reproduces the frozen protocol exactly via
scripts/screen.sh; because the seed is fixed, the min N-Fe over the 20 poses here must equal
the `per_seed_fungal["42"]` value the seed-stability run recorded — a built-in reproduction
check, printed as `matches_seed42_baseline`.

Usage:
    # dock fresh (needs vina + obabel + a prepared receptor, the pinned openafr env):
    python scripts/pose_convergence.py \
        --ligands work/candidates/cand_ligands --workdir work/candidates/dock_convergence \
        [--receptor-pdbqt work/receptor.pdbqt] [--receptor-pdb work/receptor_A.pdb] \
        [--baseline work/candidates/shortlist_confidence.json] [--json]
    # or read poses already docked into a screen dir:
    python scripts/pose_convergence.py --from-screen work/candidates/dock_convergence [...]
"""
import argparse
import json
import os
import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses, read_scores
from openafr.pose_convergence import RMSD_CUT, convergence
from scripts.shortlist_confidence import FE_CUTOFF

ROOT = pathlib.Path(__file__).resolve().parent.parent


def coordinating_pose_index(poses, fe):
    """Index of the pose with the smallest nitrogen-to-iron distance under FE_CUTOFF, or None
    (the geometry the criterion selects). Mirrors scripts/shortlist_confidence.best_distance,
    but keeps the pose index so its cluster membership can be checked."""
    best_i, best_d = None, None
    for i, atoms in enumerate(poses):
        d = min_nitrogen_iron_distance(atoms, fe)
        if d is not None and d < FE_CUTOFF and (best_d is None or d < best_d):
            best_i, best_d = i, d
    return best_i, best_d


def assess_file(pose_path, fe):
    """One candidate's convergence record from its docked .pdbqt, or None if unreadable."""
    poses = [atoms for _num, atoms in read_poses(pose_path)]
    scores = read_scores(pose_path)
    if not poses or len(scores) != len(poses):
        return None
    coord_i, coord_d = coordinating_pose_index(poses, fe)
    rec = convergence(poses, scores, coord_i)
    rec["min_n_fe"] = None if coord_d is None else round(coord_d, 3)
    return rec


def dock(ligands, workdir, receptor_pdbqt):
    """Run the frozen seed-42 screen over the prepared ligands into workdir. Deterministic."""
    os.makedirs(workdir, exist_ok=True)
    subprocess.run([str(ROOT / "scripts" / "screen.sh"), ligands, workdir],
                   check=True, env={**os.environ, "RECEPTOR": receptor_pdbqt})
    return workdir


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ligands", default="work/candidates/cand_ligands")
    ap.add_argument("--workdir", default="work/candidates/dock_convergence")
    ap.add_argument("--from-screen", metavar="SCREEN_DIR",
                    help="skip docking; read poses already docked into SCREEN_DIR")
    ap.add_argument("--receptor-pdbqt", default="work/receptor.pdbqt")
    ap.add_argument("--receptor-pdb", default="work/receptor_A.pdb")
    ap.add_argument("--baseline", default="work/candidates/shortlist_confidence.json",
                    help="seed-stability JSON for the seed-42 reproduction check")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    screen_dir = args.from_screen or dock(args.ligands, args.workdir, args.receptor_pdbqt)
    fe = iron_position(args.receptor_pdb)

    # seed-42 baseline min N-Fe from the seed-stability run, for the reproduction check.
    baseline = {}
    if args.baseline and pathlib.Path(args.baseline).is_file():
        for r in json.load(open(args.baseline))["candidates"]:
            v = (r.get("per_seed_fungal") or {}).get("42")
            if v is not None:
                baseline[r["name"]] = v

    recs = []
    for f in sorted(os.listdir(screen_dir)):
        if not f.endswith(".pdbqt"):
            continue
        name = f[:-6]
        rec = assess_file(os.path.join(screen_dir, f), fe)
        if rec is None:
            continue
        rec = {"name": name, **rec}
        base = baseline.get(name)
        if base is not None and rec["min_n_fe"] is not None:
            rec["seed42_baseline"] = base
            rec["matches_seed42_baseline"] = abs(rec["min_n_fe"] - base) <= 0.01
        recs.append(rec)
    # Sort by the reported min N-Fe (the shortlist order); unrankable last.
    recs.sort(key=lambda r: (r["min_n_fe"] is None, r["min_n_fe"] or 9e9))

    if args.json:
        print(json.dumps({"rmsd_cut": RMSD_CUT, "fe_cutoff": FE_CUTOFF,
                          "candidates": recs}, indent=2))
        return 0

    print("=" * 104)
    print("CANDIDATE POSE CONVERGENCE — within-search pose scatter, frozen seed-42 dock (#71)")
    print(f"receptor: {args.receptor_pdb}   RMSD cut {RMSD_CUT} A   "
          f"CONVERGED = >=50% of poses in the best-affinity cluster")
    print("=" * 104)
    hdr = ("candidate", "N-Fe", "poses", "top%", "clust", "radius", "gap", "coord_rk",
           "in_top", "verdict")
    print("%-16s %6s %6s %6s %6s %7s %6s %8s %7s   %s" % hdr)
    print("-" * 104)
    for r in recs:
        nfe = "-" if r["min_n_fe"] is None else "%.3f" % r["min_n_fe"]
        gap = "-" if r["energy_gap"] is None else "%.2f" % r["energy_gap"]
        ck = "-" if r["coordinating_pose_rank"] is None else str(r["coordinating_pose_rank"])
        it = "-" if r["coordinating_pose_in_top_cluster"] is None else (
            "yes" if r["coordinating_pose_in_top_cluster"] else "NO")
        print("%-16s %6s %6d %5.0f%% %6d %7.3f %6s %8s %7s   %s" % (
            r["name"], nfe, r["n_poses"], 100 * r["top_cluster_fraction"], r["n_clusters"],
            r["top_cluster_radius"], gap, ck, it, r["verdict"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
