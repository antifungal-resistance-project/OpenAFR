"""Known-answer tests for the pose-convergence logic (openafr/pose_convergence.py, #71).

These lock the pure decision rules frozen in
work/PREREGISTRATION_pose_convergence.md:
  - unaligned, index-matched heavy-atom RMSD (no superposition);
  - energy-ordered leader clustering at 2.0 A; the top cluster is the best-affinity pose's;
  - CONVERGED iff >= 50% of returned poses fall in that top cluster, else DISPERSED;
  - the coordinating (min N-Fe) pose's rank and its membership in the top cluster are
    reported alongside the verdict, never folded into it.
"""
import math

import pytest

from openafr.pose_convergence import cluster, convergence, pose_rmsd


def _pose(offset, n=4):
    """A pose of n atoms translated by `offset` along x from the origin block."""
    return [("C%d" % i, float(i) + offset, 0.0, 0.0) for i in range(n)]


def test_rmsd_unaligned_translation_is_the_offset():
    # a pure translation by d along x -> RMSD == d (no superposition removes it)
    assert pose_rmsd(_pose(0.0), _pose(1.5)) == pytest.approx(1.5)


def test_rmsd_zero_for_identical_pose():
    assert pose_rmsd(_pose(0.0), _pose(0.0)) == 0.0


def test_rmsd_atom_count_mismatch_raises():
    with pytest.raises(ValueError):
        pose_rmsd(_pose(0.0, n=4), _pose(0.0, n=5))


def test_clustering_groups_poses_within_cut():
    # three poses near origin (within 2 A) + one far away (5 A) -> two clusters
    poses = [_pose(0.0), _pose(0.5), _pose(1.0), _pose(5.0)]
    scores = [-9.0, -8.5, -8.0, -7.0]         # best-first order already
    clusters = cluster(poses, scores, cut=2.0)
    assert len(clusters) == 2
    top = clusters[0]
    assert sorted(top["members"]) == [0, 1, 2]   # the three near-origin poses
    assert top["leader"] == 0                     # led by the best-affinity pose


def test_leader_is_best_affinity_not_file_order():
    # the globally best pose is last in input order; it must still lead the top cluster
    poses = [_pose(5.0), _pose(0.5), _pose(0.0)]
    scores = [-7.0, -8.0, -9.5]
    clusters = cluster(poses, scores, cut=2.0)
    assert clusters[0]["leader"] == 2             # index of the -9.5 pose
    assert sorted(clusters[0]["members"]) == [1, 2]


def test_converged_when_majority_share_top_cluster():
    # 8 of 10 poses cluster with the winner -> CONVERGED, energy gap to the lone far cluster
    poses = [_pose(0.1 * i) for i in range(8)] + [_pose(10.0), _pose(10.2)]
    scores = [-9.0 - 0.01 * i for i in range(8)] + [-6.0, -5.9]
    rec = convergence(poses, scores, coordinating_index=0)
    assert rec["verdict"] == "CONVERGED"
    assert rec["top_cluster_pop"] == 8
    assert rec["top_cluster_fraction"] == 0.8
    assert rec["n_clusters"] == 2
    assert rec["energy_gap"] == pytest.approx(3.07, abs=0.01)  # -6.0 - (best -9.07)


def test_dispersed_when_winner_is_minority():
    # winner stands alone; the other 4 poses scatter into their own singletons -> DISPERSED
    poses = [_pose(0.0), _pose(5.0), _pose(10.0), _pose(15.0), _pose(20.0)]
    scores = [-9.0, -8.0, -7.0, -6.0, -5.0]
    rec = convergence(poses, scores, coordinating_index=0)
    assert rec["verdict"] == "DISPERSED"
    assert rec["top_cluster_fraction"] == 0.2     # 1 of 5


def test_coordinating_pose_outside_top_cluster_is_flagged():
    # search converges on a geometry (poses 0,1,2) but the coordinating pose is the far one
    poses = [_pose(0.0), _pose(0.5), _pose(1.0), _pose(9.0)]
    scores = [-9.0, -8.5, -8.0, -7.0]
    rec = convergence(poses, scores, coordinating_index=3)
    assert rec["verdict"] == "CONVERGED"                    # 3/4 share the top cluster
    assert rec["coordinating_pose_in_top_cluster"] is False  # ...but not the coordinating one
    assert rec["coordinating_pose_rank"] == 4                # worst affinity of the four


def test_single_cluster_has_no_energy_gap():
    poses = [_pose(0.1 * i) for i in range(5)]
    scores = [-9.0 - 0.1 * i for i in range(5)]
    rec = convergence(poses, scores, coordinating_index=0)
    assert rec["n_clusters"] == 1
    assert rec["energy_gap"] is None
    assert rec["coordinating_pose_in_top_cluster"] is True


def test_no_coordinating_pose_reports_none():
    poses = [_pose(0.0), _pose(0.5)]
    scores = [-9.0, -8.0]
    rec = convergence(poses, scores, coordinating_index=None)
    assert rec["coordinating_pose_rank"] is None
    assert rec["coordinating_pose_in_top_cluster"] is None
