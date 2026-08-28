"""Pose-convergence diagnostics (issue #71): did a single docking search settle on one
binding geometry, or scatter across the pocket?

Pure geometry — no I/O, no docking. Given the heavy-atom coordinates and Vina affinities of
the poses one search returned, cluster them by unaligned, index-matched RMSD and report how
concentrated the winning geometry is. Every threshold is frozen in
work/PREREGISTRATION_pose_convergence.md (sha work/PREREG_pose_convergence.sha256) BEFORE any
RMSD was read.

Why unaligned RMSD: Vina writes every mode with an identical atom order into the same
receptor frame, so a direct root-mean-square over matched atoms measures whether two poses
occupy the *same location in the pocket* — which is the question. A superposition-minimised
RMSD would call two poses "the same" even after one translated to a different sub-site, hiding
the very dispersion this diagnostic exists to catch.
"""
import math

RMSD_CUT = 2.0            # A; AutoDock pose-cluster tolerance (frozen in the pre-registration)
CONVERGED_FRACTION = 0.5  # strict majority of returned poses must share the winning geometry


def pose_rmsd(a, b):
    """Unaligned, index-matched heavy-atom RMSD (A) between two poses.

    a, b are equal-length lists of (name, x, y, z) as openafr.pdbqt.read_poses yields. Atoms
    are matched by position (Vina's constant per-mode atom order); no superposition.
    """
    if len(a) != len(b):
        raise ValueError(f"pose atom-count mismatch: {len(a)} vs {len(b)}")
    if not a:
        raise ValueError("empty pose")
    sq = sum((ax[1] - bx[1]) ** 2 + (ax[2] - bx[2]) ** 2 + (ax[3] - bx[3]) ** 2
             for ax, bx in zip(a, b))
    return math.sqrt(sq / len(a))


def cluster(poses, scores, cut=RMSD_CUT):
    """Energy-ordered leader clustering (the AutoDock convention).

    poses: list of poses (each a list of atom tuples), in the SAME order as `scores`.
    scores: Vina affinity per pose (kcal/mol), best-first is NOT assumed — sorted here.
    Returns a list of clusters, each a dict {leader, members, leader_score}, ordered by
    ascending leader affinity (best cluster first). `leader`/`members` are indices into the
    input `poses`. Processing poses best-affinity first, each pose joins the first existing
    leader within `cut`, else becomes a new leader.
    """
    order = sorted(range(len(poses)), key=lambda i: scores[i])
    clusters = []
    for i in order:
        placed = False
        for c in clusters:
            if pose_rmsd(poses[i], poses[c["leader"]]) <= cut:
                c["members"].append(i)
                placed = True
                break
        if not placed:
            clusters.append({"leader": i, "members": [i], "leader_score": scores[i]})
    return clusters


def convergence(poses, scores, coordinating_index, cut=RMSD_CUT,
                converged_fraction=CONVERGED_FRACTION):
    """Full per-candidate convergence record from one search's poses (pure).

    poses / scores: as in cluster(). coordinating_index: index (into poses) of the pose that
    gives the reported min N-Fe distance — the geometry the criterion selected — or None if
    no pose coordinated. Returns the frozen metrics + CONVERGED/DISPERSED verdict.
    """
    n = len(poses)
    if n == 0:
        return None
    clusters = cluster(poses, scores, cut)
    # The top cluster is led by the globally best-affinity pose (clusters[0], since leaders
    # are added in ascending-affinity order and the first pose processed is the best).
    top = clusters[0]
    top_pop = len(top["members"])
    top_fraction = top_pop / n
    top_radius = max((pose_rmsd(poses[m], poses[top["leader"]]) for m in top["members"]),
                     default=0.0)
    energy_gap = (clusters[1]["leader_score"] - top["leader_score"]
                  if len(clusters) > 1 else None)

    # Rank (1-based, by ascending affinity) of the coordinating pose, and whether it sits in
    # the dominant cluster.
    rank_of = {idx: r + 1 for r, idx in enumerate(sorted(range(n), key=lambda i: scores[i]))}
    if coordinating_index is None:
        coord_rank = None
        coord_in_top = None
    else:
        coord_rank = rank_of[coordinating_index]
        coord_in_top = pose_rmsd(poses[coordinating_index], poses[top["leader"]]) <= cut

    verdict = "CONVERGED" if top_fraction >= converged_fraction else "DISPERSED"
    return {
        "n_poses": n,
        "n_clusters": len(clusters),
        "top_cluster_pop": top_pop,
        "top_cluster_fraction": round(top_fraction, 3),
        "top_cluster_radius": round(top_radius, 3),
        "energy_gap": None if energy_gap is None else round(energy_gap, 3),
        "coordinating_pose_rank": coord_rank,
        "coordinating_pose_in_top_cluster": coord_in_top,
        "verdict": verdict,
    }
