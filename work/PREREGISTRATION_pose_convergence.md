# Pre-registration — per-candidate pose-convergence diagnostics (issue #71)

**Written 2026-08-27, BEFORE reading any pose-to-pose RMSD or cluster count from the
candidate dockings.** Everything below is fixed in advance; deviations are recorded as
deviations, not edited in. Inherits the run-2 protocol and the coordination criterion —
read [PREREGISTRATION_run2.md](PREREGISTRATION_run2.md) and
[PREREGISTRATION_consensus.md](PREREGISTRATION_consensus.md) first.

## Background — why this run exists

The consolidated scorecard ([RESULTS_scorecard.md](RESULTS_scorecard.md)) reports, per
candidate, the **single closest** nitrogen-to-iron distance the docking found — the best
pose. The seed-stability axis (#68) asks whether that number reproduces **across independent
searches** (different RNG seeds). Neither asks the complementary question **within one
search**: of the 20 poses Vina returns for a candidate, do they *agree* on where the ligand
sits, or does the winning pose stand almost alone while the rest scatter across the pocket?

A tight N–Fe from one pose means little if the other 19 poses cluster somewhere else — the
search did not converge, and the reported distance is a lucky outlier, not a basin. This is
exactly the failure the reliability result warned about (a single dock is a noisy estimator),
seen from *inside* a single run rather than across seeds.

This is a **per-candidate confidence axis, not a gate.** It does not justify or block a
screen. It adds one honest column per candidate to the scorecard: did the search *converge*
on the geometry we reported?

## Fixed inputs

- **Ligand set:** identical to the seed-stability, human-selectivity and resistance runs —
  `work/candidates/shortlist_focus.smi` (10 novel candidates + 6 known-azole controls), the
  same prepared 3D ligands (`work/candidates/cand_ligands/`, built by
  `scripts/prep_ligands.py`, ETKDGv3 seed 42).
- **Receptor:** `work/receptor.pdbqt` (docking) / `work/receptor_A.pdb` (heme-iron anchor) —
  the wild-type *C. albicans* 5TZ1 stand-in, the receptor the shortlist was ranked on.
  Rebuilt deterministically by `scripts/prep_receptor.py` (receptor.pdbqt sha256
  `7ed596cb906f…`, anchor unchanged).
- **Docking box + search:** the frozen `protocol.yaml` **unchanged** — center
  70.61/66.28/4.18, size 26, exhaustiveness 32, num_modes 20.
- **Seed:** **42, the single frozen protocol seed** — the one canonical search. This axis
  measures *within-search* pose scatter, so it reads the 20 poses of exactly one dock per
  candidate (not the 5-seed consensus, which is the separate #68 axis). Because docking is
  deterministic at a fixed seed, this seed-42 dock reproduces the exact poses whose min the
  seed-stability run recorded as `per_seed_fungal["42"]` — a built-in reproduction check
  (the min N–Fe over the 20 poses here must equal that stored value).

## The measure (fixed)

For each candidate, read all poses Vina returned (`openafr.pdbqt.read_poses`, heavy atoms
only) with their Vina affinities (`read_scores`, best-first file order).

- **Pose RMSD** between two poses = root-mean-square over the heavy-atom set of the
  per-atom Euclidean distance, atoms matched **by index** (Vina writes every mode with an
  identical atom order), with **NO superposition**. The poses already share the receptor
  coordinate frame; a direct (unaligned) RMSD is deliberate — it measures whether two poses
  occupy the *same location in the pocket*. A superposition-minimised RMSD would call two
  poses "the same" even if one had translated to a different sub-site, which is the very
  thing this diagnostic must catch.
- **Clustering** = energy-ordered leader clustering at `RMSD_CUT = 2.0 Å` (the standard
  AutoDock pose-cluster tolerance). Process poses in ascending affinity (best first); each
  pose joins the first already-established cluster leader within `RMSD_CUT`, otherwise it
  becomes a new leader. The **top cluster** is the one led by the rank-1 (best-affinity)
  pose.

Per candidate, frozen metrics:

- `n_poses` — poses returned (≤ 20).
- `n_clusters` — number of leader clusters at 2.0 Å.
- `top_cluster_pop` / `top_cluster_fraction` — poses in the rank-1 cluster, and that over
  `n_poses`.
- `top_cluster_radius` — max member-to-leader RMSD inside the top cluster (Å).
- `energy_gap` — affinity of the next distinct cluster's leader minus the rank-1 affinity
  (kcal/mol, ≥ 0; the energetic separation between the winning geometry and the best
  *alternative* geometry). `None` if there is only one cluster (fully converged).
- `coordinating_pose_rank` — the affinity rank (1-based) of the pose that gives the
  reported min N–Fe. Rank 1 means the best-energy pose *is* the tightest coordinator.
- `coordinating_pose_in_top_cluster` — whether that min-N–Fe pose lies within `RMSD_CUT`
  of the rank-1 pose (i.e. the geometry the criterion selected is part of the dominant
  basin, not an isolated outlier).

## Pre-committed per-candidate classification

Frozen here, before any RMSD is read. `CONVERGED_FRACTION = 0.5` — a **strict majority** of
the returned poses must share the winning geometry.

- **CONVERGED** — `top_cluster_fraction ≥ 0.5`: at least half of the returned poses cluster
  with the best-energy pose. The search settled on one basin.
- **DISPERSED** — `top_cluster_fraction < 0.5`: the winning pose is a minority geometry; the
  returned poses scatter. Lower confidence that the reported N–Fe is a robust minimum rather
  than a search outlier. This is the informative demotion this axis exists to surface.

Reported alongside, not folded into the flag: `coordinating_pose_in_top_cluster`. A candidate
can be CONVERGED yet have its *coordinating* pose sit outside the dominant cluster
(`coordinating_pose_in_top_cluster = false`) — the search converged on a geometry that is
**not** the one the N–Fe criterion selected. That disagreement is annotated verbatim; it is
a distinct, weaker concern than DISPERSED and is never silently merged into it.

## Pre-committed interpretation (do NOT overclaim)

- **Convergence is pose *reproducibility within one search*, not thermodynamic stability and
  not activity.** A CONVERGED verdict means the docking search agreed with itself on where
  the ligand sits; it says nothing about whether that pose is the true bound state, nor
  whether the molecule is active. Coordination remains necessary, not sufficient
  ([RESULTS_axial.md](RESULTS_axial.md)).
- **DISPERSED is the informative outcome.** It flags a candidate whose reported tight N–Fe
  came from a pose the search could not reproduce even once more among 19 alternatives — a
  real confidence demotion, consistent with the seed-stability finding (#68) that single-dock
  ranks can be search luck.
- **Complementary to #68, not a replacement.** Seed-stability asks *across* independent
  searches; this asks *within* one. A candidate can be seed-stable yet dispersed, or
  seed-unstable yet locally converged — both facts belong on the card.
- This axis **demotes/annotates, never promotes.** It is added to the scorecard as one
  column with a `not-yet-checked` → measured state; it does not collapse into any composite
  score and does not make any molecule a hit.

## Limitations (fixed in advance)

1. **Single frozen seed.** Within-search convergence is read from the one seed-42 dock. The
   across-seed picture is the separate #68 axis; this does not re-measure it.
2. **Unaligned, index-matched RMSD** assumes Vina's constant atom order across modes (it
   holds) and treats symmetric atoms distinctly — a ligand that is chemically the same after
   a symmetry swap can read as a higher RMSD than a symmetry-aware metric would. This is
   conservative (it can call converged poses dispersed, never the reverse) and is accepted
   to keep the metric a pure geometry read with no alignment degrees of freedom.
3. **Single conformer, rigid receptor**, coordination necessary-not-sufficient — every
   caveat of the underlying screen still binds.
4. **Controls are controls.** The 6 known azoles are here to show the metric behaves
   sensibly (a floppy azole like ketoconazole is *expected* to disperse), not to test a
   discovery.
