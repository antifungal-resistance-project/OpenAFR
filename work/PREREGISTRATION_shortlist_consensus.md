# Pre-registration — seed-stability of the candidate shortlist (issue #68)

**Written 2026-08-26, BEFORE reading any nitrogen-to-iron distance from the candidate dockings.**
Everything below is fixed in advance; deviations are recorded as deviations, not edited in.
Inherits the run-2 protocol and the consensus criterion —
read [PREREGISTRATION_run2.md](PREREGISTRATION_run2.md) and
[PREREGISTRATION_consensus.md](PREREGISTRATION_consensus.md) first.

## Background — why this run exists

The candidate shortlist (`work/RESULTS_candidate_shortlist.md`) was ranked from a screen in
which every molecule was docked **once** (single frozen seed). The reliability result
(`work/RESULTS_reliability.md`) established that a single-search minimum N-Fe is a **noisy
estimator**: flexible molecules' best coordinating pose is hard to find in one search (spread
up to ~1.1 Å across seeds), so a single draw can badly mis-rank a molecule — in the very
repurposing screen these candidates come from, the real azole clotrimazole landed at
#2774/2776. A specific rank such as "verinurad #2" therefore has no established per-molecule
reliability. This run re-docks the shortlist under the **already-validated consensus criterion**
(min N-Fe over 5 seeds) to separate candidates whose coordination is reproducible from those
whose good single-seed number was luck.

This is an **application** of the pre-registered consensus criterion, not a new criterion or a
new gate. Nothing about the ranking math changes.

## Fixed inputs

- **Ligand set:** the 10 novel survivors of the shortlist (verinurad, vatalanib,
  flucloxacillin, taranabant, lersivirine, LDN-27219, R-1479, L-838417, olprinone, nolatrexed)
  plus reference azoles present in the same top-100 (flutrimazole, ravuconazole, letrozole,
  fluconazole, voriconazole) and ketoconazole. Frozen in `work/candidates/shortlist_focus.smi`.
  The reference azoles are context, not the object of the test.
- **Receptor:** the rigid wild-type *C. albicans* CYP51 `work/receptor.pdbqt` (5TZ1 chain A +
  heme), the same receptor the shortlist was ranked on. The auris/mutant receptor is issue #69,
  a separate run.
- **Protocol:** unchanged, hash-frozen `protocol.yaml`. Only the RNG seed varies across the 5
  replicates. Seeds {42, 1, 2, 3, 4}; seed 42 is the frozen protocol seed, reproducing the
  single-seed number as the apples-to-apples baseline.

## The criterion (fixed; already validated)

For each molecule: dock under the 5 seeds, take the per-seed best (minimum) N-Fe, and summarise
with the **minimum across seeds** (consensus) — a less-biased estimate of the achievable
coordination. The run-to-run **spread** (max − min over coordinating seeds) is the confidence
signal. Coordination cutoff FE_CUTOFF = 3.0 Å, unchanged from `scripts/analyze_poses.py`.

## Pre-committed stability rule (the only new decision — frozen here)

A candidate is **STABLE** iff it produces a coordinating pose (N-Fe < 3.0 Å) in **all 5 seeds**
AND its spread ≤ **1.0 Å**. Otherwise it is **UNSTABLE** and demoted from the shortlist:
its tight single-seed distance was not reproducible, so it is not a trustworthy hypothesis on
this evidence. The 1.0 Å bound is set here, before seeing the numbers, from the reliability
result's observed active spread (~1.1 Å is the noisy end); it is not tuned to any candidate.

Computed by `scripts/shortlist_confidence.py` (SPREAD_MAX = 1.0), which reads this file.

## Pre-committed interpretation

- Candidates that are STABLE **and** land in the validated actives' iron-bound band
  (2.47–2.88 Å, `work/RESULTS_all_azoles.md`) are the reinforced shortlist.
- Candidates that are UNSTABLE are demoted with the spread reported — an honest "this rank was
  single-dock noise" outcome, which is a valid and expected result, not a failure of the run.
- The reference azoles are reported for context (do the real drugs coordinate stably here?);
  they do not gate anything.
- No post-hoc change to the stability threshold or the aggregator after seeing the numbers.

## Limitations (fixed in advance)

1. R = 5 reduces but does not eliminate search variance.
2. Same rigid single-conformer *C. albicans* receptor and ≤45-heavy-atom domain as the
   shortlist — this run tests seed-stability only, not the conformer or receptor caveats
   (issues #69, #70).
3. Stable coordination remains **necessary, not sufficient**: property-matched decoys also
   coordinate (the decoy ceiling). This run removes single-dock luck, not the ceiling.
