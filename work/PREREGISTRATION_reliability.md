# Pre-registration — reliability (pool-size-neutral) iron-approach criterion

**Written 2026-08-11, BEFORE computing the criterion on the existing poses.** Everything
below is fixed in advance. Inherits the run-2 protocol and reuses the already-docked
wild-type 5-seed poses from the consensus run (`work/screen_consensus_wt/`, no new docking).
Read [PREREGISTRATION_consensus.md](PREREGISTRATION_consensus.md) and
[RESULTS_consensus.md](RESULTS_consensus.md) first.

## Background — the specific flaw this fixes

The consensus criterion (minimum N-Fe over R=5 seeds) FAILED the wild-type gate: it improved
AUC / EF@5% / EF@10% / BEDROC but EF@1% collapsed 12.68x → 0.00x. Diagnosed cause
(`work/RESULTS_consensus.md`): `min`-over-R is an **extreme-value statistic**, so a pool of
348 decoys each getting 5 searches produces a few decoys with an unusually close pose that
win the very top over 7 actives — a **pool-size bias**, not a real geometry signal.

## The reliability criterion (fixed; chosen to remove that specific bias)

For each molecule, take the **per-seed best (minimum) N-Fe distance** — one value per seed,
R=5 values — and rank by their **median** (smaller = better).

- Rationale: the median of a fixed number of draws is a **pool-size-neutral** estimator of
  a molecule's *typical* achievable coordination. Unlike `min`, it does not reward a
  molecule (or a large class of molecules) for extra chances at a lucky extreme — it asks
  "how reliably does this molecule reach the iron," which is the property a true inhibitor
  should have and the flaw image `work/views/` suggests distinguishes drugs from look-alikes.
- Chosen because it is the textbook robust central estimator and directly neutralises the
  diagnosed extreme-value bias — NOT by comparing which statistic scores best.
- **Only the aggregation changes.** Same poses, same box/exhaustiveness/modes/seeds, same
  receptor. R=5, seeds {42,1,2,3,4}, unchanged from the consensus run.

## Primary hypothesis (H1-reliability)

On the run-2 wild-type held-out set (same 7 actives, same decoys, `work/receptor.pdbqt`),
ranking by the median-of-per-seed-minima criterion separates actives from decoys with
**AUC ≥ 0.70 AND EF@1% ≥ 5.0x** (mode C). Same bar as run 2; not lowered.

## Analysis plan

    Primary (gates):  mode C reliability — median of the 5 per-seed minimum N-Fe distances
    Reported alongside (do NOT gate): consensus min-over-seeds, and the seed-42 single search
    Metrics: AUC, EF@1/5/10%, BEDROC(alpha=20) from rdkit.ML.Scoring.
    Any active with no coordinating pose in any seed is ranked LAST, never dropped.

Graded by `scripts/validate_gate_reliability.py`, which re-checks this document's SHA-256 and
the frozen protocol hash before reporting.

## Pre-committed interpretation — and a bound on forking paths

This is the **SECOND and FINAL** aggregator tested on these wild-type poses. Trying
estimators until one passes would be p-hacking; this pre-registration bounds that explicitly:

- **PASS.** A pool-size-neutral robust estimator recovers top-1% separation. Then the
  consensus/Y132F top-1% failures were an artifact of the `min` estimator's bias, and this
  criterion becomes the candidate to carry forward — validated here, then applied elsewhere
  only under its own separate pre-registration. It does NOT retroactively change the
  consensus or Y132F results, which stand as graded.
- **FAIL.** A robust, bias-free geometric estimator *still* cannot separate real azoles from
  property-matched decoys at the top 1%. The committed conclusion is then: **on this decoy
  set, iron-approach geometry has reached its ceiling at the extreme top** — consistent with
  the decoys being genuine iron-coordinators (visible in `work/views/`), not a fixable noise
  problem. **No third aggregator will be tried on these poses.** The honest next step would
  be a different *kind* of signal (e.g. an orthogonal feature), not another way to average
  the same distances.

Whichever occurs is reported as the outcome; no post-hoc statistic invented after seeing
these numbers is substituted to manufacture a pass. The multiple-comparison history (this is
the 2nd aggregator on this data) is stated plainly in the results.

## Limitations (fixed in advance)

1. R = 5 gives a 5-point median — coarse; the estimate itself has sampling noise.
2. Reusing the same wild-type poses as the consensus run means this is a second look at one
   dataset (bounded above to exactly two aggregators, then stop).
3. Same rigid-receptor applicability domain and presumed-inactive decoys as run 2.
