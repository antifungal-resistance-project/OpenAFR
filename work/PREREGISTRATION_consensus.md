# Pre-registration — consensus (variance-reduced) iron-approach criterion

**Written 2026-08-06, BEFORE running the consensus wild-type gate.** Everything below is
fixed in advance. Deviations must be recorded as deviations, not edited in. Inherits the
run-2 protocol; read [PREREGISTRATION_run2.md](PREREGISTRATION_run2.md) first.

## Background

The run-2 criterion ranks molecules by the minimum nitrogen-to-heme-iron distance over a
single Vina search. An **exploratory** diagnostic (`.context/consensus_diagnostic.md`) on
the failed Y132F top-8 showed this single-search minimum is a **noisy estimator**: a
flexible active's best coordinating pose is hard to find in one search (efinaconazole's
per-seed min ranged 2.58–3.70 Å, spread 1.12 Å), while rigid property-matched decoys sample
their floor reproducibly (spread ~0.03–0.45 Å). A single draw can therefore mis-rank a real
drug below steady decoys.

That diagnostic was exploratory and used Y132F molecules (including efinaconazole, which is
also a wild-type held-out active). It is **not** evidence. This document fixes a consensus
method and tests it, confirmatorily, on the wild-type held-out gate.

## The consensus criterion (fixed; chosen on principle, not tuned to any metric)

For each molecule, dock it under **R independent Vina searches with different RNG seeds**,
compute the per-seed best (minimum) N-Fe distance, and rank by the **minimum across seeds**.

- Rationale: this is a less-biased estimate of the *achievable* coordination — exactly what
  the original criterion always meant to measure. Taking the min (not the mean) is chosen
  because the criterion is defined as an achievable minimum; the mean would penalise
  flexible actives for searches that miss the pose. This choice is made on definition, not
  by comparing which aggregator scores best.
- **Only the RNG seed varies.** Box (26 Å @ 70.61/66.28/4.18), exhaustiveness (32) and
  num_modes (20) stay from the hash-frozen `protocol.yaml`. This is a better estimator of
  the same quantity, not a new protocol.
- **R = 5**, seeds {42, 1, 2, 3, 4}. R is set by compute budget, not optimised against
  results. Seed 42 is the frozen protocol seed, so the seed-42 slice reproduces the
  single-search baseline for an apples-to-apples comparison graded in the same run.
- The per-molecule run-to-run **spread** (max − min) is reported as a confidence signal; it
  does not enter the ranking or the gate.

## Primary hypothesis (H1-consensus)

On the run-2 wild-type held-out set (same 7 actives, same decoys, same rigid 5TZ1 + heme
receptor `work/receptor.pdbqt`), ranking by the consensus criterion separates actives from
decoys with **AUC ≥ 0.70 AND EF@1% ≥ 5.0x** (mode C). Same bar as run 2; not lowered.

## Analysis plan

    Primary (gates):  mode C consensus — min N-Fe over the 5 seeds
    Reported alongside (does NOT gate): mode C baseline — the seed-42 single search
    Metrics: AUC, EF@1/5/10%, BEDROC(alpha=20) from rdkit.ML.Scoring.
    Any active with no usable pose in any seed is ranked LAST, never dropped.

Graded by `scripts/validate_gate_consensus.py`, which re-checks this document's SHA-256 and
the frozen protocol hash before reporting.

## Pre-committed interpretation — fixed in advance

- **PASS.** The variance-reduced criterion holds on the trusted wild-type gate. It is then
  a validated estimator, and applying it to a resistant target (e.g. re-testing Y132F)
  becomes justified — under its **own** separate pre-registration. This does **not**
  retroactively overturn the reported Y132F single-search FAIL, which stands as graded under
  its own pre-registration.
- **FAIL.** Consensus does not beat the pre-registered bar on the wild-type set. Then the
  single-search criterion was not meaningfully improved by this change, the Y132F failure is
  not attributable to seed noise in a way consensus can fix, and no consensus result may be
  reported as a pass anywhere.

Whether consensus scores higher or lower than the seed-42 baseline is reported as the
outcome either way; no post-hoc aggregator (mean, k-th pose, thresholds) invented after
seeing these numbers may be substituted to manufacture a pass.

## Limitations (fixed in advance)

1. **R = 5 is modest.** Five seeds reduce but do not eliminate search variance; residual
   noise remains, especially for the most flexible ligands.
2. **Min-over-seeds monotonically decreases with R.** It helps high-variance molecules
   most. Because decoys are property-matched (similar size/flexibility), this is treated as
   per-molecule variance reduction, not an active-vs-decoy bias — but it is a known effect.
3. Same rigid-receptor applicability domain (≤ 45 heavy atoms) and presumed-inactive decoys
   as run 2.
