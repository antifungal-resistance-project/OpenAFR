# Pre-registration — orthogonal feature: axial coordination geometry

**Written 2026-08-11, BEFORE computing the feature on any pose.** Everything below is fixed
in advance. Inherits the run-2 protocol and reuses the already-docked wild-type 5-seed poses
from the consensus run (`work/screen_consensus_wt/`, **no new docking**). Read
[PREREGISTRATION_reliability.md](PREREGISTRATION_reliability.md) and
[RESULTS_reliability.md](RESULTS_reliability.md) first — this test is the specific "different
kind of signal" that pre-registration named as the only legitimate next step.

## Background — exactly what this is allowed to test, and why now

Three looks at the wild-type poses are already on the record, and all three ranked on the
**same one number** — the N–Fe distance — differing only in how they aggregated it:

| look | ranking key | EF@1% |
|---|---|---:|
| run-2 baseline | N–Fe, single search (seed 42) | 12.68x |
| consensus | N–Fe, min over 5 seeds | 0.00x |
| reliability | N–Fe, median over 5 seeds | 0.00x |

The reliability pre-registration **closed the "aggregate the distance differently" line for
good** (2 aggregators, then stop) and stated the committed conclusion: iron-approach *distance*
has hit its ceiling at the extreme top on this decoy set, because property-matched decoys
reliably reach the iron too (visible in `work/views/`). It named exactly one honest way
forward: *"a different kind of signal (e.g. an orthogonal feature), not another way to average
the same distances."* This document is that — and, to avoid fishing, it is bounded to **one**
such feature, defined and combined mechanistically in advance.

## The orthogonal feature (fixed; genuinely not a function of N–Fe distance alone)

Azole inhibition is **axial** coordination: the drug nitrogen approaches the heme iron along
the porphyrin normal, trans to the proximal cysteine thiolate — not from the side. N–Fe
distance measures only *how close*; it is blind to *from what direction*. A property-matched
decoy can reach a short N–Fe distance with an **off-axis (equatorial)** nitrogen that is not
doing real axial coordination. The orthogonal signal is that **angle**.

Concretely, per pose:

1. **Heme normal `n`.** From `HEM A 601` in `work/receptor_A.pdb`, take the four pyrrole
   nitrogens NA, NB, NC, ND. Define `n = normalize( (NC − NA) × (ND − NB) )` — the normalized
   cross product of the two porphyrin diagonals. Deterministic; computed once from the fixed
   receptor.
2. **Coordinating nitrogen.** The ligand nitrogen **nearest the iron** — the *same* atom the
   N–Fe criterion already uses (`min_nitrogen_iron_distance`). Let `v = N − Fe`.
3. **Axial angle `θ`.** `θ = acute angle between v and n` (folded to [0°, 90°]). Because we
   fold to the acute angle and the feature below uses `sin θ`, the **orientation of `n` is
   irrelevant** — there is no distal-vs-proximal choice to tune, and no free sign.

**Combined criterion (mode E), fixed with a single untuned constant:**

    S = |v| · (1 + sin θ)

- `S = |v|` (the raw N–Fe distance) for a perfectly axial approach (θ = 0°), and `2·|v|` for a
  fully equatorial one (θ = 90°). It rewards *close AND on-axis*, and penalizes a nitrogen that
  is close but coming in from the side — the exact decoy failure mode the reliability views show.
- The weight is **λ = 1** (both terms in Å), pre-committed, **not** fit. No threshold, no tuned
  cutoff, no per-metric search. This is chosen for mechanistic form, not because it scores best —
  it cannot, since it has never been computed.

**Aggregation over seeds (carries the reliability lesson forward).** For each molecule compute
the per-seed **minimum S** (one value per seed, R = 5, seeds {42,1,2,3,4}), then rank by the
**median** of those 5 values (smaller = better). Median, not min, because the consensus run
proved `min`-over-R is an extreme-value statistic that hands the top to a large decoy pool;
the median is pool-size-neutral. Same poses, box, exhaustiveness, modes, receptor as run 2 —
**only the ranking key changes.**

## Primary hypothesis (H1-axial)

On the run-2 wild-type held-out set (same 7 actives, same decoys, `work/receptor.pdbqt`),
ranking by the median-of-per-seed-minimum **S** separates actives from decoys with
**AUC ≥ 0.70 AND EF@1% ≥ 5.0x** (mode E). Same bar as run 2; not lowered.

## Analysis plan

    Primary (gates):  mode E — median over 5 seeds of the per-seed minimum S = |v|·(1+sinθ)
    Reported alongside (do NOT gate):
      - pure axial angle: median of per-seed minimum θ (is direction alone informative?)
      - reliability distance (median N–Fe) — already graded, for reference
      - seed-42 single-search S — to see if the top-1% that single-shot distance had
        survives when the orthogonal term is added on the identical poses
    Metrics: AUC, EF@1/5/10%, BEDROC(alpha=20) from rdkit.ML.Scoring.
    Any active whose nearest-N or angle is undefined in every seed is ranked LAST, never dropped.

Graded by `scripts/validate_gate_axial.py`, which re-checks this document's SHA-256 and the
frozen protocol hash before reporting, and refuses to certify a pass if either was edited after
freezing.

## Pre-committed interpretation — and a hard bound on forking paths

This is the **first and only orthogonal feature** tested on these wild-type poses. Trying
feature after feature until one passes would be exactly the p-hacking the reliability bound was
built to stop, so it is bounded here too:

- **PASS.** An orthogonal, mechanistically-motivated axial-geometry feature recovers top-1%
  separation that distance-aggregation alone could not. Then the tool's razor-top ranking is
  rescuable — but by *geometry of coordination*, not distance — and this criterion becomes a
  candidate to carry forward, validated **only** under its own separate pre-registration
  elsewhere (Y132F included). It does **not** retroactively change the run-2, consensus, or
  reliability results, which stand as graded.
- **FAIL.** Adding the axial-direction signal *still* does not separate real azoles from
  property-matched decoys at the top 1%. The committed conclusion is then stronger than
  reliability's: **on this decoy set, the top-1% ceiling is not just a distance limitation — the
  decoys coordinate the iron with azole-like axial geometry too, so pose geometry (distance +
  direction) cannot own the extreme top here.** The defensible product remains a trustworthy
  top-~5% shortlist. **No second orthogonal feature will be fished on these poses;** a further
  attempt would require an independent dataset (a new decoy construction or new actives), not
  another feature on this one.

Whichever occurs is reported as the outcome. No statistic invented after seeing these numbers
is substituted to manufacture a pass. The multiple-comparison history is stated plainly in the
results: this is the **4th look** at the run-2 wild-type poses (baseline, consensus, reliability,
axial), and the **1st** that adds a feature orthogonal to N–Fe distance.

## Limitations (fixed in advance)

1. **Same one dataset, a 4th time.** This reuses the run-2 wild-type poses yet again. It is a
   genuinely new *question* (direction, not distance), and the reliability prereg authorized
   exactly this step — but it is still another look at one 7-active/348-decoy set, and the FAIL
   branch above forecloses further looks precisely for that reason.
2. **The normal is estimated from 4 atoms** of a rigid receptor; a real porphyrin ruffles and
   the distal pocket flexes. `n` is a fixed idealization, not the instantaneous plane.
3. **R = 5 gives a 5-point median** — coarse, with its own sampling noise (as in reliability).
4. **Nearest-N choice.** θ is measured on the single nitrogen closest to the iron, matching the
   distance criterion; a pose whose *second* nitrogen is the true axial one is scored on the
   first. This is the same atom-selection convention run 2 used, kept for consistency.
5. Same rigid-receptor applicability domain (≤ 45 heavy atoms) and presumed-inactive decoys as
   run 2.

## Reproduce (after grading)

```
conda activate openafr
# reuses work/screen_consensus_wt/ (5 seeds, no new docking):
python scripts/validate_gate_axial.py work/screen_consensus_wt   # exits 0 pass / 1 fail
```
