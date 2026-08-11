# Results — reliability criterion, and the geometry ceiling at the top

Date: 2026-08-11. Pre-registered in
[PREREGISTRATION_reliability.md](PREREGISTRATION_reliability.md) (sha256 f8aa1ff9…, verified
unmodified at grading time). Protocol hash 6eca16b1… verified unmodified. Reuses the
already-docked wild-type 5-seed poses (no new docking).

> **Outcome: FAIL — and it settles the question.** A pool-size-neutral robust estimator
> *also* cannot recover top-1% separation. Per the pre-registration this was the second and
> final aggregator on these poses; the committed conclusion is that iron-approach geometry
> has reached its ceiling at the extreme top on this decoy set. No third aggregator.

## The three aggregators, side by side (same wild-type poses)

| aggregator | AUC | **EF@1%** | EF@5% | EF@10% | BEDROC |
|---|---:|---:|---:|---:|---:|
| single search (seed 42, run-2 baseline) | 0.794 | **12.68x** | 5.63x | 4.23x | 0.325 |
| consensus — min over 5 seeds | 0.853 | **0.00x** | 8.45x | 5.63x | 0.389 |
| reliability — median over 5 seeds | 0.830 | **0.00x** | 5.63x | 4.23x | 0.310 |

Pass bar: AUC ≥ 0.70 AND EF@1% ≥ 5.0x. Reliability clears AUC but EF@1% = 0.00x → FAIL.

## What this settles

The two facts that matter, read across the table:

1. **Broad enrichment is real and robust.** AUC stays ~0.79–0.85 and EF@5% ~5.6–8.5x under
   every aggregator — and consensus/reliability actually *improve* the broad measures. The
   tool genuinely concentrates real azoles in the top ~5–10%. That holds up.

2. **The dramatic top-1% enrichment was an artifact of docking once.** The single-search
   EF@1% of 12.68x does **not** survive more thorough sampling: run five searches and combine
   them *any* way — take the best (min) or the typical (median) — and EF@1% collapses to zero.
   The reason is not noise that averaging can fix; it is that with more search the
   property-matched decoys **reliably reveal genuine iron-coordinating poses** (see
   `work/views/` — a top decoy swings an aromatic nitrogen to the iron as well as a real
   drug does). At the very top, azole-like look-alikes coordinate the iron about as well as
   azoles. Geometry alone cannot separate them there.

Reliability was the clean test that distinguished "fixable noise" from "real ceiling": a
pool-size-neutral median removed the extreme-value bias that sank consensus, and top-1%
enrichment *still* did not return. So it is the ceiling, not the estimator.

## Honest consequence for the project's headline

The run-2 headline (AUC 0.794 / **EF@1% 12.68x**) should now be read with a caveat: the AUC
is robust, but the **EF@1% figure was inflated by single-shot docking** and is not a reliable
property of the criterion. The defensible claim is *broad* enrichment (a trustworthy top-~5%
shortlist), **not** razor-top ranking. VISION.md's summary line warrants this qualification.
This does not touch the Y132F or consensus results, which stand as graded.

## Multiple-comparison disclosure (as pre-committed)

This is the **second aggregator** tested on these wild-type poses (after consensus). The
pre-registration bounded the search to exactly two and pre-committed the "ceiling" conclusion
on a fail, precisely to prevent trying estimators until one passes. Both robust aggregators
failed EF@1%; no third will be tried on these poses. A further attempt would need a
different *kind* of signal (an orthogonal feature), not another way to average N-Fe.

## Held-out actives by reliability (median N-Fe, Å)

    efinaconazole 2.735   luliconazole 2.781   fenticonazole 2.915   bifonazole 3.084
    sertaconazole 3.187   oxiconazole 3.242    butoconazole 5.441

Real azoles reliably reach ~2.7–3.2 Å — genuinely coordinating — but so do enough decoys to
own the top 1%. (butoconazole remains the consistent weak active across every method.)

## Limitations

1. R = 5 gives a 5-point median (coarse); still, both robust aggregators agree.
2. One decoy set (property-matched, presumed inactive). The ceiling is demonstrated *on this
   set*; a different decoy construction could differ, but property-matched decoys are the
   fair test.
3. Rigid receptor, ≤ 45 heavy-atom applicability domain, as run 2.

## Reproduce

```
conda activate openafr
# reuses work/screen_consensus_wt/ from the consensus run (5 seeds, no new docking):
python scripts/validate_gate_reliability.py work/screen_consensus_wt   # exits 0 pass / 1 fail
```
