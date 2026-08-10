# Results — consensus (variance-reduced) criterion, wild-type held-out gate

Date: 2026-08-10. Pre-registered in
[PREREGISTRATION_consensus.md](PREREGISTRATION_consensus.md) (sha256 1c6ea042…, verified
unmodified at grading time). Protocol hash 6eca16b1… verified unmodified.

> **Outcome: the gate FAILED on the pre-registered bar.** Consensus improved almost every
> ranking metric — but destroyed the one the gate keys on (EF@1%). Per the pre-registration
> this is reported as the outcome, unrevised; no other metric is substituted to call it a
> pass. An informative negative that exposes a real flaw in the consensus aggregator.

## What was tested

Whether estimating the iron-approach criterion with **5 independent docking searches per
molecule** (different seeds; box/exhaustiveness/modes unchanged) and ranking by the
**minimum N-Fe across seeds** ranks the run-2 wild-type held-out actives better than a
single search. Same 7 actives, same 348 decoys, same rigid 5TZ1 receptor. The seed-42 slice
is the single-search baseline (and reproduces published run-2 numbers exactly, confirming
the harness).

## The gate result

| metric | single-search baseline (seed 42) | **consensus (min over 5 seeds)** |
|---|---:|---:|
| AUC | 0.794 | **0.853** ↑ |
| **EF@1%** (gates) | **12.68x** | **0.00x** ↓↓ |
| EF@5% | 5.63x | **8.45x** ↑ |
| EF@10% | 4.23x | **5.63x** ↑ |
| BEDROC(α=20) | 0.325 | **0.389** ↑ |
| actives in top 15 | 2/7 | **3/7** ↑ |

Pass bar (both required): **AUC ≥ 0.70 AND EF@1% ≥ 5.0x**. AUC clears it easily and
*improves*; **EF@1% collapses to 0.00x**. Gate = FAIL.

## Why the top-1% collapsed while everything else improved

Consensus does help the actives — dramatically for the flexible ones (efinaconazole's best
coordination went from a poor single-search 3.966 Å to 2.701 Å). Broad enrichment (AUC,
EF@5/10%, BEDROC) all rose. But the very top of the list went to decoys:

| rank | molecule | consensus N–Fe (Å) | |
|---:|---|---:|---|
| 1 | decoy_CHEMBL7641_for_oxiconazole | 2.666 | decoy |
| 2 | decoy_CHEMBL447532_for_butoconazole | 2.673 | decoy |
| 3 | decoy_CHEMBL11793_for_butoconazole | 2.687 | decoy |
| 4 | decoy_CHEMBL6328_for_efinaconazole | 2.693 | decoy |
| 5 | decoy_CHEMBL11724_for_butoconazole | 2.698 | decoy |
| **6** | **efinaconazole** | **2.701** | **ACTIVE** |
| 7 | fenticonazole | 2.703 | ACTIVE |

Five decoys beat the best real drug — by as little as 0.035 Å. **The cause is pool-size
bias in "best-of-R" aggregation.** Taking the minimum over 5 searches is an extreme-value
statistic; run it over **348 decoys** (each getting 5 chances) versus **7 actives** and,
purely by numbers, a handful of decoys draw an unusually close pose that no active can beat
at the very tip. EF@1% looks only at the top ~3–4 molecules, so it goes to zero even as the
actives cluster strongly just below (EF@5% *rose* to 8.45x).

## Interpretation (pre-committed)

Per the pre-registration's FAIL branch: **the consensus criterion, as specified, is not
supported on the wild-type held-out gate.** Therefore:

- It is **not** a validated estimator, and it is **not** applied to the resistant target.
  The reported **Y132F single-search FAIL stands** as graded under its own pre-registration;
  consensus does not overturn it.
- No post-hoc aggregator (the pre-registration named mean, k-th pose, thresholds) invented
  after seeing these numbers is substituted to manufacture a pass. The gating metric was
  EF@1% and it failed.

## The real lesson (honest, and worth having)

The pre-registered design was chosen "on principle" (min = the achievable minimum). The
result shows that principle has a **scale flaw**: `min`-over-R is not pool-size-neutral, so
naive multi-seed consensus *inflates the extreme tail of whichever class has more members*
— here, the decoys — and damages exactly the top-1% enrichment that matters most for handing
a wet lab a short list. Consensus made the tool a **better broad ranker** (AUC, EF@5/10%,
BEDROC) and a **worse pointy-end ranker** (EF@1%). That is a specific, reusable finding
about consensus docking, not just a null result.

A pool-size-neutral estimator (e.g. a fixed per-molecule statistic that does not reward
extra draws, or a rank-based combination) is the natural next hypothesis — but only as a
**new, separately pre-registered** test, not a rescue of this one.

## Limitations (from the pre-registration, all apply)

1. R = 5 reduces but does not remove search variance.
2. `min`-over-seeds monotonically decreases with R and is pool-size-biased (the finding
   above), so it favours the larger class at the extreme tail.
3. Same rigid-receptor applicability domain and presumed-inactive decoys as run 2.

## Reproduce

```
conda activate openafr
python scripts/prep_receptor.py                       # -> work/receptor.pdbqt
cat data/ligands/actives_holdout_final.smi data/ligands/decoys_holdout.smi > /tmp/run2_all.smi
python scripts/prep_ligands.py /tmp/run2_all.smi work/run2_ligands
RECEPTOR=work/receptor.pdbqt scripts/screen_consensus.sh work/run2_ligands work/screen_consensus_wt "42 1 2 3 4"
python scripts/validate_gate_consensus.py work/screen_consensus_wt   # exits 0 on pass, 1 on fail
```
