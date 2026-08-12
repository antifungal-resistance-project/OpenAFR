# Results — orthogonal axial-geometry feature, wild-type held-out gate

Date: 2026-08-11. Pre-registered in
[PREREGISTRATION_axial.md](PREREGISTRATION_axial.md) (sha256 `3ff59c29…`, verified unmodified
at grading time). Protocol hash `6eca16b1…` verified unmodified. Graded by
`scripts/validate_gate_axial.py` on the reused wild-type 5-seed poses
(`work/screen_consensus_wt/`, seeds {42,1,2,3,4}) — **no new docking**; same 7 actives, 348
decoys, rigid 5TZ1 receptor.

> **Outcome: the gate FAILED on the pre-registered bar (EF@1% = 0.00x).** Adding an
> orthogonal *direction* signal to the iron-approach distance did not recover top-1%
> separation. Per the pre-registration this is reported as the outcome, unrevised, and it
> triggers the **stronger** committed conclusion: on this decoy set the top-1% ceiling is not
> merely a distance limitation — the property-matched decoys coordinate the iron with
> **azole-like axial geometry too**, so pose geometry (distance + direction) cannot own the
> extreme top here.

## What was tested (the first feature orthogonal to N–Fe distance)

Azole inhibition is **axial**: the drug nitrogen approaches the heme iron along the porphyrin
normal, not from the side. Every prior look ranked on the same one number — *how close* the
nitrogen gets — and the reliability pre-registration closed the "aggregate that distance
differently" line for good. This is the one authorized next step: a genuinely different
question — *from what direction*.

Per pose, on the **same** nitrogen the distance criterion uses (nearest the iron):

    v = N − Fe               |v| = the N–Fe distance
    θ = acute angle between v and the heme normal n   (folded to [0°, 90°])
    S = |v| · (1 + sin θ)     mode E — rewards close AND on-axis; λ = 1, untuned

Heme normal `n = normalize((NC − NA) × (ND − NB))` from `HEM A 601`, fixed once from the
receptor (`n = (+0.6048, −0.7591, −0.2407)`). Aggregation carries the reliability lesson:
per molecule take the per-seed **minimum S**, then rank by the **median** over 5 seeds
(pool-size-neutral, unlike `min`). Only the ranking key changed.

## The gate result

| ranking key | AUC | **EF@1%** (gates) | EF@5% | EF@10% | BEDROC(α=20) | top-15 |
|---|---:|---:|---:|---:|---:|---:|
| run-2 baseline — N–Fe, single search | 0.794 | **12.68x** | 5.63x | 4.23x | 0.325 | 2/7 |
| consensus — N–Fe, min over 5 seeds | 0.853 | **0.00x** | 8.45x | 5.63x | 0.389 | 3/7 |
| reliability — N–Fe, median over 5 seeds | 0.830 | **0.00x** | 5.63x | 4.23x | 0.310 | 2/7 |
| **axial (mode E) — S, median over 5 seeds** | **0.862** | **0.00x** | **8.45x** | **5.63x** | **0.340** | **2/7** |

Pass bar (both required): **AUC ≥ 0.70 AND EF@1% ≥ 5.0x**. AUC clears it and is the *best of
any look so far* (0.862); **EF@1% is 0.00x**. Gate = **FAIL**.

Reported alongside (do **not** gate):

| slice | AUC | EF@1% | EF@5% | BEDROC |
|---|---:|---:|---:|---:|
| axial angle **alone** — median per-seed min θ (direction only) | 0.861 | 0.00x | 5.63x | 0.240 |
| reliability distance — median per-seed min N–Fe (reference) | 0.830 | 0.00x | 5.63x | 0.310 |
| seed-42 **single-search** S — does the single-shot top-1% survive? | 0.819 | 0.00x | 8.45x | 0.320 |

**Harness check:** the reliability-distance slice reproduces the published reliability numbers
(0.830 / 0.00x / 5.63x / 4.23x / 0.310) bit-for-bit, confirming the regenerated poses are
consistent with the earlier run.

**The single-shot top-1% does not survive.** The pre-registration asked specifically whether
adding the orthogonal term to the identical seed-42 poses preserves the 12.68x the single
search once had. It does not — **seed-42 single-search S collapses to EF@1% 0.00x** as well.
The direction term, added on the same poses, pushes decoys up rather than actives.

## Why the top-1% still collapsed (concretely)

Top-1% of 355 molecules is the top ~3–4. Those slots go entirely to decoys — and not to
decoys sneaking in equatorially, but to decoys coordinating the iron **more axially than the
real drugs**:

| rank | molecule | median S | median θ | median N–Fe | |
|---:|---|---:|---:|---:|---|
| 1 | decoy_CHEMBL6328_for_efinaconazole | 2.881 | 3.4° | 2.702 | decoy |
| 2 | decoy_CHEMBL275096_for_oxiconazole | 2.928 | 2.7° | 2.790 | decoy |
| 3 | decoy_CHEMBL10329_for_oxiconazole | 3.017 | 3.9° | 2.793 | decoy |
| 4 | decoy_CHEMBL447532_for_butoconazole | 3.042 | 7.5° | 2.691 | decoy |
| 5 | decoy_CHEMBL10330_for_sertaconazole | 3.079 | 4.3° | 2.842 | decoy |
| 6 | decoy_CHEMBL10547_for_oxiconazole | 3.094 | 4.3° | 2.827 | decoy |
| **7** | **luliconazole** | **3.133** | **6.2°** | **2.781** | **ACTIVE** |

Six decoys beat the best real drug. Their axial angles (2.7°–7.5°) are as small as or
*smaller* than luliconazole's (6.2°) — they are not off-axis impostors, they are doing
azole-like axial coordination. Adding `sin θ` therefore cannot demote them; it slightly
demotes the actives whose best-S pose is a touch less axial. The seven held-out actives land
at ranks 7, 10, 18, 33, 41, 53, 201 — a *broadly* enriched list (hence AUC 0.862, EF@5%
8.45x) with a decoy-owned tip.

## Interpretation (pre-committed, FAIL branch)

Quoting the pre-registration verbatim in spirit: *"adding the axial-direction signal still does
not separate real azoles from property-matched decoys at the top 1%. The committed conclusion
is then stronger than reliability's: on this decoy set, the top-1% ceiling is not just a
distance limitation — the decoys coordinate the iron with azole-like axial geometry too, so
pose geometry (distance + direction) cannot own the extreme top here. The defensible product
remains a trustworthy top-~5% shortlist."*

Therefore:

- The axial criterion is **not** a validated top-ranking rescue and is **not** carried forward
  to the resistant target. The graded run-2, consensus, reliability, and Y132F results stand
  unchanged; this does not retroactively alter any of them.
- **No second orthogonal feature will be fished on these poses.** The pre-registration bounded
  this to exactly one feature for a reason: trying feature after feature until one passes is the
  p-hacking the reliability bound was built to stop. A further attempt requires an *independent
  dataset* — a new decoy construction or new actives — not another feature on this one.

## The real lesson (honest, and worth having)

Three looks blamed the top-1% collapse on *distance* aggregation. This look adds *direction*
and shows the collapse is not about distance at all: **property-matched decoys reach the heme
iron both close and on-axis**, indistinguishably from real azoles at the very tip. That is a
sharper, more useful statement than "our aggregator was noisy." It says the *pose geometry*
this pipeline can measure — closest-nitrogen distance and its approach angle — genuinely
saturates at the top of a hard, property-matched decoy set. What the geometry criterion buys
is **broad** enrichment (AUC ~0.86, EF@5% ~8.5x), i.e. a trustworthy top-~5% shortlist to hand
a wet lab — not razor-top ranking. The [VISION](../VISION.md) headline already carries exactly
that caveat; this result is the fourth independent confirmation of it.

What could still own the extreme top is a signal this rigid-pose pipeline does **not** see:
explicit iron–nitrogen coordination chemistry (QM/semi-empirical), induced-fit pocket
response, or desolvation — none of which is decidable from a single rigid-receptor Vina pose,
and none of which is in scope here.

## Limitations (from the pre-registration, all apply)

1. **Same one dataset, a 4th time.** A new *question* (direction, not distance), authorized by
   the reliability prereg — but still one 7-active / 348-decoy set. The FAIL branch forecloses
   further looks for exactly this reason.
2. The normal `n` is estimated from 4 atoms of a rigid receptor; a real porphyrin ruffles and
   the distal pocket flexes. `n` is a fixed idealization.
3. R = 5 gives a coarse 5-point median with its own sampling noise.
4. θ is measured on the single nitrogen nearest the iron (the distance criterion's convention);
   a pose whose *second* nitrogen is the true axial one is scored on the first.
5. Same rigid-receptor applicability domain (≤ 45 heavy atoms) and presumed-inactive decoys as
   run 2.

## Multiple-comparison history (stated plainly)

This is the **4th look** at the run-2 wild-type poses — baseline, consensus, reliability,
axial — and the **1st** that adds a feature orthogonal to N–Fe distance. Of the four, one
(the single-shot baseline) passed EF@1% and three (all with more thorough sampling, or the
orthogonal feature) did not. The pre-registered stopping rule is now in force: no further
feature or aggregator is tried on this dataset.

## Reproduce

```
conda activate openafr
python scripts/prep_receptor.py                       # -> work/receptor.pdbqt
cat data/ligands/actives_holdout_final.smi data/ligands/decoys_holdout.smi > /tmp/run2_all.smi
python scripts/prep_ligands.py /tmp/run2_all.smi work/run2_ligands
RECEPTOR=work/receptor.pdbqt scripts/screen_consensus.sh work/run2_ligands work/screen_consensus_wt "42 1 2 3 4"
python scripts/validate_gate_axial.py work/screen_consensus_wt   # exits 0 pass / 1 fail
```
