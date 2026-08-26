# Ensemble confirmation — the within-azole lift did NOT replicate (look #9)

Date: 2026-08-26
Pre-registration: `work/PREREGISTRATION_ensemble_confirm.md` (sha256 `7754050e94523a9c…`,
verified unmodified by the grader at run time)
Test set: an **independent** measured-active azole sample — N=200, seed 42, drawn deterministically
from 3,609 ChEMBL *C. albicans* measured-active azoles (`verified_actives_sample.smi`, sha
`078c247f…`), novelty-filtered <0.70 vs all 15 project actives + the 5TZ1 co-crystal ligand — vs
the **same frozen 162 azole-bearing verified inactives** as looks #6/#7/#8.
Everything else identical to look #8: exhaustive 3-conformer ensemble (5TZ1+5FSA+5V5Z, hash-pinned),
mode C, ensemble-minimum, frozen protocol. Docking: 196/200 actives posed in all 3 conformers
(4 obabel/Vina failures, ranked last — conservative). ~43 min × 3 conformers on 12 cores.

## RESULT: GATE FAILED — look #8's 0.750 was a 7-molecule sampling fluctuation

**H1 fails: ensemble-min mode C AUC 0.682 < 0.70** on the independent 200-azole sample, against
look #8's **0.750** on the 7 held-out azoles. The within-azole lift does **not** generalise to a
broad, independent population of measured-active azoles.

| ranking | pool | AUC | vs bar |
|---|---|---|---|
| **ensemble-min mode C — independent n=200** | 200 actives + 162 azole-bearing inactives | **0.682** | **FAIL** |
| ensemble-min mode C — look #8 (n=7 held-out) | 7 actives + 162 inactives | 0.750 | (did not replicate) |
| mode C — single conformer 5TZ1, look #6 (n=7) | 7 + 162 | 0.650 | FAIL |

The bootstrap CI is now **tight** — 0.647–0.715 (n=200, versus look #8's 0.683–0.825 at n=7) —
which is the entire reason this run was required: at n=7 the interval was too wide to distinguish
a real lift from noise. It is now narrow, and its **point estimate sits below the bar**. Per the
frozen analysis plan the gate is the point estimate, so this is a FAIL.

## The signal is real but sub-threshold — state both halves honestly

This is **not** "the ensemble does nothing." The permutation test is decisive (p < 0.0001,
20,000 shuffles): the ranking is far better than chance. And the ensemble *does* add a small,
real lift over the single conformer:

    per-conformer AUC (independent n=200):  5TZ1  0.638
                                            5FSA  0.666
                                            5V5Z  0.701
    marginal-conformer increment (min):     {5TZ1}            0.638
                                            {5TZ1, 5FSA}      0.676
                                            {5TZ1,5FSA,5V5Z}  0.682

So the exhaustive ensemble lifts within-azole ranking from 0.638 (5TZ1 alone) to 0.682 — a
genuine ~+0.04 — but it lands **below the 0.70 usefulness bar the project holds itself to**.
"Statistically significant" and "clears the bar" are different claims; this run has the first and
not the second.

**And the conformer that helped is itself sample-dependent.** In look #8 the lift rode on 5FSA
(posaconazole conformer, 0.783 alone). On the independent set the best single conformer is **5V5Z**
(itraconazole, 0.701) and 5FSA is middling (0.666). The "which conformer accommodates the tails"
story was partly a lucky interaction between those particular 7 held-out azoles and the 5FSA
pocket — exactly the kind of small-sample artefact an independent confirmation exists to expose.

## The novel-chemotype guard still holds — the product thesis is intact

    novel-chemotype guard (200 actives vs NON-azole inactives): AUC 0.823   (comparator 0.810)

The ensemble does **not** break what works: separating azole actives from *non-azole* chemotypes
remains strong (0.823, above the 0.810 single-conformer comparator). The failure is specifically
the *within-class* discrimination, on an independent set — the same ceiling looks #6 and #7 hit,
now confirmed to survive even the exhaustive experimental conformer ensemble.

EF/BEDROC are weak here by construction (near-balanced pool, 196 actives vs 162 inactives — the
pre-registered Limitation 4): EF@1/5/10% ≈ 1.3×, BEDROC 0.701 are not informative in either
direction; AUC is the honest, prevalence-neutral read, which is why it alone gates.

## Pre-committed interpretation (from the pre-registration, verbatim branch)

> **FAIL — AUC < 0.70.** The within-azole lift does **not** replicate. Look #8's 0.750 was
> consistent with a 7-molecule sampling fluctuation (its CI straddled the bar). **No within-class
> claim is earned**; the defensible product remains novel-chemotype triage (AUC 0.810), and the
> within-azole line closes for good on the rigid *and* the exhaustive-ensemble receptor.

## What this settles

- **The within-azole question is now closed on both receptors.** Rigid single conformer: mode C
  0.650, mode F 0.590 (looks #6/#7). Exhaustive experimental ensemble: 0.682 on an independent,
  well-powered set (look #9). Pose geometry over a fixed *C. albicans* CYP51 target — rigid or
  the full crystallographic ensemble — cannot rank working azoles above failed azole analogues to
  the project's 0.70 bar. This is a property of the question, not of a single feature or a single
  conformer.
- **Look #8 stands as graded but is superseded in interpretation.** Its 0.750 was honestly
  reported *with* the n=7 wide-CI caveat and the explicit demand for this confirmation; the
  confirmation is the answer, and it is negative. The value of look #8 was to make this test
  well-posed, not to establish a within-class result.
- **The product is unchanged and unharmed.** Novel-chemotype triage holds at 0.810–0.823 across
  rigid and ensemble receptors. The candidate shortlist (`RESULTS_candidate_shortlist.md`) stands
  on that footing; there is no within-class re-ranking to add to it.
- **Next within-class attempt, if ever, needs a different axis** — not another conformer set. A
  flexible/induced-fit *scoring* model, or an orthogonal experimental readout, would be a new,
  separately pre-registered line. This one is closed.

## Reproduce

```
conda activate openafr
python scripts/make_verified_inactives.py fetch
python scripts/make_verified_actives.py build
python scripts/make_verified_actives.py sample --n 200 --seed 42
python scripts/prep_ligands.py data/ligands/verified_actives_sample.smi work/verified_ligands
RECEPTOR=work/receptor.pdbqt      scripts/screen.sh work/verified_ligands work/screen_verified
RECEPTOR=work/receptor_5FSA.pdbqt scripts/screen.sh work/verified_ligands work/screen_5FSA
RECEPTOR=work/receptor_5V5Z.pdbqt scripts/screen.sh work/verified_ligands work/screen_5V5Z
python scripts/validate_gate_ensemble_confirm.py    # exits 1 (FAIL)
```
