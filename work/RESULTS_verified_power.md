# Verified-inactive gate at power — Limitation #2 is NOT retired: the pooled AUC falls to 0.688 at N=300

Date: 2026-09-08
Pre-registration: `work/PREREGISTRATION_verified_power.md` (sha256 `40b656a38767d9ab…`,
verified unmodified by the grader at run time; hash file `work/PREREG_verified_power.sha256`).
Frozen inputs, all re-checked by the grader before scoring:
`data/ligands/active_power_sample.smi` (sha256 `ba8d1f3cb0571701…`, N=300 measured-active azoles,
seed 7 — the same sample look #10 used), `data/ligands/verified_inactives.smi`
(sha256 `37a20eec1c026c42…`, 279 compounds MEASURED not to inhibit *C. albicans*),
`protocol.yaml` (sha256 `6eca16b1be2e5e33…`).
Compute: local Apple Silicon, `openafr` env (vina 1.2.7, obabel 3.1.0, rdkit 2025.09.5), 12 cores,
~2 h wall for 569 dockings. Grader: `scripts/validate_gate_verified_power.py`.

Single-variable widening of look #6 (`work/RESULTS_verified_inactives.md`): rigid 5TZ1, mode C
(minimum nitrogen-to-heme-iron distance), the same 279 **measured** inactives, the frozen
protocol/seed — **only the actives changed**, from 7 held-out azoles to the blind N=300 sample.
The mirror of look #10 (`work/RESULTS_active_power.md`), which retired Limitation #1 by widening
the same active side against *presumed*-inactive decoys. Here the contrast is against *measured*
non-inhibitors — the preprint's self-described "single largest caveat" (Limitation #2).

Docking: **292/300 actives and 277/279 inactives posed** (8 actives + 2 inactives were obabel/Vina
failures, ranked LAST per pre-registration — never dropped). The 8 failed actives are the *same
molecules* look #10 lost on GCP, a clean reproducibility signal.

## RESULT: GATE FAILS on the primary — and the failure is the honest, informative one

**H1 does NOT pass: mode C AUC 0.688 < 0.70**, bootstrap 95% CI **0.656–0.720**, permutation
p < 0.0001. The point estimate sits below the usefulness bar and the CI *straddles* it. Per the
pre-registered FAIL branch this is reported as prominently as a pass would have been:
**look #6's 0.716 was an n=7 fluctuation.** Once the sampling noise of seven actives is removed,
the criterion does **not** separate real azoles from compounds measured not to inhibit *to the
0.70 bar*.

| contrast (mode C) | pool | AUC | 95% CI | comparator |
|---|---|---|---|---|
| **PRIMARY — actives vs ALL measured inactives** | 292 act + 277 inact | **0.688** | **0.656 – 0.720** | look #6 (n=7) 0.716, CI 0.506–0.885 |
| S1 — vs **non-azole** measured inactives | 292 + 115 | **0.782** | — | look #6 0.810 · **look #10 0.774** |
| S2 — vs **azole-bearing** measured inactives | 292 + 162 | **0.622** | — | look #6 0.650 |
| mode A (Vina score) — reference | 292 + 277 | 0.572 | — | — |

The signal is **real** (p < 0.0001) but **sub-threshold**. This is not MARGINAL in the
pre-registered sense (that required point ≥ 0.70 with CI lower bound < 0.70); the point estimate
itself is below the bar, so it is a FAIL.

## Why it fails — the S1/S2 split reasserts itself at power (this is the finding)

The pooled 0.688 is not a mystery; it is the exact chemotype-vs-potency decomposition looks #6
and #9 already found, now nailed down with N=300 instead of 7:

1. **The product claim HOLDS and replicates.** Against **non-azole** measured inactives — the
   actual novel-chemotype triage claim — the criterion scores **AUC 0.782** (S1). That reproduces
   look #10's GCP secondary (0.774) almost exactly on freshly, locally docked poses — a clean
   GCP→local reproducibility check — and is consistent with look #6's 0.810 and look #9's 0.823
   guard. **The product axis is powered and intact.**
2. **The within-warhead ceiling is CONFIRMED real, not n=7 noise.** Against **azole-bearing**
   measured inactives — real working azoles vs *failed* azole analogues — the criterion scores
   **AUC 0.622** (S2), consistent with look #6's 0.650 and below the bar. This is the S2 branch
   the pre-registration anticipated: at power the within-azole ceiling is a property of the
   chemistry, echoing look #9's independent closure of within-class ranking. Per the stopping rule
   this does NOT reopen the within-class question.
3. **The pooled primary lands between them, dragged under the bar by benchmark composition.** The
   279 measured inactives are **58% azole-bearing (163 of 279)** — because the only compounds
   anyone measures against *C. albicans* are antifungal med-chem programs. So the pooled contrast
   is dominated by the hard S2 half, and the pooled AUC (0.688) sits between S1 (0.782) and S2
   (0.622), just below 0.70.

**The honest one-line reading:** against *measured* non-inhibitors the criterion separates real
azoles from **non-azole** chemotypes (0.782) but **cannot** separate them from **failed azole**
analogues (0.622); the pooled number is below the bar and Limitation #2 is **quantified at power,
not retired.**

## Secondaries

- **S3 — the tip.** Best true active `active_CHEMBL365796` now ranks **1 / 569** (0 inactives above
  it), versus look #6 where the best active was 8th (7 measured inactives above it). More actives
  sampled put a real drug at the very top, but the *bulk* separation is still 0.688 — the tip
  improving while the body stays sub-bar is consistent with a chemotype signal, not a potency one.
  Validated iron-bound band (2.47–2.88 Å): 134 molecules, **97 actives / 37 measured inactives** —
  at N=300 the band is now active-majority (look #6's band, on 7 actives, was inactive-dominated).
- **S4 — no size artifact.** Actives split at the median 32 heavy atoms: small half AUC 0.702,
  large half 0.672. The pooled sub-bar result is not a molecular-weight effect; both halves bracket
  it.
- **EF** is weak and non-gating by construction (near-balanced pool): primary EF@1% 1.62×.

## What this establishes, and what it does not

**Does establish:**
1. **Limitation #2 cannot be retired.** Look #6's pooled 0.716 does not survive powering; the
   powered pooled AUC against measured inactives is 0.688 [0.656, 0.720]. Any claim that the
   criterion clears the bar against *measured* non-inhibitors must be withdrawn.
2. **The product claim survives powering, on the axis it was always scoped to.** Novel-chemotype
   triage (vs non-azole measured inactives) is AUC 0.782, replicating look #10 GCP→local. The
   defensible product is unchanged and now has a second powered measurement behind it.
3. **The chemotype-not-potency mechanism is confirmed at N=300.** 0.782 (non-azole) vs 0.622
   (failed azoles) is the same split look #6 saw at n=7 and look #9 closed on an independent active
   set — now with tight support (p < 0.0001).

**Does NOT establish:**
1. Any within-azole ranking ability (S2 0.622 < 0.70; the ceiling is real and closed, look #9).
2. Anything about resistance-breaking (§4.4 untouched).

## Consequence for the preprint (§6 Limitation 2, §4.1)

Limitation #2 is **not** discharged the way Limitation #1 was. The update is a *quantification*,
not a retirement: replace "the single largest caveat … presumed inactive" with the measured
powered numbers — **pooled AUC 0.688 [0.656, 0.720] vs 279 measured non-inhibitors, sub-bar; the
separation is a chemotype axis (0.782 vs non-azole) that does not extend to ranking within the
warhead (0.622 vs failed azoles)** — and correct the standing citation of look #6's 0.716 as a
pass against measured inactives. This *strengthens* the paper's core thesis (a chemotype ceiling,
not a potency tool) while removing an over-optimistic small-sample number.

## Reproduce

```bash
conda activate openafr
python scripts/prep_receptor.py                                              # -> work/receptor.pdbqt
python scripts/prep_ligands.py data/ligands/active_power_sample.smi work/verified_power_ligands
python scripts/prep_ligands.py data/ligands/verified_inactives.smi  work/verified_power_ligands
RECEPTOR=work/receptor.pdbqt scripts/screen.sh work/verified_power_ligands work/screen_verified_power 12
python scripts/validate_gate_verified_power.py                               # exits 1 (FAIL), AUC 0.688
```
