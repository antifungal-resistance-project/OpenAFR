# Pre-registration — WIDEN the held-out active set to retire the n=7 tip (Limitation #1)

> **STATUS: FROZEN 2026-09-05, before any active/decoy label was scored.** The test-set
> sample was generated first (deterministic, label-free — its sha256 is pinned below), then
> this document is frozen by its own SHA-256. The grader
> (`scripts/validate_gate_active_power.py`) re-checks this document's hash, `protocol.yaml`,
> and the sample hash, and refuses to certify a pass if any moved. No AUC of the widened
> sample against any inactive pool has been read at freeze time.

## Why this run is required — the single most-quoted weakness

The primary held-out validation gate (Run 2, `work/RESULTS_run2_holdout.md`,
`work/PREREGISTRATION_run2.md` sha `7c89d78c…`) passed on **7 held-out azole actives vs
348 property-matched decoys**: mode C **AUC 0.794**, but with a bootstrap 95% CI of
**0.590–0.946**. Limitation #1 of both that result and the preprint (`§6.1`) states it
plainly:

> "n = 7 actives. The bootstrap 95% CI lower bound (0.590) sits **below** the 0.70 pass
> bar … a true AUC under the bar cannot be excluded. More held-out actives would tighten
> this."

This run does exactly the thing that limitation names: it re-runs the identical gate —
same receptor, same criterion, same decoy pool, same frozen protocol and seed — with the
**active side widened from 7 to a large, blind, independent sample**, so the CI becomes
tight and its lower bound can be compared against the 0.70 bar. It answers one question:
**once the sampling noise of 7 molecules is removed, does the geometry criterion's AUC
lower bound clear the project's usefulness bar?**

This is a *power* run, not a new criterion. Everything except the number of actives is
held fixed from Run 2. This is deliberately **not** the within-azole question (that is
closed for good by look #9, `work/RESULTS_ensemble_confirm.md`); the actives here are
ranked against **property-matched decoys**, the Run-2 contrast, not against azole-bearing
analogues.

## Single-variable design — everything is held fixed from Run 2 except the actives

| | |
|---|---|
| receptor | rigid **5TZ1** chain A + heme, single conformer — the product receptor, **unchanged** |
| criterion | **mode C** — `openafr.pdbqt.min_nitrogen_iron_distance` — **unchanged** |
| protocol | `protocol.yaml` unchanged — box 26 Å @ 70.61/66.28/4.18, exhaustiveness 32, 20 modes, seed 42 |
| decoys (primary) | the **same 348 property-matched decoys** as Run 2 (`data/ligands/decoys_holdout.smi`), **unchanged** |
| **actives** | **CHANGED** — a large blind sample of measured-active azoles (below), replacing the 7 held-out azoles |

Only the actives change, so any difference from Run 2's 0.794 is attributable to the
active sample size, not to method drift.

## The widened active set — built by a frozen rule, not hand-picked

Drawn from the population already built for look #9 by `scripts/make_verified_actives.py`
(`data/ligands/verified_actives.smi`, **3,609** ChEMBL *C. albicans* measured-active azoles,
sha `1ee6e45e…`; funnel and rule in `work/PREREGISTRATION_ensemble_confirm.md`). The rule,
unchanged: has-active-evidence ∧ azole-bearing ∧ ≤45 heavy atoms ∧ ≥1 aromatic N ∧
max Morgan-r2 Tanimoto < 0.70 vs all 15 project actives + the 5TZ1 co-crystal ligand. So
every molecule is genuinely distinct chemistry from the 7 Run-2 actives (no leakage) and
in the declared applicability domain.

**The frozen test set (what this run grades).** A **deterministic random sample of N = 300,
seed 7** (`make_verified_actives.py sample --n 300 --seed 7 --out
data/ligands/active_power_sample.smi`) — a seeded shuffle over the sorted population, not a
head-of-list slice; reproducible; ~43× the original n = 7 and 1.5× look #9's n = 200, for a
tighter CI. **Seed 7, not 42**, so the sample is independent of look #9's seed-42 draw
rather than a re-grade of the same molecules under a new contrast: only **18 of 300** (6%)
compounds coincide with look #9's seed-42 n=200 draw (expected ~17 from the shared 3,609
population), and look #9 graded those against azole-bearing inactives, a different contrast.

> **PINNED:** `data/ligands/active_power_sample.smi` sha256
> `ba8d1f3cb05717017ababa321c8a74b4e0122c5eae94c91ceab92915dd3c37ae` (300 lines, seed 7).
> The grader verifies this hash before scoring.

## Primary hypothesis (H1) — the only thing that gates

Ranking the pooled molecules (the N = 300 azole-active sample + the 348 property-matched
decoys) by **mode C on rigid 5TZ1** separates measured-active azoles from property-matched
decoys at **AUC ≥ 0.70** (`protocol.yaml` `min_auc`, unchanged), AND the **bootstrap 95% CI
lower bound is reported**. The retirement of Limitation #1 is graded on the **lower bound**,
not merely the point estimate — the whole reason this run exists is that at n = 7 the lower
bound (0.590) sat below the bar. Unranked actives (docking failures) are ranked **last**,
never dropped (the Run-2 anti-inflation rule).

EF is **not** part of the gate. Run 2's EF@1% = 12.68× does not survive 5-seed resampling
(→ 0.00×, preprint §3.2); EF is reported for continuity only and gates nothing here.

## Secondary — all pre-specified, none gating

1. **Product-claim contrast (the headline AUC 0.810 number).** The same N = 300 actives vs
   the **non-azole** verified inactives (`data/ligands/verified_inactives.smi`, non-azole
   subgroup) — the novel-chemotype triage claim. Comparators: 0.810 (single-conformer look #6)
   / 0.823 (look #9 guard, n = 200). Reported with its bootstrap 95% CI. This well-powers the
   *product* AUC as well as the *gate* AUC.
2. **Statistical support**: permutation p (20,000 shuffles), bootstrap 95% AUC CI, the tip
   (best decoy above the best active; percentile-rank form).
3. **EF@1/5/10%, BEDROC(20)** — reported for continuity, non-gating (seed-fragile per §3.2).
4. **Domain-edge check**: AUC vs heavy-atom count, to confirm the widened sample is not
   carried by a size artifact.

## Pre-committed interpretation

- **RETIRED — H1 AUC ≥ 0.70 AND CI lower bound ≥ 0.70.** Limitation #1 is discharged: with
  the sampling noise of n = 7 removed, the geometry criterion's held-out AUC clears the
  usefulness bar with its lower bound, not just its point estimate. The preprint §6.1 is
  updated from "a true AUC under the bar cannot be excluded" to the measured tight interval.
  The Run-2 point estimate (0.794) stands; this run reports whether it is robustly above 0.70.
- **MARGINAL — point ≥ 0.70 but CI lower bound < 0.70.** The gate still passes on the point
  estimate but a sub-bar true AUC cannot be excluded even at N = 300. Reported honestly as a
  narrowed-but-not-cleared interval; Limitation #1 is *quantified*, not retired.
- **FAIL — point estimate < 0.70.** The Run-2 0.794 does not hold on a large blind sample;
  the primary gate's point estimate was itself an n = 7 fluctuation. This would be a material
  negative result about the headline claim and would be reported as prominently as the pass
  branches — the product scope would need re-statement, not just a footnote.

## Stopping rule (binding)

**One sample, one look.** The test set is the frozen N = 300, seed 7 sample. It will **not**
be re-sampled with a different seed, N, or novelty cutoff and re-graded if it disappoints.
Receptor, criterion, decoy pool, and protocol are all unchanged from Run 2. Any further
active-power work is a new, separately pre-registered dataset.

## Limitations (fixed in advance)

1. **Whole-cell "active" ≠ CYP51 "active"** — inherited from the look #9 active construction;
   conservative (a mislabelled active ranks low, depressing AUC).
2. **Property-matched decoys are presumed-inactive, not measured** — unchanged from Run 2;
   this run widens the *active* side only. The measured-inactive contrast is the look #6/#9
   line and is reported as secondary S1, not re-litigated here.
3. **Single active record suffices for has-active-evidence** — the sample inherits it;
   conservative and representative.
4. **Rigid single conformer / ≤45 heavy-atom domain** — unchanged; the long-tail
   posaconazole/itraconazole class remains the declared blind spot.
5. **A rigid-docking false negative ranks LAST, never dropped** — conservative for the gate.

## Reproduce (docking runs on a Linux/x86 cloud VM with Vina; prep + grade are local)

```
conda activate openafr
python scripts/make_verified_actives.py build                              # -> verified_actives.smi (3609)
python scripts/make_verified_actives.py sample --n 300 --seed 7 \
       --out data/ligands/active_power_sample.smi                          # -> frozen test set (sha pinned above)
python scripts/prep_ligands.py data/ligands/active_power_sample.smi work/active_power_ligands
RECEPTOR=work/receptor.pdbqt scripts/screen.sh work/active_power_ligands work/screen_active_power
python scripts/validate_gate_active_power.py                               # exits 0 pass / 1 fail
```
