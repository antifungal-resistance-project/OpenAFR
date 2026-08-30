# External holdout — pre-registered generalization test (issue #82)

Date: 2026-08-30
Pre-registration: `work/PREREGISTRATION_external.md`
(sha256 `66df34cb...`, verified unmodified by the grading script at run time)
Protocol: `protocol.yaml` (sha256 `6eca16b1...`) — box 26 Å @ 70.61/66.28/4.18,
exhaustiveness 32, num_modes 20, seed 42, rigid 5TZ1 chain A + heme. No per-molecule tuning.
Test set: **50 never-scored *C. albicans* azole actives + 335 property-matched decoys**
(385 total). 6 molecules produced no dockable pose (2 actives, 4 decoys — a persistent
obabel→PDBQT conversion error on those structures) and are ranked **LAST**, never dropped.

## RESULT: SUPPORTED-BUT-WIDE — the criterion generalizes; point estimate clears the bar

Pre-registered bar: **mode-C AUC ≥ 0.70**, qualified by the bootstrap 95% CI per the
three-outcome interpretation fixed in advance.

| mode | AUC | EF@1% | EF@5% | BEDROC | role |
|---|---|---|---|---|---|
| **C — iron-approach distance** | **0.715** | 7.70x | **4.62x** | 0.553 | PRIMARY — point estimate passes |
| A — Vina score alone | 0.590 | 0.00x | 0.00x | 0.083 | secondary — far worse |

The geometric criterion, **frozen** and applied with no re-tuning, ranks 50 azoles it has
never seen above 335 property-matched decoys at AUC 0.715 — above the 0.70 bar, with 11 of
the top 15 molecules being actives. Overfitting of the criterion to the run-2 development
set is **not** what produced the run-2 result: it transfers to a disjoint, never-scored set.

## Statistical support (mode C)

    permutation test (20,000 label shuffles): p < 0.0001   (0 of 20,000 shuffles matched)
    bootstrap 95% CI over actives:            AUC 0.630 - 0.795
    active rank positions (of 385):           1,2,3,4,5,7,9,10,11,13,15,17,25,29,30,31,39,
                                              44,50,56,61,64,71,79,81,92,93,95,100,102,115,
                                              128,135,137,155,158,167,190,213,249,255,272,
                                              283,305,308,309,319,350,380,381

The permutation test rules out chance decisively (p < 0.0001). The bootstrap 95% CI
(0.630–0.795) **includes values below the 0.70 bar**, so — exactly as pre-registered — this
is scored **SUPPORTED-BUT-WIDE**, not a clean PASS: the point estimate and the permutation
test support generalization, but a true AUC just under the bar cannot be excluded at n = 50.

## Reading this honestly

- **It confirms generalization, and it is better powered than run 2** (n = 50 vs n = 7). The
  run-2 holdout reported AUC 0.794 with a very wide CI (0.590–0.946). This external estimate,
  0.715 (CI 0.630–0.795), sits **inside** that run-2 interval — consistent with it, and
  pointing to a true value nearer 0.70–0.75 than the 0.794 point estimate. A modest optimism
  in the small-n run-2 figure is the expected, honest reading; the *direction* (geometry
  enriches known azoles from matched decoys, far better than the docking score) replicates.
- **The two conversion-failed actives are ranked last** (positions 380, 381), so the reported
  AUC is a **conservative floor**: had they docked, the figure could only rise. The failure is
  a tool limitation (obabel PDBQT emission on those two structures), not a scoring outcome.
- **The headline contrast holds.** On this fresh set the docking program's own affinity score
  ranks the actives at AUC 0.590 with EF@5% 0.00x — better than run 2's worse-than-random
  0.471 but still far below the geometry criterion on the identical poses. For azoles against
  CYP51, geometry beats the scoring function, now on a third independent set.

## Honest limitations (fixed in the pre-registration, restated)

1. **Same curated population.** The actives are a never-*scored* subset of the same ChEMBL
   *C. albicans* azole population as `verified_actives` — a clean never-touched-by-any-gate
   holdout, **not** a different chemical world or a later time slice. A true
   temporal/prospective holdout needs a newer ChEMBL release; recorded as the stronger
   follow-up on issue #82.
2. **Decoys are presumed inactive, not measured.** Any real binder among the 335 depresses
   the AUC; any decoy-selection bias inflates it.
3. **One rigid *C. albicans* receptor**, rigid docking (no induced fit), applicability domain
   ≤ 45 heavy atoms — the same blind spots as run 2.
4. **This validates a ranking criterion, not a drug.**

## What this earns

The overfitting objection issue #82 raises is answered: the frozen geometric criterion
recovers known antifungals from property-matched decoys on a disjoint, never-scored set,
decisively above chance (p < 0.0001) and above the docking score, with a point estimate over
the pre-registered bar. The honest caveat is the width of the interval at n = 50 and the
same-population provenance — both of which the next step (a newer-ChEMBL-release temporal
holdout) would tighten.

## Reproduce

    python scripts/make_external_holdout.py actives          # offline; seed-82 disjoint sample
    python scripts/make_external_holdout.py decoys           # fetches ChEMBL
    python scripts/prep_ligands.py data/ligands/external_actives.smi work/external_ligands
    python scripts/prep_ligands.py data/ligands/external_decoys.smi work/external_ligands
    RECEPTOR=work/receptor.pdbqt scripts/screen.sh work/external_ligands work/screen_external 1
    python scripts/validate_gate_external.py work/screen_external work/receptor_A.pdb
