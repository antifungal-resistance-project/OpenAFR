# Pre-registration — the aromatic-nitrogen criterion fix, re-validated on the run-2 gate

**Written 2026-09-02, BEFORE re-docking the run-2 gate set or reading any gate distance
under either criterion.** Fixes the hypothesis, the grading, and the decision rule in
advance. Inherits the run-2 protocol and pass bar unchanged — read
[PREREGISTRATION_run2.md](PREREGISTRATION_run2.md) and
[PREREGISTRATION_coordinator.md](PREREGISTRATION_coordinator.md) first.

## Background — the defect this fixes

The frozen criterion `openafr.pdbqt.min_nitrogen_iron_distance` is the closest approach of
**any** nitrogen to the heme iron. The coordinator-identity axis
([RESULTS_coordinator.md](RESULTS_coordinator.md)) showed that on a novel library this can be
set by a **nitrile / amine / azide** nitrogen — an out-of-mechanism interaction the method
was never validated on (it was validated on an aromatic ring N donating its lone pair). The
principled fix is to rank on the **aromatic-ring nitrogen** only:
`openafr.coordinator.min_aromatic_nitrogen_iron_distance`.

## The claim being tested

Restricting the criterion to the aromatic ring nitrogen must not weaken the **validated**
result. The validated actives (run-2's 7 held-out azoles + the 8 training azoles) coordinate
through a triazole/imidazole ring nitrogen, so their distance should be ~unchanged; the
property-matched decoys include nitrile-bearing molecules (e.g.
`decoy_CHEMBL6328_for_efinaconazole`, a C#N decoy) that can only move **away** from the iron
under the restriction. So the aromatic-restricted AUC should be **≥** the any-N AUC.

**Named risk (the one active that could break this):** **luliconazole** is an active that
carries BOTH an imidazole (aromatic) and a nitrile. If its docked minimum-N pose coordinates
through the nitrile rather than the imidazole, the restriction would push luliconazole away
and could *lower* the AUC. This is a real possibility, not a rhetorical one — it is exactly
why this is measured on re-docked poses, not asserted. The gate decides on the measured AUC.

## Fixed inputs (identical to run 2 — only the criterion changes)

- **Ligands:** `data/ligands/actives_holdout_final.smi` (7 held-out azole actives) +
  `data/ligands/decoys_holdout.smi` (350 property-matched decoys) = 357, the run-2 gate set.
- **Receptor / protocol:** WT *C. albicans* 5TZ1 (`work/receptor.pdbqt` / `work/receptor_A.pdb`),
  the hash-frozen `protocol.yaml` (26 Å box, exhaustiveness 32, num_modes 20, seed 42),
  **unchanged**. Re-docked because the original run-2 poses are gitignored/gone; the protocol
  hash is verified at grade time.
- **Criteria compared, same poses:** (a) any-N =
  `openafr.pdbqt.min_nitrogen_iron_distance` (the frozen gate); (b) aromatic-N =
  `openafr.coordinator.min_aromatic_nitrogen_iron_distance` (the fix).
- **Pass bar:** the frozen run-2 bar from `protocol.yaml` — AUC ≥ 0.70 AND EF@1% ≥ 5.0x.
- **Unrankable rule (unchanged):** a molecule the criterion cannot score (no aromatic-ring N
  reaching the pocket, under the fix; no nitrogen at all, under any-N) is ranked **LAST**,
  never dropped — the run-2 rule verbatim.

## Fixed decision rule

The fix is **validated** iff BOTH hold on the measured re-dock:
1. **AUC preserved** — aromatic-N AUC ≥ any-N AUC (within float tolerance), AND
2. **Gate still passes** — aromatic-N AUC ≥ 0.70 AND aromatic-N EF@1% ≥ 5.0x.

If the fix is validated, the aromatic criterion is adopted for **novel-candidate ranking**
(the dossier / shortlist) — NOT retro-applied to the frozen gate, whose pre-registered pass
already stands on the any-N criterion; the fix is an improvement layer for scoring novel,
non-azole chemotypes, exactly where the artifact lives. If the fix is NOT validated (e.g.
luliconazole is demoted enough to drop the AUC), it is reported as a negative and the any-N
criterion stands, with the coordinator-identity axis kept as the per-candidate flag.

## What this run will and will not establish

- **Will:** measure whether the aromatic-restricted criterion preserves the validated
  actives-vs-decoys separation on the pre-registered gate set, on identical poses.
- **Will not:** re-derive the domain, re-run other gates (auris, ensemble, temporal…), edit
  `protocol.yaml`, or claim any candidate is active. Preservation of AUC is a
  criterion-quality result, not an activity result.

## Reproduce

    cat data/ligands/actives_holdout_final.smi data/ligands/decoys_holdout.smi > work/run2_all.smi
    conda run -n openafr python scripts/prep_ligands.py work/run2_all.smi work/run2_ligands
    conda run -n openafr bash scripts/screen.sh work/run2_ligands work/screen2
    conda run -n openafr python scripts/validate_gate_aromatic.py work/screen2
