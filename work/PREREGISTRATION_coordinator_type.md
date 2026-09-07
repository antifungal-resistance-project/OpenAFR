# Pre-registration — the type-aware coordinator criterion, on the run-2 gate

**Written 2026-09-05, BEFORE docking the run-2 gate set under the new criterion or reading
any type-aware distance.** Fixes the hypothesis, the grading, and the decision rule in
advance. Inherits the run-2 protocol and pass bar unchanged — read
[PREREGISTRATION_run2.md](PREREGISTRATION_run2.md),
[PREREGISTRATION_coordinator.md](PREREGISTRATION_coordinator.md), and
[PREREGISTRATION_aromatic_criterion.md](PREREGISTRATION_aromatic_criterion.md) first.

## Background — the gap this probes

Two prior results bracket the frozen any-N criterion
(`openafr.pdbqt.min_nitrogen_iron_distance`, nearest nitrogen of any type):

- **Any-N over-credits out-of-mechanism coordinators.** On the novel shortlist
  ([RESULTS_coordinator.md](RESULTS_coordinator.md)) the rank of verinurad, taranabant,
  lersivirine, olprinone was set by a **nitrile / amine / azide** N, not the aromatic ring N
  the azole mechanism was validated on.
- **Aromatic-only under-credits a real active.** Restricting to the aromatic ring N
  ([RESULTS_aromatic_criterion.md](RESULTS_aromatic_criterion.md)) *lowered* the run-2 gate
  AUC 0.794 → 0.729, because the active **luliconazole coordinates the iron through its
  nitrile in-pose** (2.764 Å; its imidazole is 5.28 Å away).

So neither endpoint is right: any-N admits amine/azide artifacts, aromatic-only throws out a
nitrile that a real drug genuinely uses. This run tests the principled middle.

## The claim being tested

A **type-aware** criterion — `openafr.coordinator.min_coordinating_nitrogen_iron_distance`
with `kinds = ("aromatic-ring", "nitrile")` — admits only nitrogen kinds with direct evidence
of type-II heme coordination on **this** system (an aromatic ring N: the azole mechanism; a
nitrile N: the mechanism luliconazole uses) and excludes sp3 **amine** and terminal **azide**
(no such evidence; the artifacts any-N rewards). The hypothesis:

> **Type-aware AUC ≥ any-N AUC (0.794)** on the run-2 held-out gate, and the gate still passes
> (AUC ≥ 0.70, EF@1% ≥ 5.0x).

Mechanism of the two contrasts it should win:
- **vs aromatic-only:** re-admitting the nitrile rescues luliconazole → AUC should climb back
  above 0.729.
- **vs any-N:** excluding amine/azide pushes away any decoy that scored through a
  non-coordinating N → AUC should not fall below 0.794, and may rise if such decoys exist.

**Named risks (measured, not asserted):**
1. **The gate decoys may not coordinate via amine/azide.** If the property-matched decoys mostly
   present aromatic or nitrile N anyway, excluding amine/azide changes nothing and type-aware
   simply reproduces any-N (0.794). That is a *null*, not a failure — it still validates the
   criterion as safe to adopt for novel-candidate ranking (where the amine/azide artifacts
   demonstrably do occur).
2. **Nitrile re-admission could pull a decoy up.** A nitrile-bearing decoy (e.g.
   `decoy_CHEMBL6328_for_efinaconazole`) that coordinates through its nitrile would rank better
   under type-aware than under aromatic-only. If enough do, type-aware could sit *below* any-N.
   The gate decides on the measured AUC.
3. **Azide topology.** A terminal azide N (N=N ~1.24 Å) sits near the nitrile bond-length
   threshold (`BOND_NITRILE` 1.25 Å) and may be misclassified `nitrile` by `classify_nitrogen`.
   R-1479 is the only azide in scope. This is a known classifier edge, reported if it bites.

## Fixed inputs (identical to run 2 — only the ranking criterion changes)

- **Ligands:** `data/ligands/actives_holdout_final.smi` (7 held-out azole actives) +
  `data/ligands/decoys_holdout.smi` (350 property-matched decoys) = 357, the run-2 gate set.
- **Receptor / protocol:** WT *C. albicans* 5TZ1 (`work/receptor_A.pdb`), the hash-frozen
  `protocol.yaml` (26 Å box, exhaustiveness 32, num_modes 20, seed 42), **read-only and
  unchanged** — `protocol.verify()` asserts the hash at grade time; this run NEVER edits it.
- **Criteria compared, same poses:** (a) any-N = `min_nitrogen_iron_distance` (the frozen
  gate); (b) aromatic-N = `min_aromatic_nitrogen_iron_distance`; (c) type-aware =
  `min_coordinating_nitrogen_iron_distance` (aromatic-ring + nitrile) — the hypothesis.
- **Pass bar:** the frozen run-2 bar from `protocol.yaml` — AUC ≥ 0.70 AND EF@1% ≥ 5.0x.
- **Unrankable rule (unchanged):** a molecule with no admissible coordinator reaching the
  pocket is ranked **LAST**, never dropped — the run-2 rule verbatim.

## Fixed decision rule

Type-aware is **adopted for novel-candidate ranking** (the dossier / shortlist scoring — NEVER
retro-applied to the frozen gate, whose pass stands on any-N, exactly as the aromatic-criterion
pre-reg ruled) iff BOTH hold on the measured re-dock:
1. **Not worse than any-N** — type-aware AUC ≥ 0.794 (within float tolerance), AND
2. **Gate still passes** — type-aware AUC ≥ 0.70 AND type-aware EF@1% ≥ 5.0x.

Beating aromatic-only alone (≥ 0.729) is reported but is **not sufficient** for adoption.
If neither adoption condition holds, it is reported as a negative and any-N stands; the
coordinator identity remains a per-candidate flag only.

## The pipeline hardening this run also delivers (independent of the outcome)

The grader writes a **per-molecule, per-nitrogen-kind distance cache**
(`work/candidates/run2_coordinator_distances.json`): for each of the 357 molecules, the minimum
Fe distance achieved by each nitrogen kind over all docked poses. That cache is a **sufficient
statistic for any kind-set criterion** — so once this dock runs, every future coordinator-criterion
idea is a *local re-grade from the cache* (`--cache`), never another dock. This retires the
"run-2 poses are gitignored/gone, so each criterion costs a fresh cloud dock" wall for this gate.

## What this run will and will not establish

- **Will:** measure whether a type-aware coordinator criterion preserves (or improves) the
  validated actives-vs-decoys separation on the pre-registered gate, on identical poses; and
  persist the per-kind cache.
- **Will not:** re-derive the domain, edit `protocol.yaml`, re-run other gates, or claim any
  candidate is active. This is a criterion-quality result, not an activity result.

## Reproduce

    # 1. one cloud dock of the run-2 gate set (exhaustiveness 32), then build the cache + grade:
    cat data/ligands/actives_holdout_final.smi data/ligands/decoys_holdout.smi > work/run2_all.smi
    conda run -n openafr python scripts/prep_ligands.py work/run2_all.smi work/run2_ligands
    conda run -n openafr bash scripts/screen.sh work/run2_ligands work/screen2
    conda run -n openafr python scripts/validate_gate_coordinator_type.py work/screen2

    # 2. thereafter, fully local — re-grade any kind-set from the cached distances, no dock:
    conda run -n openafr python scripts/validate_gate_coordinator_type.py \
        --cache work/candidates/run2_coordinator_distances.json
