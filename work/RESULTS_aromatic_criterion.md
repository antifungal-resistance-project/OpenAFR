# Results — the aromatic-nitrogen criterion fix, re-validated on the run-2 gate

Date: 2026-09-02. Pre-registered in
[PREREGISTRATION_aromatic_criterion.md](PREREGISTRATION_aromatic_criterion.md) (sha256
`9a164558…`, `work/PREREG_aromatic_criterion.sha256`, frozen before any gate distance was read
under either criterion). Criterion in [`openafr/coordinator.py`](../openafr/coordinator.py)
(`min_aromatic_nitrogen_iron_distance`, tests in
[`tests/test_coordinator.py`](../tests/test_coordinator.py)); grader
[`scripts/validate_gate_aromatic.py`](../scripts/validate_gate_aromatic.py).

> **Outcome: the fix is NOT validated. Restricting the ranking criterion to the aromatic ring
> nitrogen LOWERS the validated gate's AUC from 0.794 to 0.729** (it still clears the 0.70 bar,
> but the pre-registered requirement was that it not fall). Per the pre-registered decision
> rule, the any-N criterion **stands unchanged**; the coordinator identity is kept as a
> per-candidate flag (the dossier annotation), never a global re-ranking.

## The test

The coordinator-identity axis ([RESULTS_coordinator.md](RESULTS_coordinator.md)) showed the
frozen criterion `min_nitrogen_iron_distance` (nearest nitrogen of any type) can be set by an
out-of-mechanism nitrile/amine/azide N. The proposed principled fix was to rank on the aromatic
ring N only (`min_aromatic_nitrogen_iron_distance`). To check it does not weaken the validated
result, the run-2 held-out gate (7 azole actives vs 350 property-matched decoys,
[PREREGISTRATION_run2.md](PREREGISTRATION_run2.md)) was **re-docked** (the original poses are
gitignored/gone; 355/357 converted, the same 2 obabel failures seen on the ensemble runs, 0
actives lost) and graded under **both** criteria on identical poses — the only variable is
which nitrogen the minimum may use. Protocol hash verified unmodified.

## The result

| criterion | AUC | EF@1% | EF@5% | BEDROC | actives in top 15 |
|---|---:|---:|---:|---:|---|
| **any-N** (frozen gate) | **0.794** | 12.68x | 5.63x | 0.325 | fenticonazole, **luliconazole** |
| **aromatic-N** (the fix) | **0.729** | 12.68x | 5.63x | 0.289 | fenticonazole |

The any-N number **reproduces the original run-2 held-out result exactly** (AUC 0.794, EF@1%
12.68x — [RESULTS_run2_holdout.md](RESULTS_run2_holdout.md)), so the re-dock is faithful. The
aromatic-N restriction **drops the AUC by 0.065** and loses one active from the top 15.

## Why it fell — luliconazole, a real active, coordinates via its nitrile

The pre-registration named the exact risk: **luliconazole** is a held-out active carrying BOTH
an imidazole (aromatic) and a nitrile. Its docked pose coordinates the iron through the
**nitrile at 2.764 Å**, while its imidazole ring N sits **5.281 Å** away. So the aromatic
restriction pushes luliconazole from 2.764 → 5.281 Å, demoting a genuine active and lowering
the AUC. The other six actives are pure imidazole/triazole coordinators (aromatic-N distance
identical to any-N).

The decoys did not compensate: every run-2 decoy is property-matched to an azole active and so
carries ≥1 aromatic N (#87), and for most that aromatic N was already at or near the
iron-reaching nitrogen — restricting to it did not push the decoys away enough to offset the
luliconazole loss.

**The load-bearing lesson:** a non-aromatic (nitrile) coordinator is **not** proof of an
artifact. A bona fide CYP51 active can coordinate through one in a rigid docked pose. So
"out-of-mechanism" means *the geometric evidence is unvalidated for this coordinator type*, not
*the molecule is a non-binder*.

## Consequences (both applied)

1. **The frozen ranking criterion is unchanged.** The gate, the domain, and every downstream
   rank still use `min_nitrogen_iron_distance` (any-N). No edit to `protocol.yaml`. Adopting
   aromatic-N would have thrown away real signal (AUC 0.794 → 0.729).
2. **The coordinator identity stays a per-candidate FLAG, not a re-rank.** In the dossier
   ([RESULTS_dossier.md](RESULTS_dossier.md)) a dual-mode or out-of-mechanism coordinator now
   **blocks automatic-A** and is surfaced prominently, but **does not force priority C** —
   precisely because luliconazole shows a non-aromatic coordinator can belong to a real active.
   This caps such candidates at B (blocked from the top tier, flagged), never hard-rejecting
   them.

## What this establishes / does not

- **Establishes:** that the aromatic-restricted criterion, though mechanistically cleaner, is
  empirically worse on the validated gate (AUC −0.065), driven by a real active (luliconazole)
  that coordinates via a non-aromatic N; and that the honest use of the coordinator axis is a
  blocking flag, not a criterion swap or a hard demotion.
- **Does not:** re-open the gate's pre-registered pass (which stands on any-N at AUC 0.794),
  claim any candidate is active, or edit the frozen protocol. This is a criterion-quality
  negative, disciplined by a pre-registration that named its own failure mode in advance.

## Reproduce

    cat data/ligands/actives_holdout_final.smi data/ligands/decoys_holdout.smi > work/run2_all.smi
    conda run -n openafr python scripts/prep_ligands.py work/run2_all.smi work/run2_ligands
    conda run -n openafr bash scripts/screen.sh work/run2_ligands work/screen2
    conda run -n openafr python scripts/validate_gate_aromatic.py work/screen2

Poses are gitignored; the criterion, grader, pre-registration and this result persist.
