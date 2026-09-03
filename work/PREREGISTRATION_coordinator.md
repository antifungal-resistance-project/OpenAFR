# Pre-registration — which nitrogen coordinates the iron (coordinator identity)

**Written 2026-09-02.** Fixes the classification rule and the per-candidate readout
**before** the full-shortlist coordinator run is read. Inherits the run-2 protocol, the
coordination criterion, and the consensus (5-seed) design — read
[PREREGISTRATION_run2.md](PREREGISTRATION_run2.md),
[PREREGISTRATION_consensus.md](PREREGISTRATION_consensus.md) and
[PREREGISTRATION_pose_convergence.md](PREREGISTRATION_pose_convergence.md) first.

## Honesty disclosure — what was already seen

This axis was motivated by, and first examined on, **verinurad** during the 2026-09-02
"harden the single Priority-A candidate" pass. That exploratory look is what surfaced the
question, and its verinurad answer is therefore **not** blind. What this pre-registration
fixes blind is (a) the *rule* — which carries no free parameter tunable to any desired
answer (see below) — and (b) the run over the **other 15 focus molecules**, including the
6 known-azole controls that validate the classifier. The verinurad numbers are reported as
what they are: a confirmatory re-run of an exploratory finding, not a blind test.

## Background — why this run exists

The validated criterion is `openafr.pdbqt.min_nitrogen_iron_distance`: the closest approach
of **any** nitrogen in a pose to the heme iron. It is element-only by design — it does not
ask *which* nitrogen. But the method was validated on the **azole mechanism**: an aromatic
sp2 ring nitrogen (imidazole / triazole / pyridine) donating its lone pair to the iron. A
molecule that instead points a **nitrile** (C#N), an **amine**, or an **amide** nitrogen at
the iron can post an identical tight N–Fe distance while coordinating through chemistry the
method never validated. Its geometry rank would then be an **out-of-mechanism artifact**.

The applicability domain already requires ≥1 aromatic N (mechanistic scope, #87). This axis
tests whether the aromatic N is the nitrogen that **actually** sets the reported distance, or
whether a non-aromatic nitrogen wins the minimum and inflates the rank.

## Fixed classification rule (no free parameter)

For a nitrogen atom in a docked pose, count its bonded heavy-atom neighbours (any heavy atom
within **1.85 Å** — longer than every C–N/C–C/C–O/C–S single bond, shorter than any
non-bonded contact) and take the shortest such bond:

- **nitrile** — exactly 1 neighbour, shortest bond **< 1.25 Å** (a C#N triple bond, ~1.16 Å).
- **amine** — exactly 1 neighbour, shortest bond ≥ 1.25 Å (a single C–N, ~1.47 Å).
- **aromatic-ring** — 2+ neighbours, shortest bond < 1.45 Å (aromatic C–N, ~1.34 Å).
- **amide-or-other** — anything else (2+ neighbours all longer, or an isolated N).

Only **aromatic-ring is in-mechanism**. The one threshold that separates the classes,
**1.25 Å**, sits in the wide empty gap between a triple bond (~1.16 Å) and every single or
aromatic C–N bond (≥ 1.34 Å): there is no value in that gap that changes any real molecule's
label, so it is a physical constant, not a knob. The rule is pure pose geometry — no
atom-order bookkeeping through obabel, no external atom labels. Implemented in
`openafr/coordinator.py`, locked by `tests/test_coordinator.py`.

## Fixed inputs

- **Ligand set:** `work/candidates/shortlist_focus.smi` — the identical 16 molecules
  (10 novel candidates + 6 known-azole controls) used by the seed-stability (#68),
  human-selectivity (#73), pose-convergence (#71) and dossier runs.
- **Receptor / protocol:** WT *C. albicans* 5TZ1 (`work/receptor.pdbqt` /
  `work/receptor_A.pdb`), the hash-frozen `protocol.yaml` (26 Å box, exhaustiveness 32,
  num_modes 20), **unchanged**.
- **Seeds:** the frozen consensus set {42, 1, 2, 3, 4}, matching #68.

## Fixed readout (per candidate, decided in advance)

For each molecule, over all 20 poses × 5 seeds:

1. **coordinator kind** — the class of the nitrogen achieving the min N–Fe in the frozen
   seed-42 dock (the pose the scorecard/dossier reports), plus the modal kind across seeds.
2. **in-mechanism distance** — the best (smallest) Fe distance reached by an
   **aromatic-ring** nitrogen, over all poses/seeds. This is the distance the molecule would
   post if the criterion were restricted to the validated mechanism.
3. **gap** = in-mechanism distance − reported (any-N) distance. A large positive gap means
   the reported rank is carried by a non-aromatic nitrogen.
4. **verdict (annotation, never a gate; never promotes):**
   - **in-mechanism** — the reported coordinator is aromatic-ring (rank is azole-like).
   - **out-of-mechanism** — the reported coordinator is nitrile/amine/amide AND the
     in-mechanism distance falls outside the validated actives' band (2.47–2.88 Å): the
     rank is an artifact and the molecule does not coordinate azole-like within the band.
   - **dual-mode** — reported coordinator is non-aromatic BUT the aromatic N still reaches
     inside the band (2.47–2.88 Å) in some pose: the rank is inflated, yet an azole-like
     binding mode is available; re-rank on the in-mechanism distance.

## Controls (validate the classifier, not discoveries)

The 6 known azoles (fluconazole, voriconazole, ravuconazole, letrozole, flutrimazole,
ketoconazole) coordinate through a **triazole or imidazole ring nitrogen** by definition.
The classifier **must** call their reported coordinator `aromatic-ring` (in-mechanism). If it
labels a known azole out-of-mechanism, the classifier — not the molecule — is wrong. Controls
are exempt from any demotion; their role is to exercise the machinery.

## What this run will and will not establish

- **Will:** label, per candidate, whether the reported iron-coordinating nitrogen is the
  validated aromatic mechanism or an out-of-mechanism nitrile/amine/amide; recover the correct
  aromatic-N call on all 6 azole controls; and, for any out-of-mechanism candidate, report the
  in-mechanism distance behind its rank.
- **Will not:** re-run the gate, re-derive AUC, promote any molecule, or relax any standing
  caveat (single rigid conformer, coordination necessary-not-sufficient, per-dock variance).
  A `dual-mode` or `in-mechanism` label is **not** an activity claim.

## Reproduce

    conda run -n openafr python scripts/prep_receptor.py
    conda run -n openafr python scripts/prep_ligands.py \
        work/candidates/shortlist_focus.smi work/candidates/coord_ligands
    conda run -n openafr bash scripts/screen_consensus.sh \
        work/candidates/coord_ligands work/candidates/dock_coordinator "42 1 2 3 4"
    conda run -n openafr python scripts/coordinator_identity.py \
        --docks work/candidates/dock_coordinator \
        --ligands work/candidates/shortlist_focus.smi --json
