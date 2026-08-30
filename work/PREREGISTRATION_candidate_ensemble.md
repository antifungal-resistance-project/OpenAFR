# Pre-registration — candidate ensemble-consistency (issue #70)

**Written 2026-08-29, BEFORE reading any nitrogen-to-iron distance from the 5FSA / 5V5Z
candidate dockings.** Everything below is fixed in advance; deviations are recorded as
deviations, not edited in. Inherits the run-2 protocol, the consensus criterion, and the
within-azole ensemble look — read [PREREGISTRATION_run2.md](PREREGISTRATION_run2.md),
[PREREGISTRATION_consensus.md](PREREGISTRATION_consensus.md) and
[PREREGISTRATION_ensemble.md](PREREGISTRATION_ensemble.md) first.

## Background — why this run exists

The consolidated scorecard ([RESULTS_scorecard.md](RESULTS_scorecard.md)) ranks and vets its
10 novel candidates against a **single frozen *C. albicans* conformer** (5TZ1). The scorecard
itself prints "candidate-path ensemble/flexible docking (#70)" as `not-yet-checked` — the
ensemble machinery (`openafr/ensemble.py`, `scripts/prep_receptor_ensemble.py`) was built and
validated on the **held-out azole set** (looks #8/#9, [RESULTS_ensemble.md](RESULTS_ensemble.md))
but was never turned on the **candidate path** — the molecules we would actually forward.

A candidate that coordinates the heme iron across **every experimentally-observed CYP51 channel
conformation** is a stronger hypothesis than one that only fits the single snapshot it was ranked
against. This run applies the *same exhaustive crystallographic ensemble* the ceiling work used
to each shortlist candidate and reports, per molecule, an ensemble-min distance and a
conformer-consistency verdict.

This is a **per-candidate confidence axis, not a gate.** It does not justify or block a new
screen, and it does not re-open within-azole ranking (look #9 closed that on this exact
ensemble). It adds one honest column per candidate to the scorecard: does the iron-coordinating
geometry that placed a molecule on the shortlist survive being scored against the *other*
experimental conformers, or is it a single-snapshot artifact?

## Fixed inputs

- **Ligand set:** identical to the seed-stability, human-selectivity and resistance runs —
  `work/candidates/shortlist_focus.smi` (sha256 `39d8fe0f9a95…`; 10 novel candidates + 6
  known-azole controls). Same prepared 3D ligands (`work/candidates/cand_ligands/`,
  `scripts/prep_ligands.py`).
- **Ensemble (EXHAUSTIVE, not selected):** the three and only three *C. albicans* CYP51
  crystal structures in the PDB, exactly the set the ceiling work used —
  - `work/receptor_A.pdb` (5TZ1) sha256 `a1a8410868d1…` — the project anchor / WT baseline conformer;
  - `work/receptor_5FSA_A.pdb` (5FSA, posaconazole co-crystal, 2.86 Å) sha256 `9531a1531234…`;
  - `work/receptor_5V5Z_A.pdb` (5V5Z, itraconazole co-crystal, 2.9 Å) sha256 `93a261010162…`.

  The two new conformers are turned into docking receptors by
  `scripts/prep_receptor_ensemble.py 5FSA 5V5Z` (deterministic Kabsch superposition onto the
  5TZ1 frame + `obabel -xr -p 7.4`), producing `work/receptor_5FSA.pdbqt`
  (sha256 `b00970fb96cb…`) and `work/receptor_5V5Z.pdbqt` (sha256 `160b2cbe5b4a…`). Because the
  set is every available structure, there is **no conformer to select**.
- **Docking box + search:** the frozen `protocol.yaml` **unchanged** — center 70.61/66.28/4.18,
  size 26, exhaustiveness 32, num_modes 20. Every conformer shares the 5TZ1 frame (the box does
  **not** move), and each conformer's N–Fe distance is measured against **its own** heme iron.
- **Seeds:** {42, 1, 2, 3, 4} — a 5-seed consensus min per conformer, apples-to-apples with the
  WT fungal consensus already in `work/candidates/shortlist_confidence.json`.

## 5TZ1 baseline (no re-dock)

The 5TZ1 per-conformer value is the **already-frozen** `per_seed_fungal` consensus in
`work/candidates/shortlist_confidence.json` — the same 16 ligands, the same 5 seeds, the same
frozen protocol, receptor `work/receptor_A.pdb`. It is **not** re-docked here; re-docking it
would change the anchor conformer and break the controlled comparison. Only 5FSA and 5V5Z are
newly docked.

## Self-docking guard (fixed; a no-op for this ligand set, stated for honesty)

The ceiling work excluded a co-crystal azole from the conformer it was solved in (posaconazole
from 5FSA, itraconazole from 5V5Z). **None of the 16 shortlist molecules is posaconazole or
itraconazole**, so the guard excludes nothing here. It is still applied by name (via
`openafr.ensemble.ensemble_rows`'s `exclude` map) so that if the ligand set ever changes, no
molecule is silently scored against the structure it was co-crystallized in.

## The measure (fixed)

Per molecule, per conformer, `conformer_consensus` = MIN nitrogen-to-iron distance over the
seeds that produced a coordinating pose (< FE_CUTOFF = 3.0 Å), with `n_coord/n_seeds` — the
**exact** `consensus()` logic in `scripts/shortlist_confidence.py`, unchanged. A conformer
"coordinates" for a molecule iff **≥ MIN_COORD_SEEDS = 3 of 5 seeds** reach the iron there
(the same bar the resistance axis uses).

Then, over the three conformers:

    ensemble_min  = MIN conformer_consensus over conformers where the molecule had any
                    coordinating pose  (openafr.ensemble.aggregate, mode="min")
    ensemble_mean = mean of the same    (diagnostic, reported, never a verdict — prereg
                    stopping rule from the ceiling work)
    n_conf_coord  = how many of the 3 conformers coordinate (≥3/5 seeds)

## Pre-committed per-candidate classification

Frozen here, before any 5FSA/5V5Z distance is read. Based only on `n_conf_coord`:

- **ROBUST** — coordinates in **3/3** conformers. The shortlist geometry holds across every
  experimentally-observed channel conformation: the strongest form of the hypothesis.
- **CONSISTENT** — coordinates in **2/3** conformers. Holds in a majority of conformers.
- **CONFORMER-FRAGILE** — coordinates in exactly **1/3** conformers. The coordination that
  shortlisted the molecule appears in only one snapshot — a single-conformer artifact. This is
  the informative demotion this run exists to surface.
- **NON-COORDINATING** — coordinates in **0/3** conformers. Never reliably reaches the iron on
  the ensemble (expected only for a molecule already seed-unstable on 5TZ1).

`ensemble_min` is reported for every molecule regardless of verdict.

## Pre-committed interpretation (do NOT overclaim)

- **This axis measures robustness, not potency.** Look #9 already showed that the ensemble does
  **not** lift within-azole *ranking* to the 0.70 bar ([RESULTS_ensemble_confirm.md](RESULTS_ensemble_confirm.md));
  a ROBUST verdict therefore means only that a molecule's coordination geometry is reproduced
  across conformers, **not** that it binds more tightly or defeats any azole. It is a
  single-snapshot-artifact filter, nothing more.
- **The CONFORMER-FRAGILE / NON-COORDINATING verdict is the informative outcome.** It flags a
  shortlist molecule whose coordination rests on one crystal snapshot — exactly the artifact the
  ensemble exists to catch.
- **Coordination remains necessary-not-sufficient.** Reaching the iron in three rigid crystal
  conformers is still not proof of binding; every caveat of the underlying screen binds.
- This axis **demotes/annotates, never promotes.** It is added to the scorecard as one column
  with a `not-yet-checked` → measured state; it does not collapse into any composite score and
  does not make any molecule a hit.

## Pre-committed scorecard demotion

A **CONFORMER-FRAGILE** or **NON-COORDINATING** novel candidate earns one scorecard concern
("coordination is conformer-fragile / does not generalize across CYP51 conformers"). Controls
are exempt from the *concern* (their conformer behavior validates the axis rather than striking a
discovery), but their verdict is still reported.

## Limitations (fixed in advance)

1. **Crystallographic conformers only.** The ensemble is the three rigid experimental *C.
   albicans* CYP51 structures. MD-relaxed / induced-fit conformers (mentioned in #70) need an MD
   stack not in the pinned environment and remain out of scope — this is a rigid-ensemble
   robustness check, not full receptor flexibility.
2. **Three conformers is a small ensemble.** n_conf_coord is one of {0,1,2,3}; the verdict is
   coarse by construction. It distinguishes "single-snapshot" from "generalizes," not fine
   gradations.
3. ***C. albicans*, not *C. auris*.** These are the *C. albicans* structures (no experimental
   *C. auris* CYP51 exists); resistance-pocket robustness is the separate Y132F axis
   ([RESULTS_candidate_resistance.md](RESULTS_candidate_resistance.md)).
4. **Coordination necessary-not-sufficient; single rigid pose per conformer.** Every caveat of
   the underlying screen still binds.
