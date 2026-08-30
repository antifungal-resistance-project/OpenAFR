# Pre-registration — candidate resistance-pocket retention, add-atom mutants (issue #77)

**Written 2026-08-30, BEFORE reading any nitrogen-to-iron distance from the *C. auris*
K143R or F126L candidate dockings.** Everything below is fixed in advance; deviations are
recorded as deviations, not edited in. This is the direct sibling of
[PREREGISTRATION_candidate_resistance.md](PREREGISTRATION_candidate_resistance.md) (the
Y132F run) — read it first. It inherits that run's criterion, bar and interpretation
**unchanged**; the only differences are the two receptors and how they were built.

## Background — why this run exists

The Y132F retention axis ([RESULTS_candidate_resistance.md](RESULTS_candidate_resistance.md))
closed the dominant global *C. auris* substitution, but the scorecard still prints
"auris K143R/F126L mutant retention (#77)" as `not-yet-checked` — the two clade-specific
substitutions the mission also has to survive: **K143R** (clade I) and **F126L** (clade III).
Unlike Y132F (a pure para-hydroxyl *deletion*, built with zero modeling by
`scripts/mutate_receptor.py`), both of these **add** side-chain atoms, so they cannot be built
by atom deletion — they need a side-chain repacker. That repacker did not exist in the pinned
environment, which is exactly why these two mutants were deferred.

## What changed since the Y132F run — the pinned repacker (declared honestly)

This run introduces `scripts/repack_mutant.py` and pins `pymol-open-source=3.1.0` in a
**separate** build-time env, `environment-repack.yml` (it does not co-solve with the main env's
rdkit/vina pins, and the repacker is build-time only — its output is a committed, sha-frozen
receptor PDB, so nothing downstream imports pymol). The two mutant receptors are built by
PyMOL's mutagenesis wizard, which is a **modeling step Y132F did not need**. To keep it inside
the reproducibility contract:

- **Deterministic rotamer choice, frozen here.** For the fixed backbone frame the wizard
  enumerates the backbone-dependent rotamer library and the script selects the
  **minimum-strain** rotamer (fewest steric clashes with the fixed pocket) — one reproducible
  choice, not a judgement call. Re-running the build yields a byte-identical receptor PDB
  (`tests/test_repack_mutant.py`).
- **The mutation stays the only variable.** The build is REFUSED unless every atom outside the
  mutated residue is bit-identical to `work/receptor_A.pdb` and the heme iron has not moved
  (verified in `repack_mutant.py`; observed worst non-mutated-atom delta 0.000 Å, Fe unchanged).
- **Labeled MODELED-ROTAMER wherever it surfaces.** These mutants are never conflated with the
  exact-geometry Y132F/V125A deletions. The added modeling assumption is a *limitation*
  (below), carried on the card and in the results.

## Fixed inputs

- **Ligand set:** identical to the Y132F run — the same prepared 3D ligands regenerated
  deterministically from `work/candidates/shortlist_focus.smi` (ETKDGv3 seed 42), 10 novel
  candidates + 6 known-azole controls.
- **Receptors (built BEFORE this run's grading, shas frozen here):**
  - `work/receptor_auris_K143R.pdb` — `scripts/repack_mutant.py K143R`; heavy-atom Lys→Arg,
    minimum-strain rotamer.
  - `work/receptor_auris_F126L.pdb` — `scripts/repack_mutant.py F126L`; heavy-atom Phe→Leu,
    minimum-strain rotamer.
  - Each converted with `scripts/prep_receptor.py --mutant`. Both share 5TZ1's coordinate
    frame, so the docking box does **not** move.
- **Docking box + search:** the frozen `protocol.yaml` **unchanged** — center 70.61/66.28/4.18,
  size 26, exhaustiveness 32, num_modes 20. The **only** variable versus the wild-type candidate
  dock is the mutated side chain.
- **Seeds:** {42, 1, 2, 3, 4} — a 5-seed consensus min, apples-to-apples with the WT fungal
  consensus already in `work/candidates/shortlist_confidence.json`.

## Wild-type baseline (no re-dock)

Identical to the Y132F run: the WT candidate baseline is the **already-frozen** `fungal_consensus`
in `work/candidates/shortlist_confidence.json` (same 16 ligands, same 5 seeds, same frozen
protocol, receptor `work/receptor_A.pdb`). It is **not** re-docked. Retention is the shift of
each mutant consensus relative to that fixed WT consensus.

## The measure (fixed — identical to the Y132F run)

Per molecule, per mutant, `mutant_consensus` = MIN nitrogen-to-iron distance over the seeds that
produced a coordinating pose (< FE_CUTOFF = 3.0 Å), with `n_coord/n_seeds` and `spread` — the
**exact** `consensus()` logic in `scripts/shortlist_confidence.py`, unchanged. Then:

    retention_shift = mutant_consensus − fungal_consensus   (Å; + = moved AWAY from iron)

## Pre-committed per-candidate classification (frozen, identical bar to the Y132F run)

`RETENTION_TOL = 0.5 Å`, `MIN_COORD_SEEDS = 3` — reused verbatim via
`scripts/candidate_resistance.classify` (the frozen function the Y132F tests already pin).

- **RETAINED** — coordinates the mutant iron in **≥ 3/5 seeds** AND `|retention_shift| ≤ 0.5 Å`.
- **SHIFTED** — coordinates in **≥ 3/5 seeds** but `|retention_shift| > 0.5 Å` (direction
  reported; SHIFTED-closer is never silently praised).
- **LOST** — coordinates in **< 3/5 seeds** on the mutant pocket. The informative demotion.

## Pre-committed interpretation (do NOT overclaim — same as Y132F)

- **Retention ≠ resistance-breaking.** The held-out Y132F transfer test showed the geometric
  criterion is nearly insensitive to these pocket edits, so broad retention is close to a
  tautology, not evidence a molecule *defeats* resistance. A **RETAINED** verdict means only:
  the coordination geometry that shortlisted the molecule is **not destroyed** by the
  substitution. The **LOST** verdict is the more informative outcome.
- This axis **demotes/annotates, never promotes.** Per-mutant verdicts are added to the
  scorecard; they never collapse into a composite score and never make a molecule a hit.

## Limitations (fixed in advance)

1. **Modeled rotamer, not observed geometry.** The K143R/F126L side chains are placed by the
   repacker's minimum-strain rotamer, not solved crystallographically. This is a strictly
   weaker structural claim than the Y132F deletion. A different (still low-strain) rotamer could
   change a borderline verdict; the reported verdict is conditioned on the frozen rotamer choice.
2. **Rigid backbone, single conformer, everything else frozen.** This is not an induced-fit /
   flexible-receptor model — only the one mutated side chain differs from wild type. Real
   clinical resistance also involves backbone/H-bond rearrangement and co-occurring
   efflux/expression changes not captured here.
3. **Point substitutions only.** TR-type tandem-duplication / promoter changes and non-ERG11
   loci are copy-number/expression events, not CDS pocket point mutations, and remain out of
   structural scope for this axis (noted, not silently dropped).
4. **Controls are controls.** The 6 known azoles are clinically resistance-defeated; ranking
   them well on a mutant pocket is a robustness demonstration, never a resistance claim.
