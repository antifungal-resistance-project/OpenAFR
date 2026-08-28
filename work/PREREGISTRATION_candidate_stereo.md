# Pre-registration — per-candidate stereochemistry axis (issue #85)

**Written 2026-08-27, BEFORE reading any nitrogen-to-iron distance from the per-isomer
dockings.** Everything below is fixed in advance; deviations are recorded as deviations, not
edited in. Inherits the run-2 protocol and the consensus criterion — read
[PREREGISTRATION_run2.md](PREREGISTRATION_run2.md) and
[PREREGISTRATION_consensus.md](PREREGISTRATION_consensus.md) first.

## Background — why this run exists

The consolidated scorecard ([RESULTS_scorecard.md](RESULTS_scorecard.md)) docks every
candidate as **one supplied SMILES string**. Issue #85 flags the risk: if a SMILES omits or
arbitrarily fixes a stereocenter, the rank refers to *some* isomer, and a different isomer of
the same molecule can fit the pocket very differently (and have very different real-world
activity/tox). Before trusting a rank we must know whether it is **stereochemistry-explicit**
— whether the shortlist SMILES pins a single defined isomer — and, where it does not, whether
the iron-coordinating geometry survives across all enumerated isomers or depends on which one
was docked.

This is a **per-candidate confidence axis, not a gate.** It does not justify or block a screen.
It adds one honest column per candidate to the scorecard: is this rank tied to a single
defined isomer, and (if not) is the coordination geometry isomer-robust or isomer-dependent?

## Fixed inputs

- **Molecule set:** identical to every other scorecard axis —
  `work/candidates/shortlist_focus.smi` (10 novel candidates + 6 known-azole controls).
- **Enumeration (deterministic, pre-committed):** RDKit
  `FindMolChiralCenters(includeUnassigned=True, useLegacyImplementation=False)` for
  tetrahedral centers and `FindPotentialStereoBonds` for double-bond (E/Z) stereo; distinct
  isomers from the **unassigned** centers via `EnumerateStereoisomers(onlyUnassigned=True,
  unique=True, maxIsomers=16)`. `stereo_explicit` ≡ **zero** unassigned tetrahedral centers
  **and** zero unassigned stereo bonds in the supplied SMILES. Produced by
  `scripts/enumerate_stereo.py` (network-free, docking-free — a pure function of the input
  SMILES, so it is not a gradeable result and may be computed before this freeze).
- **Isomer embedding:** each enumerated isomeric SMILES embedded to one 3D conformer with
  RDKit ETKDGv3, fixed seed `0xF00D`, MMFF-optimized, written as `<candidate>__isoNN.sdf`.
  isoNN indices follow canonical-SMILES sort order, so they are stable across runs.
- **Receptor:** the wild-type *C. albicans* docking receptor `work/receptor.pdbqt`
  (`scripts/prep_receptor.py`, anchor `work/receptor_A.pdb`) — the **same** receptor the
  shortlist rank was measured on. This axis asks whether the rank changes with the isomer, so
  the receptor must be the ranking receptor, not a mutant.
- **Docking box + search:** the frozen `protocol.yaml` **unchanged** — center
  70.61/66.28/4.18, size 26, exhaustiveness 32, num_modes 20. The **only** variable versus the
  shortlist dock is the ligand stereochemistry.
- **Seeds:** {42, 1, 2, 3, 4} — a 5-seed consensus min, apples-to-apples with the fungal
  consensus already in `work/candidates/shortlist_confidence.json`.

## The measure (fixed)

Per isomer, `iso_consensus` = MIN nitrogen-to-iron distance over the seeds that produced a
coordinating pose (< FE_CUTOFF = 3.0 Å), with `n_coord/n_seeds` and `spread` — the **exact**
`consensus()` logic in `scripts/shortlist_confidence.py`, unchanged. Per candidate:

    iso_spread = max(iso_consensus) − min(iso_consensus)   over isomers that coordinated

## Pre-committed per-candidate classification

Frozen here, before any per-isomer distance is read.

- **EXPLICIT** — `stereo_explicit` is true: the supplied SMILES already pins a single defined
  isomer. The shortlist rank is stereochemistry-explicit; **no isomer docking is performed or
  needed.** This is the expected majority outcome and is *not* a strike — it means the axis has
  nothing to demote.
- For a molecule that is **not** `stereo_explicit`, all isomers are docked and:
  - **ROBUST** — **every** enumerated isomer coordinates the iron (each ≥ 3/5 seeds) AND
    `iso_spread ≤ ISO_TOL`. The coordination geometry the rank rests on does not depend on
    which isomer was docked.
  - **ISOMER-DEPENDENT** — the isomers disagree: at least one coordinates and at least one
    fails (< 3/5 seeds), OR `iso_spread > ISO_TOL`. The rank is sensitive to a stereochemistry
    the SMILES left unspecified — the informative flag this axis exists to surface.

`ISO_TOL = 1.0 Å`, chosen equal to the `SPREAD_MAX = 1.0 Å` stable-consensus tolerance
(#68): an isomer spread within the run-to-run noise band already tolerated for a "stable"
candidate is treated as materially unchanged; beyond it, the isomer choice matters.

## Pre-committed interpretation (do NOT overclaim)

- **EXPLICIT is the honest headline, not a pass.** A stereo-explicit rank is unambiguous *as
  to which isomer was scored*; it says nothing about whether that isomer is the bioactive one,
  or whether the SMILES's fixed configuration is the correct one. It only removes the specific
  failure #85 names: a rank silently averaged over, or pinned to the wrong member of, an
  unspecified stereo set.
- **The controls are expected to be the ambiguous ones and that is by design.** Racemic /
  underspecified azole drugs (e.g. ketoconazole's *cis* pair supplied without configuration)
  are the natural positive controls for the isomer-docking machinery — they exercise the
  ROBUST/ISOMER-DEPENDENT branch that the (stereo-explicit) novel candidates do not reach.
  Ranking a control's isomers is a demonstration the axis works, never a candidate claim.
- **This axis demotes/annotates, never promotes.** It adds one column with a
  `not-yet-checked` → measured state; it does not collapse into any composite score and does
  not make any molecule a hit.

## Limitations (fixed in advance)

1. **Enumeration is over UNASSIGNED centers only.** A molecule whose SMILES *asserts* a
   configuration is taken at its word — this axis does not second-guess a supplied `@`/`@@`
   or test the enantiomer of an already-specified center. It closes the "unspecified stereo"
   hole, not the "supplied the wrong configuration" hole.
2. **One embedded conformer per isomer**, then the frozen search — the usual single-conformer,
   rigid-receptor, coordination-necessary-not-sufficient caveats all still bind.
3. **Tautomers/protonation are a separate axis (#84)** and are out of scope here; `-p 7.4`
   protonation is applied at dock time exactly as for every other candidate, unchanged.
4. **maxIsomers cap = 16.** No shortlist molecule approaches it (the ambiguous ones enumerate
   to 2 and 4), but a combinatorially stereo-rich input would be capped and flagged.
