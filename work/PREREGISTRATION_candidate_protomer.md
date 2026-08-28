# Pre-registration — per-candidate protonation / tautomer axis (issue #84)

**Written 2026-08-27, BEFORE reading any nitrogen-to-iron distance from the per-microstate
dockings.** Everything below is fixed in advance; deviations are recorded as deviations, not
edited in. Inherits the run-2 protocol and the consensus criterion — read
[PREREGISTRATION_run2.md](PREREGISTRATION_run2.md) and
[PREREGISTRATION_consensus.md](PREREGISTRATION_consensus.md) first, and it is the twin of the
stereochemistry axis [PREREGISTRATION_candidate_stereo.md](PREREGISTRATION_candidate_stereo.md)
(#85) — same machinery, different degree of freedom.

## Background — why this run exists

The consolidated scorecard ([RESULTS_scorecard.md](RESULTS_scorecard.md)) docks every candidate
as **one supplied SMILES string**, protonated once by `obabel -p 7.4` at dock time. That fixes
two chemistry degrees of freedom silently: the **protonation microstate** (which ionizable
groups carry a charge at physiological pH) and the **tautomer** (where a mobile hydrogen sits).
A molecule can adopt more than one of these in solution near pH 7.4, and a different microstate
can present a different hydrogen-bond donor/acceptor pattern to the pocket — including at the
very nitrogen the coordination criterion measures — and so dock differently. The stereo axis
(#85) explicitly deferred this to #84 (its limitation 3). Before trusting a rank we must know
whether it is **microstate-stable** — whether the coordinating geometry survives the protonation
and tautomer states the SMILES + `-p 7.4` left unexamined.

This is a **per-candidate confidence axis, not a gate.** It does not justify or block a screen.
It adds one honest column per candidate to the scorecard: is this rank tied to a single relevant
microstate at pH 7.4, and (if not) is the coordination geometry microstate-robust or
microstate-dependent?

## Fixed inputs

- **Molecule set:** identical to every other scorecard axis —
  `work/candidates/shortlist_focus.smi` (10 novel candidates + 6 known-azole controls).
- **Reference microstate (the docked scorecard form):** `obabel -p 7.4` of the supplied SMILES —
  the exact protonation the scorecard dock used. It is always a member of the enumerated set;
  the enumerator asserts this.
- **Protonation enumeration (deterministic, pinned-env):** the distinct microspecies returned by
  `obabel -p {6.4, 7.4, 8.4}` — Open Babel's own built-in pKa model (the same model the dock
  relies on), sampled across the physiological band pH 7.4 ± 1.0. A group that does not change
  protonation anywhere in 6.4–8.4 is treated as fixed; a group that titrates within that band is
  a genuine second microstate. **No external pKa predictor is added** (that would break the
  pinned environment); the axis therefore closes the "protonation obabel itself sees titrating
  near pH 7.4" hole, not the "obabel's pKa model is wrong" hole (limitation 1).
- **Tautomer enumeration (deterministic, pinned-env):** RDKit
  `MolStandardize.rdMolStandardize.TautomerEnumerator` (default transform set), `MaxTautomers=16`,
  applied to **each** protonation microstate; every distinct enumerated tautomer is kept and
  docked. No energy/score threshold is applied to prune tautomers: docking *every* enumerated
  tautomer — including ones RDKit scores as minor — is the conservative choice, so a ROBUST
  verdict means robust even against implausible tautomers, and a MICROSTATE-DEPENDENT verdict
  records *which* tautomer diverged (its SMILES) so a human can judge whether that tautomer is
  physiologically populated (limitation 4).
- **Microstate set:** the distinct canonical SMILES over the (protonation × tautomer) product,
  capped at `MAX_STATES = 16`. `protomer_explicit` ≡ the set has **exactly one** member (the
  reference): protonation is stable across 6.4–8.4 **and** there is a single tautomer. Produced by
  `scripts/enumerate_protomers.py` (network-free, docking-free — a pure function of the input
  SMILES + obabel/RDKit, so it is not a gradeable result and may be computed before this freeze).
- **Microstate embedding:** each microstate SMILES embedded to one 3D conformer with RDKit
  ETKDGv3, fixed seed `0xF00D`, MMFF-optimized, written as `<candidate>__stNN.sdf`. stNN indices
  follow canonical-SMILES sort order, so they are stable across runs. Docked exactly as every
  candidate is (`-p 7.4` is re-applied at dock time, unchanged — it is a no-op on an
  already-correct protonation and keeps this axis byte-identical to the main pipeline).
- **Receptor:** the wild-type *C. albicans* docking receptor `work/receptor.pdbqt`
  (`scripts/prep_receptor.py`, anchor `work/receptor_A.pdb`) — the **same** receptor the shortlist
  rank was measured on. This axis asks whether the rank changes with the microstate, so the
  receptor must be the ranking receptor, not a mutant.
- **Docking box + search:** the frozen `protocol.yaml` **unchanged** — center 70.61/66.28/4.18,
  size 26, exhaustiveness 32, num_modes 20. The **only** variable versus the shortlist dock is the
  ligand protonation/tautomer state.
- **Seeds:** {42, 1, 2, 3, 4} — a 5-seed consensus min, apples-to-apples with the fungal consensus
  already in `work/candidates/shortlist_confidence.json`.

## The measure (fixed)

Per microstate, `state_consensus` = MIN nitrogen-to-iron distance over the seeds that produced a
coordinating pose (< FE_CUTOFF = 3.0 Å), with `n_coord/n_seeds` and `spread` — the **exact**
`consensus()` logic in `scripts/shortlist_confidence.py`, unchanged. Per candidate:

    state_spread = max(state_consensus) − min(state_consensus)   over microstates that coordinated

## Pre-committed per-candidate classification

Frozen here, before any per-microstate distance is read. Mirrors the stereo axis exactly.

- **EXPLICIT** — `protomer_explicit` is true: at pH 7.4 the molecule has a single relevant
  microstate (protonation stable across 6.4–8.4 AND one tautomer). The shortlist rank is
  microstate-explicit; **no extra docking is performed or needed.** Expected for rigid,
  non-ionizable, tautomer-locked molecules; it is *not* a strike — it means the axis has nothing
  to demote.
- For a molecule that is **not** `protomer_explicit`, all microstates are docked and:
  - **ROBUST** — **every** enumerated microstate coordinates the iron (each ≥ 3/5 seeds) AND
    `state_spread ≤ STATE_TOL`. The coordination geometry the rank rests on does not depend on
    which protonation/tautomer was docked.
  - **MICROSTATE-DEPENDENT** — the microstates disagree: at least one coordinates and at least one
    fails (< 3/5 seeds), OR `state_spread > STATE_TOL`. The rank is sensitive to a protonation or
    tautomer the pipeline fixed silently — the informative flag this axis exists to surface.

`STATE_TOL = 1.0 Å`, chosen equal to the `SPREAD_MAX = 1.0 Å` stable-consensus tolerance (#68)
and to the stereo axis's `ISO_TOL` — a microstate spread within the run-to-run noise band already
tolerated for a "stable" candidate is treated as materially unchanged; beyond it, the microstate
choice matters. `MIN_COORD_SEEDS = 3` (of 5), identical to the stereo axis.

## Pre-committed interpretation (do NOT overclaim)

- **EXPLICIT is the honest headline, not a pass.** A microstate-explicit rank is unambiguous *as
  to which protonation/tautomer was scored*; it says nothing about whether that microstate is the
  bioactive one. It only removes the specific failure #84 names: a rank silently tied to one
  protonation/tautomer when others are populated at pH 7.4.
- **A control that turns out MICROSTATE-DEPENDENT exercises the machinery, not a candidate claim.**
  As with stereo, the controls are the positive-control workout for the docking half.
- **This axis demotes/annotates, never promotes.** It adds one column with a `not-yet-checked`
  → measured state; it does not collapse into any composite score and does not make any molecule a
  hit.

## Limitations (fixed in advance)

1. **Protonation enumeration is only as good as Open Babel's pKa model.** The axis tests the
   microstates obabel itself flips across pH 6.4–8.4; it does not second-guess that model, nor
   enumerate microstates from a group whose true pKa obabel places outside the band. It closes the
   "protonation obabel sees titrating near physiological pH" hole, not the "obabel's pKa is wrong"
   hole. No external pKa tool is added, to keep the pinned environment intact.
2. **One embedded conformer per microstate**, then the frozen search — the usual single-conformer,
   rigid-receptor, coordination-necessary-not-sufficient caveats all still bind.
3. **Stereochemistry (#85) is a separate, already-measured axis** and is out of scope here; each
   microstate inherits the supplied SMILES's stereochemistry unchanged.
4. **Every enumerated tautomer is docked without an energy filter.** This is deliberately
   conservative (ROBUST-against-everything is the strong statement), but it means a
   MICROSTATE-DEPENDENT verdict may be driven by a high-energy tautomer that is not meaningfully
   populated; the diverging microstate's SMILES is recorded so that judgement can be made
   downstream rather than hidden. `MAX_STATES = 16` caps a combinatorially rich input.
