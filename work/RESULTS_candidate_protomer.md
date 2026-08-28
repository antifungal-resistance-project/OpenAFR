# Results — per-candidate protonation / tautomer axis (#84)

Date: 2026-08-28. Pre-registered in
[PREREGISTRATION_candidate_protomer.md](PREREGISTRATION_candidate_protomer.md)
(sha256 854757d9…, frozen BEFORE any per-microstate distance was read;
`work/PREREG_candidate_protomer.sha256`). Protocol hash unchanged from run 2. Twin of the
stereochemistry axis [RESULTS_candidate_stereo.md](RESULTS_candidate_stereo.md) (#85) — same
machinery, protonation/tautomer state instead of stereoisomer.

> **This is a per-candidate confidence axis, not a gate.** It annotates the scorecard; it never
> promotes. It answers issue #84: is each shortlist rank tied to a single **relevant** microstate
> at pH 7.4, or could a protonation/tautomer the pipeline fixed silently (`obabel -p 7.4` + the
> supplied-SMILES tautomer) be making or breaking the iron-coordinating geometry the rank rests on?

## What was tested

Every one of the 16 scorecard molecules was enumerated straight from its supplied SMILES
(`work/candidates/shortlist_focus.smi`, `scripts/enumerate_protomers.py`, deterministic,
docking-free) into the set of protonation/tautomer microstates plausibly populated at pH 7.4:

- **Protonation** — the distinct microspecies Open Babel returns at pH **6.4, 7.4, 8.4**
  (`obabel -p`, its own built-in pKa model), i.e. the states that titrate across the physiological
  band pH 7.4 ± 1. The pH-7.4 microspecies is the *reference* — the exact form the scorecard dock
  used. No external pKa predictor is added, so the pinned environment is unchanged.
- **Tautomers** — RDKit `TautomerEnumerator` (default transforms, MaxTautomers=16) of each
  protonation microstate; every distinct enumerated tautomer is kept and docked (no energy pruning
  — the conservative choice: ROBUST then means robust even against implausible tautomers).

`protomer_explicit` ≡ the (protonation × tautomer) product collapses to a single microstate. For
the molecules that were **not** protomer-explicit, every microstate was embedded (ETKDGv3, fixed
seed) and docked into the **same wild-type receptor the shortlist rank was measured on**
(`work/receptor.pdbqt`), under the frozen protocol, 5 seeds {42,1,2,3,4} — the **only** variable
versus the shortlist dock is the ligand protonation/tautomer state. Per-microstate consensus N–Fe
uses the frozen `consensus()` logic unchanged.

Verdicts (frozen): **EXPLICIT** = single relevant microstate, no dock needed; **ROBUST** = every
microstate coordinates (≥3/5 seeds) and the per-microstate consensus spread ≤ 1.0 Å;
**MICROSTATE-DEPENDENT** = the microstates disagree (some coordinate, some don't) or the spread > 1.0 Å.

## The result

**Protonation is stable across the physiological band for 15 of 16 molecules** — only
flucloxacillin (a carboxylic acid) titrates between pH 6.4–8.4. **Tautomers are the live
variable.** 6 of 16 molecules are protomer-explicit (single microstate, no dock); the other 10
were docked per microstate (47 microstates total).

### Novel candidates (10)

| candidate | states | verdict | note |
|---|:--:|---|---|
| verinurad | 1 | **EXPLICIT** | single microstate |
| lersivirine | 1 | **EXPLICIT** | single microstate |
| L-838417 | 1 | **EXPLICIT** | single microstate |
| vatalanib | 2 | **ROBUST** | both tautomers coordinate, spread 0.131 Å |
| taranabant | 2 | **ROBUST** | both tautomers coordinate, spread 0.087 Å |
| nolatrexed | 6 | **ROBUST** | all 6 tautomers 5/5, spread 0.114 Å |
| olprinone | 6 | **ROBUST** | all 6 tautomers ≥4/5, spread 0.207 Å |
| R-1479 | 3 | **MICROSTATE-DEPENDENT** | one tautomer 5/5, the docked reference only 2/5, a third 0/5 |
| LDN-27219 | 3 | **MICROSTATE-DEPENDENT** | only the reference tautomer coordinates (4/5); the other two 0/5 |
| flucloxacillin | 16 | **MICROSTATE-DEPENDENT** | no microstate reaches 3/5 on re-embed (reference 1/5) |

### Controls (known azoles — validate the axis, not discoveries)

| control | states | verdict | note |
|---|:--:|---|---|
| flutrimazole | 1 | **EXPLICIT** | single microstate |
| letrozole | 1 | **EXPLICIT** | single microstate |
| fluconazole | 1 | **EXPLICIT** | single microstate |
| ravuconazole | 3 | **ROBUST** | all 3 tautomers 5/5, spread 0.176 Å |
| voriconazole | 4 | **MICROSTATE-DEPENDENT** | 2 tautomers coordinate 5/5, 2 fail (1/5, 2/5) |
| ketoconazole | 2 | **MICROSTATE-DEPENDENT** | both states 1/5 — its known floppy-ligand non-convergence |

## What this shows

- **The two zero-concern hits survive the axis without a dock.** **lersivirine and L-838417**
  are both protomer-explicit — no group titrates across pH 6.4–8.4 and each has a single
  tautomer, so the rank the scorecard puts forward refers to the one relevant microstate at
  physiological pH. The specific failure #84 names does not touch them.
- **The axis is not vacuous — it demotes two novel candidates on a genuine tautomer dependence.**
  **LDN-27219** coordinates the iron *only* in its reference tautomer (4/5 seeds); its two
  alternative tautomers never reach it (0/5). **R-1479** is worse-behaved still: one enumerated
  tautomer coordinates 5/5 while the docked reference manages only 2/5 and a third fails outright
  — the geometry that shortlisted it lives on a *different* tautomer than the one docked. In both
  cases *which* mobile-hydrogen position is present decides whether the molecule looks like a
  coordinator, exactly the risk this axis exists to catch. Each gains the concern "rank depends on
  protonation/tautomer state at pH 7.4."
- **flucloxacillin exposes a conformer-fragility the geometry axis hid.** Re-embedding its 16
  microstates and re-docking, **none reach the 3/5 coordination bar** (the reference protonation
  only 1/5) — its scorecard 5/5 at 2.517 Å did not reproduce on a fresh single-conformer embed.
  flucloxacillin is a flexible β-lactam; the honest reading is that its coordinating pose is
  conformer- and tautomer-fragile rather than that a single alternative tautomer breaks it. It
  already carried 3 concerns (structural alert, a prior *cell*-active measurement, coordinator-
  fragile domain); this reinforces caution and takes it to 4.
- **The controls exercise the axis as designed.** **ravuconazole is ROBUST** across its three
  tautomers (all 5/5, spread 0.176 Å). **voriconazole and ketoconazole are MICROSTATE-DEPENDENT**
  — voriconazole because two of its four tautomers fail to coordinate, ketoconazole because it is
  the known 1/5 floppy-ligand non-convergence already isolated as a search-power fault
  ([RESULTS_human_control.md](RESULTS_human_control.md)) and seen again on the stereo axis. Both
  are control-exempt: the point of a control is to work the machinery, not to be struck.
- **The scorecard gains a `protomer` column.** Novel candidates read `explicit`, `robust/Nst`, or
  `state-dep/Nst`. The MICROSTATE-DEPENDENT flag is a demotion for a *novel* candidate but
  control-exempt. Removing #84 from the `not-yet-checked` block leaves three open axes
  (ensemble/flexible #70, other auris mutants #77, broad CYP panel #74).

## What this establishes / does not

- **Does:** measure, per candidate, whether the iron-coordinating geometry survives the
  protonation and tautomer degrees of freedom the pipeline fixed silently; confirm the two
  zero-concern hits (lersivirine, L-838417) are single-microstate at pH 7.4; surface a real
  tautomer dependence in LDN-27219 and R-1479 and a conformer fragility in flucloxacillin; and
  demonstrate on the controls that the microstate-docking machinery separates a robust molecule
  (ravuconazole) from state-dependent ones (voriconazole, ketoconazole).
- **Does not:** second-guess Open Babel's pKa model (it tests the microstates obabel itself flips
  across pH 6.4–8.4, not whether that model is right); weight the enumerated tautomers by
  population (a MICROSTATE-DEPENDENT verdict may be driven by a tautomer that is not meaningfully
  populated — the diverging microstate's SMILES is recorded so that judgement can be made
  downstream); or relax any caveat of the underlying screen (single conformer per microstate,
  rigid receptor, coordination necessary-not-sufficient). Stereochemistry is the separate,
  already-measured axis #84's twin (#85).

## Reproduce

    conda run -n openafr python scripts/prep_receptor.py               # -> work/receptor.pdbqt
    conda run -n openafr python scripts/enumerate_protomers.py         # audit + microstate SDFs
    conda run -n openafr bash -c 'RECEPTOR=work/receptor.pdbqt \
        scripts/screen_consensus.sh work/candidates/protomer_ligands \
        work/candidates/dock_protomer "42 1 2 3 4"'                    # dock multi-state molecules
    conda run -n openafr python scripts/candidate_protomer.py --json \
        > work/candidates/shortlist_protomer.json                      # grade
    conda run -n openafr python work/candidates/build_scorecard.py     # fold into the scorecard

The enumeration + grader are committed; the per-seed dock tree
(`work/candidates/dock_protomer/`) is regenerated by the screen and gitignored like the other
dock trees.
