# Consolidated per-candidate confidence scorecard (#88)

Date: 2026-08-27 (resistance-retention axis added same day —
[RESULTS_candidate_resistance.md](RESULTS_candidate_resistance.md), #69/#77)

The capstone of the candidate-confidence batch (#65–#87): one honest, at-a-glance
card per candidate that rolls **every axis already measured** into the single
artifact you actually reason over before choosing what to send to a wet lab.
Confidence had been living across eight separate outputs — geometry, seed
stability, human-CYP51 selectivity, ADMET, synthesizability, applicability-domain
fragility, assay precedent — so a molecule could look good because we had only
looked at one axis. This forces all of them onto the same row.

> **No composite score.** A single blended number is exactly what this card must
> *not* produce — it would hide the one thing the scorecard exists to show: where a
> candidate is strong versus unproven. Each axis is reported on its own, and every
> axis the project has **not** yet measured is printed as `not-yet-checked`, never
> passed over in silence, so a card can never look complete by omission.

Builder: [`work/candidates/build_scorecard.py`](candidates/build_scorecard.py) —
a pure merge of frozen per-axis outputs; it re-docks and re-measures nothing.
Data: [`work/candidates/shortlist_scorecard.json`](candidates/shortlist_scorecard.json).

## The scorecard

16 molecules: 10 novel repurposing hypotheses (geometry-ranked) + 6 known-azole
controls that validate the axes rather than test a discovery.

### Novel candidates

| candidate | N–Fe (Å, 5-seed) | seeds | sel-margin (Å) | ADMET | SA band | domain | Y132F | precedent | concerns |
|---|---:|:--:|---:|---|---|---|---|---|:--:|
| verinurad | 2.417 | 5/5 | +0.51 | clean | easy | coordinator-fragile | RETAINED | no-precedent | 1 |
| vatalanib | 2.479 | 5/5 | +1.82 | hERG? | easy | clean | **LOST** | no-precedent | 2 |
| flucloxacillin | 2.517 | 5/5 | +2.66 | alert | moderate | coordinator-fragile | RETAINED | cell tested-active | 3 |
| taranabant | 2.564 | 5/5 | +2.05 | Ro5×2, insoluble | easy | coordinator-fragile | RETAINED | no-precedent | 2 |
| lersivirine | 2.572 | 5/5 | +2.07 | clean | easy | clean | RETAINED | no-precedent | **0** |
| L-838417 | 2.580 | 5/5 | +0.36 | clean | easy | clean | RETAINED | no-precedent | **0** |
| LDN-27219 | 2.603 | **2/5** | +0.64 | alert | easy | clean | **LOST** | no-precedent | 3 |
| R-1479 | 2.621 | 5/5 | +0.46 | Veber, PAINS, alert | moderate | clean | RETAINED | no-precedent | 1 |
| nolatrexed | 2.677 | 5/5 | **−0.06** | clean | easy | clean | RETAINED | no-precedent | 1 |
| olprinone | 2.680 | **4/5** | +1.32 | clean | easy | clean | RETAINED | no-precedent | 1 |

### Controls (known azoles — should look like azoles, not like discoveries)

| control | N–Fe (Å) | seeds | sel-margin | ADMET | Y132F | precedent | concerns |
|---|---:|:--:|---:|---|---|---|:--:|
| flutrimazole | 2.371 | 5/5 | −0.06 | alert | RETAINED | no-precedent | 2 |
| ravuconazole | 2.653 | 5/5 | +2.20 | clean | RETAINED | record-inconclusive | 0 |
| voriconazole | 2.657 | 5/5 | +0.67 | clean | RETAINED | **enzyme tested-active** | 1 |
| letrozole | 2.705 | 5/5 | −0.15 | clean | RETAINED | record-inconclusive | 1 |
| fluconazole | 2.711 | 5/5 | +0.06 | clean | RETAINED | **enzyme tested-active** | 1 |
| ketoconazole | 2.921 | **1/5** | +0.75 | PAINS, hERG? | RETAINED | — | 2 |

> **Controls RETAINING is the honesty check, not a resistance claim.** These azoles are
> clinically defeated by Y132F, yet they still coordinate the mutant iron — because the
> criterion is nearly insensitive to Y132F ([RESULTS_auris_Y132F.md](RESULTS_auris_Y132F.md)).
> So RETAINED means "geometry not destroyed," never "beats resistance." See the resistance
> results doc for the full caveat.

## What the aggregated view shows that no single axis did

- **Two candidates are clean on every axis we have measured: lersivirine and
  L-838417** — they coordinate the fungal iron reproducibly (5/5 seeds), open a
  positive human-selectivity window, carry no ADMET alert, sit in the easy SA band,
  are neither size- nor coordinator-fragile in the domain, and **both RETAIN
  coordination in the *C. auris* Y132F resistance pocket** (5/5 seeds there too). They
  are the two the scorecard puts forward with zero open concerns — *within the measured
  axes*.
- **The resistance axis surfaced a demotion no other axis caught: vatalanib.** It is
  5/5 seed-stable, one of the tightest coordinators (2.479 Å) and has the second-widest
  selectivity window (+1.82 Å) on the *C. albicans* stand-in — yet it collapses to
  **2/5 coordinating in the Y132F pocket (LOST)**. The geometry that shortlisted it is
  fragile in the resistant target the mission actually cares about; only docking the
  candidate against Y132F (#69/#77) exposed it. (Read the honest caveat: RETAINED ≠
  resistance-breaking; the informative direction of this axis is the LOST verdict.)
- **Selectivity and geometry can disagree with each other.** vatalanib,
  flucloxacillin and taranabant have the widest human-selectivity margins
  (+1.8 to +2.7 Å) yet each carries a separate liability the margin hides —
  vatalanib a hERG risk factor, flucloxacillin a prior *cell*-active measurement
  and a structural alert, taranabant two Ro5 violations plus predicted insolubility.
  A shortlist ranked on selectivity alone would have promoted all three unqualified.
- **The instability demotions survive aggregation.** LDN-27219 (2/5) and olprinone
  (4/5) remain the only sub-5/5 novel candidates — their single-dock ranks were
  search luck, consistent with the seed-stability finding (#68).
- **nolatrexed is the one geometrically-real, seed-stable hit that fails
  selectivity** (−0.06 Å: reaches the human iron as readily as the fungal one) — a
  demotion only the human counter-screen (#73) could produce, now visible on the
  same row as its clean ADMET/domain profile.
- **The controls behave as a negative control on the scorecard itself.**
  ketoconazole surfaces its known 1/5 non-convergence (a floppy-ligand global-search
  artifact already isolated as a search-power, not receptor, fault in
  [RESULTS_human_control.md](RESULTS_human_control.md)); fluconazole/voriconazole
  correctly carry the only `enzyme tested-active` precedent flags. A scorecard that
  hid these would be the broken one.

## Axes deliberately left `not-yet-checked` (open issues)

Printed on every card so no molecule looks fully vetted:

- **Ensemble / flexible-receptor docking on the candidate path** (#70) — the
  ensemble was validated on the held-out set, never per candidate.
- **Cluster-spread / RMSD pose convergence** (#71) — seed spread is a proxy; true
  pose-cluster RMSD is unmeasured.
- **auris K143R / F126L mutant retention** (#77) — the dominant Y132F substitution is
  now measured (the `Y132F` column above;
  [RESULTS_candidate_resistance.md](RESULTS_candidate_resistance.md)), but the clade-I
  K143R and clade-III F126L mutants need a side-chain repacker not in the pinned
  environment and remain unchecked.
- **Protonation / tautomer enumeration at pH 7.4** (#84) and **all-stereoisomer
  docking** (#85) — each candidate was docked as one supplied SMILES.
- **Broader human CYP off-target panel** (CYP3A4/2C9/2D6/2C19/1A2, #74) — only the
  CYP51 counter-screen exists.

## What this establishes / does not

- **Does:** produce the decision artifact — every measured confidence axis on one
  row per candidate, with honest `not-yet-checked` states — and confirm that when
  all axes are seen together, **lersivirine and L-838417** are the two novel hits
  with no open concern on any measured axis, while the geometry order alone would
  have mixed them in with seed-unstable (LDN-27219, olprinone), non-selective
  (nolatrexed), ADMET-flagged (R-1479, taranabant) and resistance-fragile (vatalanib)
  molecules.
- **Does not:** call any molecule a hit, collapse the axes into a rank, or imply the
  `not-yet-checked` axes would pass. Every caveat of the underlying screen still
  binds (single rigid *C. albicans* conformer; coordination necessary, not
  sufficient; `no-precedent` ≠ untested).

## Reproduce

    conda run -n openafr python work/candidates/build_scorecard.py

Reads only committed per-axis outputs (`shortlist_confidence.json`,
`shortlist_admet.tsv`, `shortlist_synth.json`, `domain_sensitivity.json`,
`top100_precedent.tsv`, `shortlist_resistance.json`, `work/repurposing/ranked.tsv`);
writes `work/candidates/shortlist_scorecard.json`. No network, no docking.
