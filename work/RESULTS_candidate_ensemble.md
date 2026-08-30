# Results — candidate ensemble-consistency (*C. albicans* CYP51 crystallographic ensemble)

Date: 2026-08-29. Pre-registered in
[PREREGISTRATION_candidate_ensemble.md](PREREGISTRATION_candidate_ensemble.md)
(sha256 65c3db8f3283…, frozen BEFORE any 5FSA/5V5Z distance was read;
`work/PREREG_candidate_ensemble.sha256`). Protocol hash unchanged from run 2.

> **This is a per-candidate confidence axis, not a gate.** It demotes/annotates the
> scorecard; it never promotes. It closes the first `not-yet-checked` line the scorecard
> printed on every card (#70): are the shortlist candidates' iron-coordinating geometries a
> property of the single 5TZ1 conformer they were ranked in, or do they hold across every
> *C. albicans* CYP51 channel conformation we have a crystal for?

## What was tested

The 16 scorecard molecules (10 novel candidates + 6 known-azole controls) were docked into
the **exhaustive crystallographic *C. albicans* CYP51 ensemble** — the three and only three
such structures in the PDB: **5TZ1** (the project anchor), **5FSA** (posaconazole co-crystal),
**5V5Z** (itraconazole co-crystal). The two new conformers were prepared by
`scripts/prep_receptor_ensemble.py` (Kabsch-superposed onto the 5TZ1 frame so the frozen box
lands in each pocket; each N–Fe measured against that conformer's own iron) and docked under
the frozen protocol, 5 seeds {42,1,2,3,4}. The 5TZ1 values are the already-frozen
`fungal_consensus` in `work/candidates/shortlist_confidence.json` (**not re-docked**).

Per molecule: `ensemble_min` = MIN consensus over conformers where it had any coordinating
pose; `n_conf_coord` = how many conformers it coordinates in (≥3/5 seeds each). Verdicts
(frozen): **ROBUST** = 3/3 conformers coordinate; **CONSISTENT** = 2/3; **CONFORMER-FRAGILE**
= exactly 1/3 (single-snapshot artifact); **NON-COORDINATING** = 0/3.

Self-docking guard applied but a no-op here — none of the 16 molecules is posaconazole or
itraconazole, so nothing was excluded from its own co-crystal.

## The result

| candidate | WT N–Fe (5TZ1) | 5FSA (seeds) | 5V5Z (seeds) | ens-min | conf | verdict |
|---|---:|---:|---:|---:|:--:|---|
| flutrimazole* | 2.371 (5/5) | 2.778 (5/5) | 2.759 (5/5) | 2.371 | 3/3 | ROBUST |
| verinurad | 2.417 (5/5) | 2.727 (5/5) | 2.667 (4/5) | 2.417 | 3/3 | ROBUST |
| **vatalanib** | 2.479 (5/5) | 2.725 (5/5) | 2.811 (5/5) | 2.479 | 3/3 | **ROBUST** |
| flucloxacillin | 2.517 (5/5) | 2.920 (5/5) | — (0/5) | 2.517 | 2/3 | CONSISTENT |
| taranabant | 2.564 (5/5) | 2.717 (5/5) | 2.848 (5/5) | 2.564 | 3/3 | ROBUST |
| **lersivirine** | 2.572 (5/5) | 2.720 (5/5) | 2.596 (5/5) | 2.572 | 3/3 | **ROBUST** |
| L-838417 | 2.580 (5/5) | — (0/5) | 2.965 (3/5) | 2.580 | 2/3 | CONSISTENT |
| **LDN-27219** | 2.603 (2/5) | — (0/5) | 2.944 (4/5) | 2.603 | **1/3** | **CONFORMER-FRAGILE** |
| R-1479 | 2.621 (5/5) | 2.794 (5/5) | 2.814 (3/5) | 2.621 | 3/3 | ROBUST |
| ravuconazole* | 2.653 (5/5) | 2.697 (5/5) | 2.953 (3/5) | 2.653 | 3/3 | ROBUST |
| voriconazole* | 2.657 (5/5) | 2.722 (5/5) | 2.731 (5/5) | 2.657 | 3/3 | ROBUST |
| nolatrexed | 2.677 (5/5) | 2.741 (5/5) | 2.833 (5/5) | 2.677 | 3/3 | ROBUST |
| **olprinone** | 2.680 (4/5) | — (0/5) | — (0/5) | 2.680 | **1/3** | **CONFORMER-FRAGILE** |
| letrozole* | 2.705 (5/5) | 2.809 (5/5) | 2.820 (5/5) | 2.705 | 3/3 | ROBUST |
| fluconazole* | 2.711 (5/5) | 2.785 (5/5) | 2.802 (5/5) | 2.711 | 3/3 | ROBUST |
| ketoconazole* | 2.921 (1/5) | 2.686 (4/5) | 2.797 (4/5) | 2.686 | 2/3 | CONSISTENT |

*controls (known azoles). Tally: **11 ROBUST, 3 CONSISTENT, 2 CONFORMER-FRAGILE, 0 NON-COORDINATING.**
(5FSA/5V5Z distances shown are the 5-seed consensus min; "—" = no conformer seed reached the iron.)

## What this shows that no earlier axis did

- **olprinone is a genuinely new demotion.** It coordinated the 5TZ1 iron in 4/5 seeds
  (2.680 Å) and RETAINED on Y132F — but on **both** other crystal conformers it never reaches
  the iron (0/5 on 5FSA *and* 5V5Z). Its coordination is an artifact of the single snapshot it
  was ranked in; **no prior axis — geometry, seed stability, selectivity, ADMET, Y132F —
  surfaced this.** CONFORMER-FRAGILE, one new scorecard concern (now 2).
- **LDN-27219's demotion piles on, consistently.** Already seed-unstable (2/5 on 5TZ1) and
  Y132F-LOST, it now coordinates reliably in only 5V5Z (4/5; 5TZ1 2/5, 5FSA 0/5) →
  CONFORMER-FRAGILE. Same molecule failing the same way on yet another axis — not a surprise,
  a confirmation.
- **lersivirine is now the single strongest clean hit.** Of the two novel candidates that were
  clean on every prior axis, **lersivirine is ROBUST (3/3, all seeds in every conformer)**
  while **L-838417 is CONSISTENT (2/3 — it never coordinates 5FSA, 0/5)**. CONSISTENT is not a
  demotion (both keep zero concerns), but the ensemble is the first axis to *separate* the two
  former co-front-runners: lersivirine's coordination holds across every experimental
  conformation, L-838417's holds in two of three.
- **Conformer-robustness and resistance-robustness are orthogonal — vatalanib proves it.**
  vatalanib is **ROBUST across the WT crystal ensemble (3/3)** yet **LOST in the Y132F pocket**.
  A molecule can be insensitive to which wild-type conformation you dock it against and *still*
  be fragile to the resistance substitution the mission targets. The two axes measure different
  things; neither subsumes the other.
- **The ensemble does not collapse the shortlist.** 11/16 molecules are ROBUST and none is
  NON-COORDINATING — consistent with the held-out finding that the geometric criterion is only
  weakly conformer-sensitive ([RESULTS_ensemble.md](RESULTS_ensemble.md)). The axis is a
  single-snapshot-artifact *filter*, and it fired on exactly the two molecules that were already
  the shakiest coordinators, plus sharpened one (olprinone) that had looked fine on geometry.

## Why ROBUST is not a potency claim (do NOT overclaim)

The within-azole ceiling work already ran this exact ensemble and showed it does **not** lift
ranking of working azoles above failed ones to the 0.70 bar (look #9,
[RESULTS_ensemble_confirm.md](RESULTS_ensemble_confirm.md), ensemble-min AUC 0.682). So a
**ROBUST** verdict means only that a molecule's coordination geometry is reproduced across
every experimentally-observed *C. albicans* channel conformation — a robustness/artifact check,
**not** evidence it binds more tightly or defeats any azole. The known-azole controls being
mostly ROBUST is the clearest proof: those are the drugs *C. auris* already defeats, yet their
coordination is conformer-robust. **The informative direction of this axis is the
CONFORMER-FRAGILE verdict**, which flags a shortlist molecule resting on one crystal snapshot.

- **ketoconazole's ensemble "improvement" is not signal.** Its ens-min (2.686 Å) is tighter
  than its 5TZ1 value (2.921 Å) because the 5TZ1 baseline was a single poorly-converged seed
  (1/5, a known floppy-ligand global-search artifact —
  [RESULTS_human_control.md](RESULTS_human_control.md)); the other two conformers simply recover
  the coordination that one bad search on 5TZ1 missed (4/5 each). Reported, not leaned on —
  and a useful demonstration that the ensemble smooths a per-conformer search artifact rather
  than inventing a binding claim.

## Limitations (all from the pre-registration, all still bind)

1. **Crystallographic conformers only.** Three rigid experimental structures; MD-relaxed /
   induced-fit flexibility (mentioned in #70) needs an MD stack not in the pinned environment
   and remains out of scope. This is a rigid-ensemble robustness check, not full receptor
   flexibility.
2. **Three conformers is a coarse ensemble** — the verdict is one of {0,1,2,3}/3; it separates
   "single-snapshot" from "generalizes," not fine gradations.
3. ***C. albicans*, not *C. auris*.** These are *C. albicans* structures (no experimental *C.
   auris* CYP51 exists); resistance-pocket robustness is the separate Y132F axis
   ([RESULTS_candidate_resistance.md](RESULTS_candidate_resistance.md)).
4. **Coordination necessary-not-sufficient; single rigid pose per conformer.** Every caveat of
   the underlying screen still binds.

## Reproduce

```
conda activate openafr
python scripts/prep_receptor_ensemble.py 5FSA 5V5Z
python scripts/prep_ligands.py work/candidates/shortlist_focus.smi work/candidates/cand_ligands
RECEPTOR=work/receptor_5FSA.pdbqt scripts/screen_consensus.sh \
    work/candidates/cand_ligands work/candidates/dock_ens_5FSA
RECEPTOR=work/receptor_5V5Z.pdbqt scripts/screen_consensus.sh \
    work/candidates/cand_ligands work/candidates/dock_ens_5V5Z
python scripts/candidate_ensemble.py --json > work/candidates/shortlist_ensemble.json
python scripts/candidate_ensemble.py    # human-readable table
python work/candidates/build_scorecard.py   # fold into the consolidated scorecard
```
