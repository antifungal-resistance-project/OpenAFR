# Results — candidate resistance-pocket retention, add-atom mutants (*C. auris* K143R, F126L)

Date: 2026-08-30. Pre-registered in
[PREREGISTRATION_candidate_resistance_addatom.md](PREREGISTRATION_candidate_resistance_addatom.md)
(sha256 913cde54…, frozen — with both mutant receptor shas — BEFORE any K143R/F126L distance
was read; `work/PREREG_candidate_resistance_addatom.sha256`). Protocol hash unchanged from run 2.
The direct sibling of the Y132F run ([RESULTS_candidate_resistance.md](RESULTS_candidate_resistance.md));
read it first — the criterion, bar and interpretation are inherited unchanged.

> **This is a per-candidate confidence axis, not a gate.** It demotes/annotates the scorecard;
> it never promotes. It closes the last two mission-shaped holes on every card (#77): do the
> shortlist candidates keep their iron-coordinating geometry in the two *clade-specific*
> resistance pockets — **K143R** (clade I) and **F126L** (clade III) — or only in the
> *C. albicans* stand-in and the dominant Y132F pocket already tested?

## What is new versus the Y132F run — a pinned, honestly-labeled repacker

Y132F is a pure para-hydroxyl **deletion**: its mutant side chain is a strict subset of the
wild-type atoms, so `scripts/mutate_receptor.py` builds it with **zero modeling**. K143R
(Lys→Arg) and F126L (Phe→Leu) both **add** side-chain atoms, so they cannot be built that way
— they need a side-chain repacker, which is exactly why they were deferred. This run introduces
one, declared in full:

- **`scripts/repack_mutant.py`** + a pinned `pymol-open-source=3.1.0` in the separate build-time
  env `environment-repack.yml` (kept out of the main env, which it does not co-solve with; the
  repacker is build-time only). The new side chain is placed by PyMOL's mutagenesis wizard,
  selecting the **minimum-strain**
  (fewest-clash) backbone-dependent rotamer — one deterministic choice, re-running yields a
  byte-identical receptor PDB.
- **The mutation stays the only variable.** The build is refused unless every atom outside the
  mutated residue is bit-identical to `work/receptor_A.pdb` and the heme iron is unmoved
  (observed: worst non-mutated-atom delta **0.000 Å**, Fe fixed at 64.163/71.316/2.711).
- **Labeled MODELED-ROTAMER everywhere it surfaces** (scorecard column, concern text, JSON),
  never conflated with the exact-geometry Y132F/V125A deletions. It is a **weaker structural
  claim** — see the honesty caveat below.

## What was tested

The 16 scorecard molecules (10 novel candidates + 6 known-azole controls) were docked into the
**K143R** and **F126L** *C. auris* pockets (`work/receptor_auris_{K143R,F126L}.pdbqt`, each a
single repacked side chain on the validated 5TZ1 chain-A scaffold — same heme, same iron, same
box) under the frozen protocol, 5 seeds {42,1,2,3,4}. The **only** variable versus the wild-type
candidate dock is the one mutated side chain. The wild-type baseline is the already-frozen
`fungal_consensus` in `work/candidates/shortlist_confidence.json` (not re-docked). Verdicts are
the frozen `classify` rule reused verbatim: **RETAINED** = coordinates the mutant iron ≥3/5 seeds
AND |shift| ≤ 0.5 Å; **SHIFTED** = ≥3/5 but |shift| > 0.5 Å; **LOST** = <3/5 seeds reach the iron.

## The result

### K143R (clade I) — Lys143→Arg

| candidate | WT N–Fe | K143R N–Fe | coord | shift (Å) | verdict |
|---|---:|---:|:--:|---:|---|
| flutrimazole* | 2.371 | 2.362 | 5/5 | −0.009 | RETAINED |
| verinurad | 2.417 | 2.425 | 5/5 | +0.008 | RETAINED |
| vatalanib | 2.479 | 2.436 | 5/5 | −0.043 | RETAINED |
| flucloxacillin | 2.517 | 2.540 | 5/5 | +0.023 | RETAINED |
| taranabant | 2.564 | 2.526 | 5/5 | −0.038 | RETAINED |
| lersivirine | 2.572 | 2.562 | 5/5 | −0.010 | RETAINED |
| L-838417 | 2.580 | 2.558 | 5/5 | −0.022 | RETAINED |
| LDN-27219 | 2.603 | 2.610 | 3/5 | +0.007 | RETAINED |
| R-1479 | 2.621 | 2.610 | 5/5 | −0.011 | RETAINED |
| ravuconazole* | 2.653 | 2.786 | 3/5 | +0.133 | RETAINED |
| voriconazole* | 2.657 | 2.689 | 5/5 | +0.032 | RETAINED |
| nolatrexed | 2.677 | 2.691 | 5/5 | +0.014 | RETAINED |
| olprinone | 2.680 | 2.687 | 5/5 | +0.007 | RETAINED |
| letrozole* | 2.705 | 2.719 | 5/5 | +0.014 | RETAINED |
| fluconazole* | 2.711 | 2.688 | 3/5 | −0.023 | RETAINED |
| ketoconazole* | 2.921 | 2.815 | 1/5 | — | LOST |

*controls. Tally: **15 RETAINED, 1 LOST** — and the single LOST is ketoconazole, whose WT
baseline was already a 1/5 floppy-ligand global-search artifact (see the Y132F run and
[RESULTS_human_control.md](RESULTS_human_control.md)), not a receptor effect. **No novel
candidate is demoted by K143R.** The clade-I substitution barely perturbs the criterion
(every RETAINED molecule shifts <0.15 Å) — the same near-insensitivity the Y132F run reported.

### F126L (clade III) — Phe126→Leu

| candidate | WT N–Fe | F126L N–Fe | coord | shift (Å) | verdict |
|---|---:|---:|:--:|---:|---|
| flutrimazole* | 2.371 | 2.498 | 5/5 | +0.127 | RETAINED |
| verinurad | 2.417 | 2.833 | 5/5 | +0.416 | RETAINED |
| vatalanib | 2.479 | 2.508 | 5/5 | +0.029 | RETAINED |
| **flucloxacillin** | 2.517 | 4.077 | **0/5** | — | **LOST** |
| taranabant | 2.564 | 2.616 | 5/5 | +0.052 | RETAINED |
| lersivirine | 2.572 | 2.558 | 4/5 | −0.014 | RETAINED |
| L-838417 | 2.580 | 2.609 | 5/5 | +0.029 | RETAINED |
| LDN-27219 | 2.603 | 2.592 | **1/5** | — | **LOST** |
| R-1479 | 2.621 | 2.696 | 5/5 | +0.075 | RETAINED |
| ravuconazole* | 2.653 | 2.864 | **2/5** | — | **LOST** |
| voriconazole* | 2.657 | 2.645 | 5/5 | −0.012 | RETAINED |
| nolatrexed | 2.677 | 2.689 | 5/5 | +0.012 | RETAINED |
| olprinone | 2.680 | 2.593 | 5/5 | −0.087 | RETAINED |
| letrozole* | 2.705 | 2.652 | 5/5 | −0.053 | RETAINED |
| fluconazole* | 2.711 | 2.686 | 5/5 | −0.025 | RETAINED |
| ketoconazole* | 2.921 | 2.824 | 5/5 | −0.097 | RETAINED |

*controls. Tally: **13 RETAINED, 3 LOST**. The clade-III pocket is the more discriminating of
the two — the bulky Phe126→Leu opening shifts several molecules by >0.1 Å (verinurad +0.416 Å,
the largest RETAINED move in the whole axis), and it produces the one genuinely new demotion.

## What this shows that no earlier axis did

- **flucloxacillin is a new, F126L-specific demotion.** It was RETAINED on Y132F (4/5) and on
  K143R (5/5), sat mid-shortlist on geometry, and passed seed-stability — yet in the **F126L**
  pocket it coordinates in **0/5** seeds, its closest nitrogen falling to 4.08 Å. The clade-III
  substitution destroys its coordination entirely. **No prior axis — geometry, seed stability,
  Y132F, K143R, ensemble — surfaced this**; it is exactly the failure the mission cares about
  (solid on the stand-in, fragile in a resistant target), and it is now a scorecard concern
  (flucloxacillin 3 → 4 concerns).
- **Resistance-robustness is mutation-specific, not a single property.** vatalanib is LOST on
  Y132F but RETAINED on both K143R and F126L; flucloxacillin is the mirror image (RETAINED on
  Y132F/K143R, LOST on F126L). A candidate can survive one clade's pocket and collapse in
  another — a molecule is only as good as its *worst* resistant pocket, which is precisely why
  the panel had to be completed rather than stopping at the dominant Y132F.
- **The two front-runners survive all three pockets.** lersivirine (Y132F 5/5, K143R 5/5, F126L
  4/5 — all RETAINED) and L-838417 (5/5 on all three) keep their zero/one-concern standing. The
  add-atom pockets open no new concern on either, so the scorecard's clean hits stay clean.
- **The other two F126L LOSTs are not new signal.** LDN-27219 was already seed-unstable
  everywhere (2/5 WT, 1/5 Y132F, 1/5 F126L) — same instability, not a resistance effect;
  ravuconazole is a control azole (2/5), exempt by construction.

## The honesty caveat (carried from the Y132F run, plus the modeling caveat — do NOT overclaim)

**RETAINED ≠ resistance-breaking.** As on Y132F, the criterion is nearly insensitive to these
pocket edits, so broad retention is close to a tautology — a RETAINED verdict means only that the
coordination geometry which shortlisted a molecule is **not destroyed** by the substitution, a
robustness check, never evidence a molecule *defeats* resistance. The known-azole controls
mostly RETAINING (they are clinically resistance-defeated) is the clearest proof. **The
informative direction is the LOST verdict** — and flucloxacillin/F126L is the new one.

**Additionally: these are MODELED-ROTAMER pockets, a strictly weaker structural claim than
Y132F.** The K143R/F126L side chains are placed by the repacker's minimum-strain rotamer, not
solved crystallographically. A different (still low-strain) rotamer could move a borderline
verdict; every verdict here is conditioned on the frozen rotamer choice. flucloxacillin's
collapse is large (2.52 → 4.08 Å, 5/5 → 0/5) and so is unlikely to be a rotamer artifact, but it
is reported as a modeled-pocket demotion, not an observed-structure fact.

## Limitations (all from the pre-registration, all still bind)

1. **Modeled rotamer, not observed geometry** (the new caveat above).
2. Rigid backbone, single conformer, only the one mutated side chain differs from wild type — no
   induced-fit repacking, no efflux/expression changes.
3. **Point substitutions only.** TR-type tandem-duplication / promoter (copy-number) resistance
   is an expression event, not a CDS pocket point mutation, and stays out of structural scope
   (printed as the one remaining `not-yet-checked` line on the scorecard).
4. **Controls are controls.** Ranking the known azoles well on a mutant pocket is a robustness
   demonstration, never a resistance claim.

## Reproduce

```
# Build the two add-atom mutants with the pinned repacker (deterministic min-strain rotamer).
# The repacker lives in its own build-time env (pymol does not co-solve with the main pins):
conda env create -f environment-repack.yml   # once
conda activate openafr-repack
pymol -cq scripts/repack_mutant.py -- K143R -o work/receptor_auris_K143R.pdb
pymol -cq scripts/repack_mutant.py -- F126L -o work/receptor_auris_F126L.pdb
# Everything else runs in the main env:
conda activate openafr
python scripts/prep_receptor.py --mutant work/receptor_auris_K143R.pdb work/receptor_auris_K143R.pdbqt
python scripts/prep_receptor.py --mutant work/receptor_auris_F126L.pdb work/receptor_auris_F126L.pdbqt
python scripts/prep_ligands.py work/candidates/shortlist_focus.smi work/candidates/cand_ligands
RECEPTOR=work/receptor_auris_K143R.pdbqt scripts/screen_consensus.sh \
    work/candidates/cand_ligands work/candidates/dock_auris_K143R
RECEPTOR=work/receptor_auris_F126L.pdbqt scripts/screen_consensus.sh \
    work/candidates/cand_ligands work/candidates/dock_auris_F126L
python scripts/candidate_resistance_addatom.py K143R --out work/candidates/shortlist_resistance_K143R.json
python scripts/candidate_resistance_addatom.py F126L --out work/candidates/shortlist_resistance_F126L.json
python work/candidates/build_scorecard.py --out work/candidates/shortlist_scorecard.json
```
