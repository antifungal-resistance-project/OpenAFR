# Ligand set provenance + active-set expansion spec (T1)

This file documents how the active sets were built and specifies the disciplined
expansion called for by the 2026-08-21 eng review + outside voice
(`docs/designs/refocus-resistant-cyp51-triage.md`, task T1).

## Current sets (as of 2026-08-21)

### Training actives — `actives.smi` (n=8)
Classic azole antifungals, used to define the geometric criterion:
fluconazole, voriconazole, itraconazole, posaconazole, ketoconazole, clotrimazole,
miconazole, isavuconazole.

### Held-out actives — `actives_holdout_final.smi` (n=7)
Blind test set, fixed in `work/PREREGISTRATION_run2.md` and reused unchanged by the
Y132F transfer test (`work/PREREGISTRATION_auris_Y132F.md`). Ten azoles were retrieved;
three (econazole 0.79, tioconazole 0.76, terconazole) were removed as near-duplicates of
training actives (Morgan r=2 Tanimoto ≥ 0.70). Remaining 7:
efinaconazole, luliconazole, sertaconazole, oxiconazole, bifonazole, fenticonazole,
butoconazole. (`actives_holdout.smi` is the pre-filter list of 10.)

### Decoys — `decoys.smi`, `decoys_holdout.smi`
Property-matched presumed-inactives (see `scripts/make_decoys.py`).

## The problem T1 fixes

Two independent findings converged on the active set:

1. **Every active is a defeated azole.** Re-discovering drugs that resistant *C. auris*
   already beats does not demonstrate the tool can find resistance-*breakers* — the actual
   goal of any wet-lab shortlist. Mutant AUC 0.805 ≈ WT 0.794: the criterion is nearly
   *insensitive* to the resistance mutation, so "broad enrichment transfers to Y132F" is
   close to a tautology.
2. **n=7 is too thin.** EF@5% on 7 actives / 350 decoys counts ~18 molecules; one active
   moves it ~1.4x. AUC and the permutation test have large variance at 7 positives.

## Expansion design (differential-activity, not "more azoles")

The discriminating test: add CYP51 inhibitors that **retain measured activity against
azole-resistant *C. auris*** (the newer tetrazoles), so a genuinely resistance-aware
criterion should rank the still-active compounds **above** the defeated azoles *in the
mutant pocket*. That converts the validation from "rank azoles vs decoys" (tautological on
the mutant) to "rank still-active vs defeated, in the resistant pocket" (a real claim).

**Candidate additions — SMILES VERIFIED from PubChem 2026-08-21 (isomeric/stereo form,
copied from the authoritative record, not AI-generated):**

- **Oteseconazole (VT-1161)** — PubChem **CID 77050711**, C23H16F7N5O2, MW 527.4.
  `C1=CC(=CC=C1C2=CN=C(C=C2)C([C@](CN3C=NN=N3)(C4=C(C=C(C=C4)F)F)O)(F)F)OCC(F)(F)F`
- **VT-1598 (quilseconazole)** — PubChem **CID 91886002**, C22H14F7N5O2, MW 513.4.
  `C1=CC(=CC=C1C2=CN=C(C=C2)C([C@](CN3C=NN=N3)(C4=C(C=C(C=C4)F)F)O)(F)F)OC(F)(F)F`

Both carry a **tetrazole** metal-binding group (`CN3C=NN=N3`) rather than the
imidazole/triazole of every current active — genuinely distinct chemistry that should clear
the Morgan-Tanimoto novelty filter. The two differ only at the ether tail (oteseconazole
–OCH2CF3 vs VT-1598 –OCF3).

- (Optional) additional agents only if a primary-literature MIC vs a resistant *C. auris*
  isolate exists.

**Activity evidence vs *C. auris* (literature search 2026-08-21 — verify numbers against the
primary papers before freezing):**

- **VT-1598 — SUPPORTED. Use this one.** Wiederhold et al., "The Fungal Cyp51-Specific
  Inhibitor VT-1598 Demonstrates In Vitro and In Vivo Activity against *Candida auris*,"
  *Antimicrob Agents Chemother* 2018, PMID 30530603 (doi 10.1128/aac.02233-18). Reported
  modal MIC **0.25 µg/mL**, range **0.03–8 µg/mL** across 100 clinical *C. auris* isolates,
  activity against azole-resistant isolates + in vivo burden reduction. This is the
  differential-activity evidence the design needs. (Numbers are from the abstract/summary;
  confirm exact values in the paper before the pre-registration freezes.)
- **Oteseconazole (VT-1161) — NOT SUPPORTED for this use; DROP.** Its strong MIC data is vs
  *C. albicans*/*C. glabrata* (incl. fluconazole-resistant) from vulvovaginal-candidiasis
  studies; no *C. auris*-specific MIC surfaced. Including it as a "resistance-retained" active
  would be the exact unsupported assumption this design guards against. Structure stays on
  record (CID 77050711) but it is NOT added as an active without a real *C. auris* source.

**Net:** add **VT-1598 only** as the resistance-retained active. (It also resolves the
0.823 self-similarity issue — with oteseconazole dropped there is no near-duplicate pair.)
One tetrazole is thin on its own for lifting n; consider additional distinct-scaffold
resistance-active chemotypes (e.g. a non-CYP51 agent as a positive-control active only if
the geometric criterion is not expected to rank it — decide deliberately in the prereg).

## Discipline constraints (do NOT shortcut)

1. **New pre-registration, not an edit of the frozen holdout.** The existing holdout is
   sha-pinned and already graded. Any expanded set defines a NEW pre-registered run
   (`work/PREREGISTRATION_run3.md` or similar) with its own frozen list + hash. Never add
   compounds to `actives_holdout_final.smi` and re-grade an existing run — that is the
   p-hacking move D5 already rejected.
2. **Novelty filter vs BOTH sets.** Each new active must pass the same Morgan r=2 Tanimoto
   < 0.70 filter against training AND the existing holdout, so it is genuinely new chemistry.
3. **Verified inputs only.** Each new compound needs: a PubChem CID, the canonical isomeric
   SMILES copied from that authoritative record (NOT AI-generated), and a primary-literature
   citation for its activity against a resistant strain. Record CID + citation inline here.
   A wrong SMILES silently poisons the run — treat this like the receptor hash gate.
4. **Docking is gated on the cloud Linux/x86 VM** (vina/obabel absent on osx-arm64), same as
   every other screen. Curation + pre-registration here are local; the run is not.

## Status

- [x] Current sets documented (this file).
- [x] Verified + recorded CID + isomeric SMILES for oteseconazole (CID 77050711) and
      VT-1598 (CID 91886002) from PubChem.
- [x] Literature check for MIC vs resistant *C. auris*: VT-1598 SUPPORTED (PMID 30530603,
      modal MIC 0.25 µg/mL); oteseconazole NOT supported for *C. auris* → dropped. Add
      VT-1598 only. (Confirm exact MICs against the paper before freezing the prereg.)
- [ ] Write the expanded-run pre-registration with the differential-activity hypothesis and bar.
- [x] Novelty-filtered vs the 15 existing actives (Morgan r=2, 2048-bit, run under
      `conda run -n openafr`, rdkit 2025.09.5, 2026-08-21):
      oteseconazole max Tanimoto **0.333** (nearest fluconazole); VT-1598 **0.358**
      (nearest fluconazole). Both PASS (< 0.70) — genuinely novel vs the azole panel.
      **Caveat:** oteseconazole vs VT-1598 = **0.823** (same tetrazole scaffold), so by the
      project's own 0.70 rule the two are near-duplicates of EACH OTHER — put one in the
      holdout and the other in training, not both in the same blind set.
- [ ] Dock + grade on the cloud VM (AUC/EF@5%; do NOT lean on EF@1%).
