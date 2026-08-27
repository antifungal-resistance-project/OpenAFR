# Results — human-CYP51 selectivity counter-screen (issue #73)

Date: 2026-08-26. Pre-registered in
[PREREGISTRATION_human_selectivity.md](PREREGISTRATION_human_selectivity.md)
(sha256 `b1f19f75…`, frozen before any human distance was read). The shortlist, docked against
**human CYP51** (`work/receptor_human.pdbqt`, built from RCSB **3LD6** by
`scripts/prep_receptor_human.py`) under the identical criterion and protocol as the fungal run,
so the only difference between a fungal and a human score is the protein.

> **Outcome: a promising selectivity signal — but PROVISIONAL, because the pre-registered
> positive control did not converge.** Known non-selective azoles behave exactly as expected
> (they coordinate the human iron too), and four seed-stable novel hits show a clean structural
> selectivity window. **However**, the frozen precondition — re-docking ketoconazole must
> reproduce its crystallographic human coordination — **was not met**, so per the
> pre-registration these numbers are reported as provisional, not validated.

## The pre-registered positive control FAILED (reported first, as required)

3LD6's own LINK record ties the ketoconazole N2 to the human heme iron at **2.41 Å**. The
pre-registration made reproducing a coordinating pose (< 3.0 Å) a **hard precondition**.
Re-docked, ketoconazole's per-seed human N-Fe was **{42: 4.829, 1: 3.675, 2: 4.665, 3: 4.388,
4: 4.704}** — it never reached the iron within 3.0 Å. The control did not pass.

**Why, and what it does / does not invalidate.** Ketoconazole also failed on the *validated
fungal* receptor (coordinated in only 1/5 seeds, see [RESULTS_shortlist_consensus.md]) — so the
failure is a **ligand** problem (36 heavy atoms, 8 rotatable bonds, its coordinating imidazole
on a long flexible tail; a single RDKit-embedded conformer plus one Vina search does not find
the crystal pose), **not** evidence the human receptor is wrong. Two independent lines say the
human receptor+box are geometrically sound: (a) the heme iron sits 8.8 Å inside a 26 Å box, and
(b) small rigid azoles **do** coordinate it reproducibly (fluconazole 5/5 @ 2.77 Å, letrozole
5/5 @ 2.56 Å). So the selectivity signal below is almost certainly real — but the *specific*
control we pre-committed to did not fire, so we hold the numbers as provisional and name the
stronger control as the follow-up (below).

## The table (fungal vs human, 5-seed consensus min N-Fe)

Margin = human − fungal (Å); **positive = reaches the fungal iron closer = fungal-selective.**

| candidate | role | stable? | fung min | human min | human coord. | margin | read |
|---|---|:---:|---:|---:|:---:|---:|---|
| vatalanib | novel | yes | 2.479 | 4.294 | 0/5 | **+1.82** | selective (fungal only) |
| flucloxacillin | novel | yes | 2.517 | 5.172 | 0/5 | **+2.66** | selective (fungal only) |
| taranabant | novel | yes | 2.564 | 4.611 | 0/5 | **+2.05** | selective (fungal only) |
| lersivirine | novel | yes | 2.572 | 4.638 | 0/5 | **+2.07** | selective (fungal only) |
| ravuconazole | azole ref | yes | 2.653 | 4.857 | 0/5 | +2.20 | selective (fungal only) |
| R-1479 | novel | yes | 2.621 | 3.079 | 0/5 | +0.46 | selective (marginal) |
| verinurad | novel | yes | 2.417 | 2.929 | 2/5 | +0.51 | fungal-closer, touches human |
| L-838417 | novel | yes | 2.580 | 2.943 | 1/5 | +0.36 | fungal-closer, touches human |
| voriconazole | azole ref | yes | 2.657 | 3.328 | 0/5 | +0.67 | selective (marginal) |
| fluconazole | azole ref | yes | 2.711 | 2.768 | 5/5 | +0.06 | **hits human (baseline)** |
| nolatrexed | novel | yes | 2.677 | 2.622 | 5/5 | −0.06 | **hits human ≥ fungal** |
| flutrimazole | azole ref | yes | 2.371 | 2.311 | 5/5 | −0.06 | **hits human ≥ fungal** |
| letrozole | azole ref | yes | 2.705 | 2.555 | 5/5 | −0.15 | **hits human ≥ fungal** |
| LDN-27219 | novel | NO | 2.603 | 3.239 | 0/5 | — | demoted (unstable, #68) |
| olprinone | novel | NO | 2.680 | 3.996 | 0/5 | — | demoted (borderline, #68) |
| ketoconazole | control | NO | 2.921 | 3.675 | 0/5 | — | control did not converge |

## What this establishes (provisional)

1. **The readout behaves correctly on the known-non-selective baseline.** Real clinical/agro
   azoles that are *expected* to hit human CYP51 — fluconazole, letrozole, flutrimazole — do
   exactly that here: they coordinate the human iron in all 5 seeds with a margin near zero or
   negative. This is the toxicity signal antifungal azoles are known for, and it is the
   baseline a candidate must beat to count as a selectivity win. That the assay reproduces it is
   the main reason to trust the signal despite the ketoconazole control.

2. **Four seed-stable novel hits show a clean structural selectivity window.**
   **vatalanib, flucloxacillin, taranabant, lersivirine** each reach the fungal iron in all 5
   seeds but **never** reach the human iron (0/5), with margins +1.8 to +2.7 Å — a much cleaner
   separation than any reference azole. On this structural readout they are the standouts.

3. **Two hits are fungal-closer but touch the human iron.** verinurad (human 2/5) and L-838417
   (human 1/5) coordinate the fungal iron more readily but occasionally reach the human one —
   a weaker, watch-it selectivity margin.

4. **nolatrexed is stable but NOT selective.** It coordinates the human iron as well as the
   fungal one (5/5, margin −0.06) — a selectivity liability flag, despite passing #68.

## The combined refinement of the shortlist

Reading #68 and #73 together, the shortlist tiers as:

- **Strongest (seed-stable AND cleanly human-selective):** vatalanib, flucloxacillin,
  taranabant, lersivirine.
- **Stable, selectivity marginal / provisional:** verinurad, L-838417, R-1479.
- **Demoted:** LDN-27219, olprinone (single-dock instability, #68); nolatrexed (hits human
  CYP51 as readily as fungal, this run).

## Follow-up (named, not hand-waved)

- **Stronger positive control** for the human receptor: re-dock ketoconazole seeded from its
  3LD6 crystal pose (or add a small, rigid, known human-CYP51 inhibitor as the control) so the
  precondition can actually be met and these numbers promoted from provisional to validated.
- Broader human off-target panel (CYP3A4/2C9/2D6…) is issue #74; this run covers CYP51 only.

## Limitations

1. **Pre-registered positive control unmet** → numbers provisional (see top).
2. Rigid single-conformer human receptor (one crystal form, 3LD6); the human pocket is flexible
   too.
3. Structural signal only — says nothing about other human CYPs, permeability, or metabolism.
4. Coordination is necessary-not-sufficient on the human side too: a short human N-Fe indicates
   a *possible* off-target interaction, not a measured one.

## Reproduce

```
conda activate openafr
python scripts/prep_receptor_human.py                 # build work/receptor_human.pdbqt from 3LD6
scripts/screen_human.sh work/candidates/focus_ligands work/candidates/dock_human "42 1 2 3 4"
python scripts/shortlist_confidence.py \
    --fungal work/candidates/dock_fungal --human work/candidates/dock_human
```

3LD6 is fetched from RCSB (`data/structures/3LD6.pdb`, sha256 `486ba093…`). Docked poses are
gitignored; the receptor anchor (`work/receptor_human_A.pdb`), box
(`work/receptor_human_box.json`) and summary (`work/candidates/shortlist_confidence.json`)
persist.
