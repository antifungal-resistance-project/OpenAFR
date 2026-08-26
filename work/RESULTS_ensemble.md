# Ensemble receptor — the within-azole question, re-opened (look #8)

> **SUPERSEDED IN INTERPRETATION 2026-08-26 by look #9 (`work/RESULTS_ensemble_confirm.md`).**
> This result stands exactly as graded (AUC 0.750 on the 7 held-out azoles), but it was reported
> *with* an explicit n=7 wide-CI caveat and an explicit demand for a pre-registered confirmation
> on an independent active set. That confirmation ran on an independent 200-azole sample and
> **FAILED at AUC 0.682 < 0.70** — the 0.750 did not replicate and was consistent with the
> 7-molecule sampling fluctuation this document's own CI (0.683–0.825) flagged. The within-azole
> line is now CLOSED on both the rigid and the ensemble receptor. Read look #9 for the settled
> conclusion; the value of this look was to make that test well-posed, not to establish a result.

Date: 2026-08-25
Pre-registration: `work/PREREGISTRATION_ensemble.md` (sha256 `8041b38ef37dde7a…`,
verified unmodified by the grader at run time)
Ensemble: the **exhaustive** set of experimental *C. albicans* CYP51 crystal structures —
5TZ1 (2.0 Å, VT-1161), 5FSA (2.86 Å, posaconazole), 5V5Z (2.9 Å, itraconazole). Receptor
anchors hash-pinned; all three verified at run time. 5FSA/5V5Z Cα-superposed onto the 5TZ1
frame (RMSD 1.82 / 1.06 Å); the three heme irons land within 0.3 Å, so the frozen box seats
identically in each pocket.
Criterion: mode C (min ligand-N-to-heme-Fe distance), **unchanged**, aggregated by the
pre-committed **ensemble-minimum** over the 3 conformers. Protocol unchanged, hash-verified —
box 26 Å @ 70.61/66.28/4.18, exhaustiveness 32, num_modes 20, seed 42.
Docking: 292/294 molecules in **every** conformer (the same 2 inactives —
`inactive_CHEMBL4551289_for_fenticonazole`, `inactive_CHEMBL469908_for_efinaconazole` —
failed obabel/Vina conversion identically in all three, so the pool stays matched; both are
inactive, no active was lost). ~99 min wall × 3 conformers on 12 cores, Apple Silicon.

## RESULT: GATE PASSED — the within-azole ceiling was a single-conformer artefact

**H1 passes: ensemble-min mode C AUC 0.750 ≥ 0.70**, on the *identical* 169-molecule pool
(7 held-out azoles + 162 azole-bearing measured inactives) where the single rigid 5TZ1
conformer scored 0.650 (mode C, look #6) and 0.590 (mode F, look #7). A real experimental
conformer ensemble ranks working azoles above failed azole analogues where one rigid pose
could not. **The two prior FAILs were properties of the rigid single conformer, not of the
geometric criterion** — exactly the leading explanation both pre-registrations committed to in
advance.

| ranking | pool | AUC | vs bar |
|---|---|---|---|
| **ensemble-min mode C (5TZ1+5FSA+5V5Z)** | 7 actives + 162 azole-bearing inactives | **0.750** | **PASS** |
| ensemble-mean mode C (diagnostic) | same | 0.714 | (does not gate) |
| mode C — single conformer 5TZ1 (look #6) | same | 0.650 | FAIL |
| mode F — channel engagement (look #7) | same | 0.590 | FAIL |

The novel-chemotype **guard** — the pre-committed check that the ensemble does not *break* the
thing that already works — not only holds but **improves**: AUC **0.949** vs non-azole
inactives (comparator 0.810, single conformer). So the ensemble does not merely make every
molecule coordinate more easily; the held-out actives benefit more than the non-azole decoys.
This is the clean PASS branch of the pre-committed interpretation, not the "guard collapsed" one.

## Where the lift comes from — one conformer, honestly

The "ensemble" benefit is not three conformers contributing democratically. It is almost
entirely the **posaconazole-bound 5FSA** conformation:

    per-conformer AUC on the gate pool:   5TZ1  0.650   (reproduces look #6 exactly — sanity check PASSES)
                                          5FSA  0.783
                                          5V5Z  0.727
    marginal-conformer increment (min):   {5TZ1}            0.650
                                          {5TZ1, 5FSA}      0.751   <- 5FSA does the work
                                          {5TZ1,5FSA,5V5Z}  0.750   <- 5V5Z is ~neutral

The 5TZ1-alone column landing on 0.650 is the built-in pipeline sanity check (it must
reproduce look #6, and it does — the ensemble code is not silently re-scoring). Mechanistically
the result is sensible: posaconazole is the longest-tailed azole, so its co-crystal channel is
the most *open*, and it is that open channel that accommodates the held-out azoles' tails. The
honest headline is therefore narrower than "an ensemble lifts the ceiling": **a single
wide-channel experimental conformation lifts it, and the exhaustive ensemble captures that
conformation.**

## The caveats that keep this from being "within-class triage solved"

**Early enrichment is poor.** The pass is a *broad-ranking* AUC pass, not top-of-list
concentration:

    EF@1%  0.00x     EF@5%  0.00x     EF@10%  1.42x     BEDROC(20)  0.054
    best active: fenticonazole, rank 14/169 (8.28th percentile), 13 inactives above it

So the ensemble separates actives from failed analogues *on average*, but does **not** pack them
at the very top of the list — the first true active is only at rank 14. A wet-lab shortlist
built from the top few would still be dominated by failed analogues. AUC moved; the tip did not.

**n = 7, wide CI.** Bootstrap 95% CI over the actives is **0.683 – 0.825** — the lower bound sits
*just below* the 0.70 bar, exactly the small-count fragility pre-registered as Limitation 1. The
permutation test is significant (p = 0.0118, 20,000 shuffles), so the ranking is not noise, but
the point estimate is not comfortably clear of the bar.

**n = 15 powered form: AUC 0.760** (self-pair excluded and raw are both 0.760 — posaconazole and
itraconazole draw their ensemble-min from a non-self conformer anyway, so there is no detectable
self-docking inflation). Consistent with the gate, not gating.

**Holo bias (Limitation 2).** All three conformers are azole-bound, so the channel is already in
azole-accommodating states — a thumb on the scale *toward* the actives. A PASS must be read
against that; it is why the honest next step is confirmation, not a wet-lab claim.

## Pre-committed interpretation (from the pre-registration, verbatim branch)

> **PASS — AUC ≥ 0.70, AND the novel-chemotype guard still ≥ ~0.78.** A real experimental
> conformer ensemble lifts within-azole ranking above the single-rigid-conformer ceiling
> *without* breaking novel-chemotype triage. This is the concrete evidence that the mode C and
> mode F FAILs were **induced-fit / single-conformer artefacts**, not failures of the geometric
> criterion. The within-class triage product re-opens, still at n=7 with a wide CI; **a
> separately pre-registered confirmation on an *independent* active set is required before any
> within-class wet-lab claim.**

## What this changes, and what it does not

- **Does:** overturn the look #6/#7 conclusion that the within-azole ceiling was intrinsic. It
  was not — it was the single rigid conformer. Adding the exhaustive experimental ensemble lifts
  mode C from 0.650 to 0.750 and *raises* novel-chemotype enrichment to 0.949. Every candidate on
  the repurposing shortlist (`RESULTS_candidate_shortlist.md`) can now be re-scored under the
  ensemble, hardening it against the "single rigid conformer" caveat that qualified all of them.
- **Does not:** claim within-class azole triage is solved or ready for a lab. Early enrichment is
  poor (first active at rank 14), the CI straddles the bar, the lift rides on one holo conformer,
  and n = 7. Per the pre-registration, the required next step is a **separately pre-registered
  confirmation on an independent active set** — not a wet-lab shortlist claimed within-class, and
  not another feature on these poses.

## Reproduce

```
conda activate openafr
python scripts/prep_receptor.py                          # 5TZ1 -> work/receptor.pdbqt
python scripts/prep_receptor_ensemble.py 5FSA 5V5Z       # align + prep the 2 new conformers
python scripts/prep_ligands.py data/ligands/actives_holdout_final.smi work/verified_ligands
python scripts/prep_ligands.py data/ligands/verified_inactives.smi    work/verified_ligands
python scripts/prep_ligands.py data/ligands/actives.smi               work/verified_ligands
RECEPTOR=work/receptor.pdbqt      scripts/screen.sh work/verified_ligands work/screen_verified
RECEPTOR=work/receptor_5FSA.pdbqt scripts/screen.sh work/verified_ligands work/screen_5FSA
RECEPTOR=work/receptor_5V5Z.pdbqt scripts/screen.sh work/verified_ligands work/screen_5V5Z
python scripts/validate_gate_ensemble.py                 # exits 0 (PASS)
```
</content>
