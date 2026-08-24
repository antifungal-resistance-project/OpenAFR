# Verified inactives — the criterion tested against compounds MEASURED not to inhibit

Date: 2026-08-23
Pre-registration: `work/PREREGISTRATION_verified_inactives.md` (sha256 `d2926f7761496906…`,
verified unmodified by the grading script at run time)
Test set: 7 held-out azoles + **279 verified inactives** (284 ranked; 2 inactives failed to
dock, both reported, neither counted as a decoy)
Protocol: unchanged and hash-verified — box 26 Å @ 70.61/66.28/4.18, exhaustiveness 32,
num_modes 20, seed 42, rigid 5TZ1 chain A + heme. Criterion unchanged (mode C).
Host: Apple Silicon, `openafr` env, ~2.5 h wall for 286 dockings on 12 cores.

> **Erratum (not corrected in place, on purpose).** The pre-registration quotes the inactive
> set's sha256 inline as `37a20eec1c026c42c…`; the true value is
> `37a20eec1c026c42d3c8efb53c1dd27cda2f9057ca541cfc1afed21a35faea6b` (one transcribed character
> differs). The frozen document is **not** edited to fix it — a pre-registration amended after
> the result is known is worth nothing, and its own sha256 is what the grading script verifies.
> The authoritative, machine-checkable hashes are in `work/PREREG_verified_inactives.sha256`,
> which verifies clean (`shasum -a 256 -c`) for all three frozen files. Nothing about the run
> depends on the typo.

## RESULT: GATE PASSED — and the pre-specified subgroup is the real finding

**H1 passes: AUC 0.716 ≥ 0.70.** The geometric criterion does separate real azoles from
compounds *measured* not to inhibit *C. albicans*. **The published separation is not an
artifact of decoy construction** — which discharges Limitation 2 of the preprint, the one it
called "the single largest caveat."

But the pass is carried almost entirely by one half of the benchmark, and that was
pre-specified, not discovered afterwards:

| ranking | pool | AUC | EF@5% | BEDROC | actives in top 15 |
|---|---|---|---|---|---|
| **mode C — iron-approach distance** | 7 actives + 277 inactives | **0.716** | 2.70x | 0.157 | 1/7 |
| ├─ vs **non-azole** inactives | 7 + 115 | **0.810** | 4.98x | 0.487 | 4/7 |
| └─ vs **azole-bearing** inactives | 7 + 162 | **0.650** | 2.68x | 0.108 | 1/7 |
| mode A — Vina score alone | 7 + 277 | 0.436 | 0.00x | 0.002 | 0/7 |

**Against molecules that share the azole warhead but were measured not to work, the criterion
scores AUC 0.650 — below its own pass bar.** Against non-azole chemotypes it scores 0.810.
Most of what the criterion measures is *"is this an azole-like molecule?"*, not *"does this
azole work?"*

EF@1% is 0.00x and is not interpretable: the top 1% of 284 molecules is a **two-molecule
slice**, whose achievable values (0, 13.62, 27.24, …) were computed and published in the
pre-registration precisely so this number could not be read either way.

## Statistical support (descriptive; the gate is the AUC point estimate)

    permutation test (20,000 label shuffles): p = 0.0247      (run 2: 0.0028)
    bootstrap 95% CI over actives           : AUC 0.506 - 0.885   (run 2: 0.590 - 0.946)
    active rank positions (of 284)          : 8, 23, 43, 58, 72, 134, 240

As in run 2, the CI lower bound sits **below** the 0.70 bar. The point estimate passes and the
permutation test still rules out chance, but a true AUC under the bar cannot be excluded at
n = 7. This run is a *weaker* pass than run 2 on every one of these numbers.

## H2 — the tip, now measured rather than presumed

    best true active : fenticonazole, rank 8 of 284  (2.82nd percentile)
    inactives above it: 7
    run-2 comparator : best active rank 3 of 355 (0.85th percentile), 2 presumed decoys above

    rank  molecule                                  N-Fe
       1. inactive_CHEMBL476215_for_oxiconazole     2.378 A
       2. inactive_CHEMBL10290_for_efinaconazole    2.617 A
       3. inactive_CHEMBL275780_for_efinaconazole   2.629 A
       4. inactive_CHEMBL30105_for_butoconazole     2.644 A
       5. inactive_CHEMBL10471_for_efinaconazole    2.646 A
       6. inactive_CHEMBL305539_for_efinaconazole   2.665 A
       7. inactive_CHEMBL10631_for_efinaconazole    2.689 A
       8. fenticonazole                             2.703 A   <- first real drug
       9. inactive_CHEMBL273435_for_efinaconazole   2.707 A
      10. inactive_CHEMBL276017_for_efinaconazole   2.709 A

Composition of the tip: the top 10 is **1 active and 9 azole-bearing measured inactives**; the
top 20 is 1 active, 17 azole-bearing inactives, 2 non-azole. These are not look-alike molecules
assumed to be inactive. They are compounds someone synthesised, tested against *C. albicans*,
and measured **not** to work — and they coordinate the heme iron *more closely than six of the
seven real drugs do*.

**The single most useful number in this run:** of the 39 molecules landing inside the validated
iron-bound band (2.47–2.88 Å, the band `score_molecule.py` reports against), **33 are
azole-bearing compounds measured not to inhibit, and only 2 are true actives.**

## What this establishes, and what it does not

**Does establish:**
1. The decoy ceiling is **real, not an artifact of decoy construction**. §4.1 of the preprint
   worried that mechanism-matched decoys were built to be indistinguishable. They were not
   built at all here — they were measured, and they are still indistinguishable at the tip.
   The bind is a fact about the chemistry, not about our decoy generator.
2. Broad enrichment survives an independent dataset (AUC 0.716, p = 0.0247), so the top-5%
   shortlist claim stands on data the criterion has never seen, from a source that had no
   contact with this project.
3. The mechanism the enrichment actually exploits is **chemotype**, not potency: 0.810 against
   non-azoles versus 0.650 against failed azoles.

**Does NOT establish:**
1. Any ability to rank *within* the azole class — the practical medicinal-chemistry question.
   AUC 0.650 there is below the bar, and the honest reading of a shortlist produced by this
   tool is that it will be crowded with azole-like molecules that have already been tried.
2. Anything about resistance-breaking (§4.4 stands untouched).
3. Anything at the razor top: a two-molecule 1% slice cannot support an EF@1% claim in either
   direction.

## Head-to-head with run 2 (the only differences are the decoys)

| | run 2 (presumed decoys) | this run (measured inactives) |
|---|---|---|
| non-actives | 348, presumed inactive, 27.4 % azole-bearing | 277, **measured** inactive, 58.4 % azole-bearing |
| AUC (mode C) | 0.794 | 0.716 |
| permutation p | 0.0028 | 0.0247 |
| bootstrap 95 % CI | 0.590 – 0.946 | 0.506 – 0.885 |
| best active percentile | 0.85th (rank 3/355) | 2.82nd (rank 8/284) |
| verdict | PASS | PASS |

The degradation is expected and was predicted in the pre-registration: this benchmark is
**harder by construction**, because the only compounds anyone measures against *C. albicans*
are antifungal programs, so more than half of them carry the warhead. The comparison is
harder-benchmark, not like-for-like, and the drop should not be read as instability in the
method.

## Limitations (carried from the pre-registration, unchanged)

1. A whole-cell MIC is not a CYP51 assay; a compound can inhibit the enzyme and still fail the
   cell assay through permeability or efflux. Compounds with measured CYP51 enzyme activity are
   excluded where ChEMBL has the data (816 records), but coverage is thin. The residual error is
   conservative — a false inactive depresses measured separation, it cannot inflate it.
2. Composition is not matched to the presumed set (58.4 % vs 27.4 % azole-bearing).
3. ChEMBL inactivity is heterogeneous across labs, media, inoculum, endpoints and decades.
4. n = 7 actives — every small-count caveat carries over, including the CI above.
5. One docking program, one conformer, one rigid receptor, single seed.
6. Publication bias in what gets tested and reported inactive.

**One look, no re-filtering.** Per the pre-registration, the inactive set will not be rebuilt
with different thresholds and re-run. This result stands as graded.

## Reproduce

```
conda activate openafr
python scripts/make_verified_inactives.py fetch     # ~88k ChEMBL records -> work/chembl/
python scripts/make_verified_inactives.py build     # -> data/ligands/verified_inactives.smi
python scripts/prep_receptor.py                     # work/receptor.pdbqt sha256 7ed596cb906f…
python scripts/prep_ligands.py data/ligands/actives_holdout_final.smi work/verified_ligands
python scripts/prep_ligands.py data/ligands/verified_inactives.smi   work/verified_ligands
scripts/screen.sh work/verified_ligands work/screen_verified
python scripts/validate_gate_verified.py work/screen_verified
```

Poses (`work/screen_verified/`) are gitignored; the two docking failures
(`inactive_CHEMBL4551289_for_fenticonazole`, `inactive_CHEMBL469908_for_efinaconazole`)
reproduce.
