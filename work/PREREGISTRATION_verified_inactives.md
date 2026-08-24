# Pre-registration — verified inactives: an independent decoy construction

**Written 2026-08-23, BEFORE a single molecule of this set was docked.** The inactive set was
built first (it has to exist to be described), but nothing has been docked against it, no pose
has been read, and no metric has been computed. Everything below is fixed in advance.

## Background — why this look is permitted at all

The published result (`work/PREPRINT_geometry_ceiling.md`) closed its own dataset. §5:

> **The stopping rule is now in force: no further criterion will be tried on these poses.**
> Any further attempt requires an independent dataset — new actives or a new decoy
> construction — not another feature on this one.

Five looks are on the ledger (run-2 baseline, Y132F transfer, consensus, reliability, axial),
all against the same 7 held-out azoles and the same 348 **presumed**-inactive decoys. This is
look **#6**, and the first that changes the *data* rather than the criterion. It is a new
decoy construction, which §5 explicitly names as legitimate.

It also addresses the paper's own stated priority. §4.1, on the decoy-construction bind:

> The only clean exit is **experimentally confirmed inactives** — compounds measured not to
> inhibit — rather than computationally presumed ones. We do not have those; obtaining them
> for CYP51 is, in our view, the most useful next contribution anyone could make to this
> benchmark.

And Limitation 2: *"Decoys are presumed inactive, not verified. … This is the single largest
caveat."* This run removes that caveat, at the price of a different one (stated in full below).

## The dataset (already built and frozen — recorded here so it cannot be silently re-filtered)

Built by `scripts/make_verified_inactives.py` from a ChEMBL bioactivity pull; the evidence rule
and every filter live in `openafr/inactives.py` and are locked by `tests/test_inactives.py`.

**Evidence rule.** Per activity record against *C. albicans*:

| verdict | rule |
|---|---|
| inactive | MIC `> X` µg/mL with X ≥ 32 (tested to a real ceiling, nothing seen), or MIC `=` X ≥ 64 µg/mL, or an explicit "Not Active"/"Inactive" comment |
| active | MIC ≤ 8 µg/mL (or ≤ 10 µM), an "Active" comment, or pChEMBL ≥ 5 |
| ambiguous | measured, but in the 8–64 µg/mL band, or censored below 32 |
| ignore | no usable quantitative or textual verdict |

A compound qualifies as a **verified inactive** only under *consistent inactivity*: ≥ 1
`inactive` record and **zero** `active` and **zero** `ambiguous` records, across every assay it
appears in — including the *C. albicans* CYP51 enzyme target (CHEMBL1780), so a compound with
measured enzyme activity that merely fails a whole-cell assay is excluded, not counted inactive.

**Provenance and counts, fixed:**

| item | value |
|---|---|
| ChEMBL pull, whole-cell *C. albicans* | 87,879 records, sha256 `c97f0dca3191ce3b…` |
| ChEMBL pull, *C. albicans* CYP51 (CHEMBL1780) | 816 records, sha256 `fa2b027f6cf23e40…` |
| pulled (UTC) | 2026-08-24T03:34:00Z |
| distinct compounds ever tested | 43,004 |
| verdicts | 9,270 active · 8,908 **verified-inactive** · 5,306 ambiguous · 19,520 no-usable-evidence |
| verified inactives with a structure | 8,871 |
| after filters (funnel below) | **279** |
| `data/ligands/verified_inactives.smi` | sha256 `37a20eec1c026c42c…` (full value in the results) |

**Filter funnel** (8,871 in → 279 kept): 3,569 rejected for no aromatic nitrogen, 4,011 for no
property match, 378 for a full per-active quota, 348 for exceeding the ≤ 45 heavy-atom
applicability domain, 286 for Tanimoto ≥ 0.35 against an active, 0 unparseable, 0 duplicate.

**The single-variable swap.** Every filter that made the presumed-inactive set honest is held
byte-identical to `scripts/make_decoys.py`: property matching (MW ± 25, cLogP ± 1.5, rotatable
bonds ± 2, HBD ± 1, HBA ± 2), the ≥ 1 aromatic-nitrogen mechanism requirement, Tanimoto < 0.35
novelty against *every* active (training and held-out), ≤ 50 per active. The only thing that
changes is where "inactive" comes from. Achieved match, verified inactives vs the 7 actives:
MW 379.9 vs 392.5, cLogP 4.8 vs 5.7, **heavy atoms 26.3 vs 25.4**, rotatable bonds 4.6 vs 5.4.
The heavy-atom agreement is the load-bearing one — docking score correlates with heavy-atom
count at r = −0.91 on this system, so a size-mismatched set would "enrich" for the wrong reason.

## The one thing that is NOT like-for-like, stated before the result

You can only measure what somebody tested, and the compounds tested against *C. albicans* are
overwhelmingly antifungal medicinal-chemistry programs. So this set is chemically different from
a set drawn from all of ChEMBL, and the difference is measurable — computed before docking:

| set | fraction carrying an azole (imidazole/triazole/tetrazole) ring |
|---|---|
| held-out actives (n=7) | 100 % |
| presumed-inactive decoys (n=350) | 27.4 % |
| **verified inactives (n=279)** | **58.4 %** |

**More than half of the verified inactives are failed azole analogues** — molecules carrying the
same iron-coordinating warhead as the drugs, from real optimization campaigns, measured not to
work. This makes the benchmark **substantially harder** than the presumed one, and it changes
what a FAIL would mean (see the interpretation branches). It is a confound in any direct
comparison to the run-2 numbers, and it is inherent to the exercise: a verified inactive can
only come from a compound someone bothered to test.

## What is held fixed (nothing about the method changes)

| | |
|---|---|
| criterion | **mode C** — minimum ligand-nitrogen-to-heme-iron distance, unchanged since run 2 |
| protocol | `protocol.yaml` sha256 `6eca16b1be2e5e33…` — box 26 Å @ (70.61, 66.28, 4.18), exhaustiveness 32, 20 modes, **seed 42** |
| receptor | `work/receptor_A.pdb` sha256 `a1a8410868d1b9ee…` → `work/receptor.pdbqt` sha256 `7ed596cb906f82ae…` |
| actives | `data/ligands/actives_holdout_final.smi` (n=7) sha256 `0d964c885d51f6c6…`, unchanged since run 2 |
| ligand prep | `scripts/prep_ligands.py`, ETKDGv3 seed 42, MMFF optimization |
| docking | `scripts/screen.sh`, one Vina process per ligand, frozen seed |

Single seed 42, matching the **run-2 baseline** exactly, which makes run 2 the direct
comparator: same criterion, same protocol, same receptor, same actives, same seed — different
decoy provenance. (The run-2 baseline is AUC 0.794 with EF@1% 12.68x on 348 presumed decoys;
those numbers stand as graded and are not recomputed here.)

## Primary hypothesis (H1) — the only thing that gates

Ranking the pooled 286 molecules (7 actives + 279 verified inactives) by mode C separates
actives from measured non-inhibitors with **AUC ≥ 0.70** — the `protocol.yaml` `min_auc`,
unchanged, not lowered.

**AUC alone gates. EF does not, and here is the arithmetic reason, computed in advance** (the
preprint's own practical recommendation 2: *"compute the achievable values first; if the metric
is Bernoulli, it is not a gate"*). With 286 molecules and 7 actives:

| metric | top slice | achievable values (0…7 actives in the slice) |
|---|---|---|
| EF@1% | **2 molecules** | 0, 13.62, 27.24, 40.86, then saturated — effectively Bernoulli |
| EF@5% | 14 molecules | 0, 2.72, 5.45, 8.17, 10.90, 13.62, 16.34, 19.07 |
| EF@10% | 28 molecules | 0, 1.41, 2.82, 4.23, 5.64, 7.04, 8.45, 9.86 |

EF@1% cannot be a gate on a two-molecule slice, and is reported here only because run 2 reported
it. **No EF value, at any fraction, will be used to declare a pass or a failure.**

## Secondary, all pre-specified, none gating

1. **H2 — the tip.** On presumed decoys, six outranked the best true active. The pool-size-neutral
   statistic is the **percentile rank of the best-ranked true active**; the raw count of
   inactives above it is reported alongside (a raw count is not comparable across pool sizes).
2. **EF@1/5/10%, BEDROC(α=20)**, and **mode A** (rank by Vina score), as in run 2.
3. **Two pre-specified subgroups**, because the composition finding above makes them the
   interesting question — and pre-specifying them is what stops them from being post-hoc slices:
   - **azole-bearing inactives (n=163)** — "can geometry tell a working azole from a failed
     azole analogue?" This is the medicinal-chemistry question and the hardest form of the test.
   - **non-azole inactives (n=116)** — the comparison group.
   Both are reported with AUC. Neither gates, and no subgroup result will be promoted to the
   headline if the primary fails.

## Analysis plan

    Primary (gates):  mode C, pooled 286 molecules, AUC >= 0.70
    Reported:         EF@1/5/10%, BEDROC(alpha=20), mode A, H2 tip percentile,
                      the two pre-specified subgroups
    Unrankable actives (no nitrogen-bearing pose in any mode) are ranked LAST, never dropped.
    Metrics from rdkit.ML.Scoring, via openafr/scoring.py — the same code that graded run 2.

Graded by `scripts/validate_gate_verified.py`, which re-checks this document's SHA-256 and the
frozen protocol hash before reporting and refuses to certify a pass if either moved.

## Pre-committed interpretation

- **PASS (AUC ≥ 0.70).** The separation published against presumed-inactive decoys is **not an
  artifact of decoy construction** — the criterion also separates real azoles from compounds
  *measured* not to inhibit, on a set where 58% carry the azole warhead. Limitation 2 of the
  preprint is discharged, and the top-5% shortlist claim is supported on an independent dataset.
  This does **not** license any resistance-breaker claim (§4.4 stands), and does not reopen the
  top-1% claim, which the two-molecule slice above cannot address either way.
- **FAIL (AUC < 0.70).** The criterion cannot separate real azoles from measured non-inhibitors.
  This is the **stronger negative**: the ceiling is not a property of how the decoys were built,
  and the shortlist claim must be re-scoped or withdrawn. The composition caveat is then
  load-bearing and must be reported with it — the honest statement would be *"geometry does not
  separate working azoles from failed azole analogues"*, which is narrower than *"geometry does
  not separate actives from inactives"*, and the pre-specified non-azole subgroup is what
  distinguishes those two readings.

**One look. No re-filtering.** If the result disappoints, the inactive set will **not** be
rebuilt with different thresholds, a different MIC cutoff, a stricter consistency rule, or a
different assay selection, and re-run. The evidence rule was fixed in `openafr/inactives.py`
before this document and is locked by tests. Any future change to it produces a *separately
pre-registered* dataset, disclosed as such, with this result left standing as graded.

## Limitations (fixed in advance)

1. **A whole-cell MIC is not a CYP51 assay.** A compound can inhibit the enzyme yet fail the
   cell assay through permeability or efflux, which would make it a false inactive here. Where
   ChEMBL has direct CYP51 data we exclude such compounds, but coverage is thin (816 records).
   The residual error is **conservative** — a false inactive depresses measured separation, it
   cannot inflate it.
2. **Composition is not matched to the presumed set** (58.4% vs 27.4% azole-bearing), for the
   structural reason given above. The comparison to run 2 is therefore *harder-benchmark*, not
   like-for-like.
3. **ChEMBL inactivity is heterogeneous**: different labs, media, inoculum, endpoints, decades.
   Consistency across records is required, but the records are not a single controlled screen.
4. **n = 7 actives**, unchanged from run 2 — every small-count caveat from the preprint's
   Limitation 1 carries over, including a bootstrap CI on AUC that may straddle the bar.
5. **One docking program, one conformer, one rigid receptor** (5TZ1, ≤ 45 heavy-atom domain),
   single seed. Nothing here revisits the induced-fit or resistance-breaker questions.
6. **Publication bias in what gets tested.** Compounds tested against *C. albicans* and found
   inactive are more likely to be published as part of a series than as singletons; the set
   inherits whatever that selection does.

## Reproduce

```
conda activate openafr
python scripts/make_verified_inactives.py fetch          # ~88k ChEMBL records -> work/chembl/
python scripts/make_verified_inactives.py build          # -> data/ligands/verified_inactives.smi
python scripts/prep_receptor.py
python scripts/prep_ligands.py data/ligands/actives_holdout_final.smi work/verified_ligands
python scripts/prep_ligands.py data/ligands/verified_inactives.smi   work/verified_ligands
scripts/screen.sh work/verified_ligands work/screen_verified
python scripts/validate_gate_verified.py work/screen_verified        # exits 0 pass / 1 fail
```
