# Repurposing screen — first novel-library run

Date: 2026-08-05

The first time the pipeline has been pointed at a **novel library** rather than its own
validation set. This is the pipeline's intended use per [VISION.md](../VISION.md): point it
at a candidate library → dock once under the frozen protocol → get a short, evidence-backed
ranked shortlist. It is a demonstration of the *path*, not a drug.

> **Every rank here is a hypothesis for a wet lab, never a hit.** No molecule below has been
> tested against any fungus. The list is meaningful *only* because the validation gate passed
> first (below), and only within the declared applicability domain.

## Precondition: the gate was re-verified PASS in this environment

Screening a novel library is only justified once `validate_gate2.py` passes under the frozen
protocol. Before this run, the run-2 held-out validation was re-docked and re-graded from
scratch in this environment. It reproduced the published result **bit-for-bit**:

| mode | AUC | EF@1% | role |
|---|---|---|---|
| **C — iron-approach distance** | **0.794** | **12.68x** | PRIMARY — **PASS** |
| A — Vina score alone | 0.471 | 0.00x | secondary — worse than random |

Pre-registration (`sha256 7c89d78c…`) and protocol (`sha256 6eca16b1…`) both verified
unmodified at grading time. See [RESULTS_run2_holdout.md](RESULTS_run2_holdout.md) for the
original run.

## Library

- **Source:** the Broad Institute **Drug Repurposing Hub**, sample table dated 2020-03-24
  (`repurposing_samples_20200324.txt`, https://clue.io/repurposing). For non-commercial use —
  consistent with this repository's PolyForm Noncommercial license.
- **Provenance chain (nothing dropped silently):**

  | stage | count | tool |
  |---|---:|---|
  | unique-named compounds in the hub | 6798 | dedupe by `pert_iname` |
  | in applicability domain (≤ 45 heavy atoms **and** ≥ 1 aromatic N) | 2797 | `scripts/filter_library.py` |
  | 3D conformer embedded | 2787 | `scripts/prep_ligands.py` (10 embedding failures, logged) |
  | **docked and ranked** | **2776** | `scripts/screen.sh` protocol (11 docking failures, logged) |

  The 10 embedding failures are mostly cinchona alkaloids (quinine, quinidine and
  relatives) plus two metal complexes; all are named in `work/repurposing/prep.log`.

- **Protocol:** identical to the frozen gate — 26 Å box @ (70.61, 66.28, 4.18),
  exhaustiveness 32, 20 modes, seed 42, rigid 5TZ1 chain A + heme. No per-molecule tuning.
  For throughput the 2776 docks were run with Vina `--cpu 1` across parallel processes; this
  was verified to produce **byte-identical** poses to `scripts/screen.sh` (Vina is
  deterministic given the seed and exhaustiveness, independent of thread count), so the
  frozen protocol is unchanged.

## The run validated itself: known antifungals re-surfaced at the top

The library was screened blind — the pipeline was never told which of the 2776 drugs are
antifungals. The ranking criterion, on its own, pulled known azole antifungals to the top of
the list. This is an **embedded positive control** on a genuinely novel library:

| rank (of 2776) | azole antifungal | N–Fe (Å) | percentile |
|---:|---|---:|---|
| 1 | flutrimazole | 2.371 | top 0.04% |
| 9 | ravuconazole | 2.653 | top 0.3% |
| 22 | fosfluconazole (fluconazole prodrug) | 2.717 | top 0.8% |
| 31 | fluconazole | 2.740 | top 1.1% |
| 57 | voriconazole | 2.799 | top 2.1% |
| 73 | luliconazole | 2.817 | top 2.6% |
| 89 | croconazole | 2.838 | top 3.2% |

**Seven known azole antifungals land in the top 3.2%**; the criterion also ranks the
heme-coordinating CYP inhibitor **letrozole at #23**. This is consistent with the held-out
gate result and is the sanity check that the screen behaved. (One outlier, clotrimazole, is
discussed in the limitations.)

## Output: the top novel hypotheses

The known azoles above are controls, not discoveries. The *repurposing hypotheses* are the
high-ranking drugs approved or developed for **unrelated** indications — the molecules whose
geometry unexpectedly brings a nitrogen as close to the heme iron as a real azole does:

| rank | drug | usual target / class | N–Fe (Å) | Vina |
|---:|---|---|---:|---:|
| 2 | verinurad | URAT1 inhibitor (gout) | 2.417 | -7.3 |
| 3 | vatalanib | VEGFR tyrosine-kinase inhibitor | 2.479 | -11.0 |
| 4 | flucloxacillin | β-lactam antibiotic | 2.558 | -7.8 |
| 5 | taranabant | CB1 inverse agonist | 2.586 | -10.5 |
| 6 | lersivirine | HIV NNRTI | 2.599 | -7.1 |
| 7 | LDN-27219 | transglutaminase-2 inhibitor | 2.603 | -9.3 |
| 8 | R-1479 (balapiravir) | HCV polymerase inhibitor | 2.621 | -7.4 |
| 12 | epacadostat | IDO1 inhibitor (oncology) | 2.674 | -8.5 |
| 14 | nolatrexed | thymidylate synthase inhibitor | 2.692 | -9.4 |
| 19 | apatinib | VEGFR2 inhibitor | 2.711 | -10.7 |

Full detail — per-candidate iron-coordinating pose counts, the demoted docking score, and
provenance — is in [`work/repurposing/REPORT_repurposing.md`](repurposing/REPORT_repurposing.md).
The complete ranking of all 2776 is `work/repurposing/ranked.tsv`.

Several top hits are N-heterocyclic kinase inhibitors and azole-adjacent scaffolds — chemistry
that can coordinate a heme iron. That is exactly what the criterion measures, and exactly why
none of this is a hit yet: coordinating the iron in a static pose is necessary, not sufficient.

## Honest limitations

1. **A rank is a hypothesis, not a hit.** Nothing here was tested against any fungus. The
   whole point of the project is that this list makes a wet-lab yes/no *more likely than
   random*, not that it is a yes.
2. **The controls also expose the noise floor.** Clotrimazole — a known azole antifungal —
   ranks **#2774 of 2776** (best pose 17.97 Å from the iron): docking simply failed to find a
   coordinating pose for it in 20 modes. A criterion that ranks a real drug dead last on one
   molecule is a criterion with real per-molecule variance; the top ranks are enriched, not
   infallible.
3. **SMILES taken as the hub provides them.** Salts, stereochemistry and protonation are as
   supplied; a single conformer was docked per molecule. No manual curation.
4. **Same receptor limitation as the gate.** One rigid *C. albicans* 5TZ1 structure, not
   *C. auris*, no induced fit. The applicability domain (≤ 45 heavy atoms) was enforced, but
   the receptor is still a stand-in for the resistant target that matters.
5. **Decoys were presumed inactive in the gate; drugs are presumed nothing here.** This screen
   ranks; it does not label. Iron coordination is one mechanistic filter, not a binding assay.

## Reproduce

```
conda env create -f environment.yml && conda activate openafr
python scripts/prep_receptor.py                          # -> work/receptor.pdbqt (hash-checked)
# gate must pass first:
scripts/run_screen2.sh && python scripts/validate_gate2.py work/screen2
# then the novel screen:
curl -sL https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_samples_20200324.txt -o samples.txt
#   -> build SMILES<TAB>name, then:
python scripts/filter_library.py work/repurposing/library_all.smi -o work/repurposing/library_in_domain.smi
python scripts/prep_ligands.py work/repurposing/library_in_domain.smi work/repurposing_ligands
scripts/screen.sh work/repurposing_ligands work/repurposing_screen
python scripts/rank_candidates.py  work/repurposing_screen -o work/repurposing/ranked.tsv
python scripts/report_candidates.py work/repurposing_screen -n 25 -o work/repurposing/REPORT_repurposing.md
```

The filtered library inputs (`library_all.smi`, `library_in_domain.smi`,
`library_excluded.tsv`), the full ranking (`ranked.tsv`) and the report are versioned under
`work/repurposing/`. The pose files and prepared ligands (~180 MB) are regenerated, not
versioned, per the repository's docking-intermediates convention.
