# ADMET / structural-alert profile of the candidate shortlist (#75)

Date: 2026-08-26

The chemistry-realism axis the shortlist was missing. The geometry criterion ranks a
molecule by how close its nitrogen comes to the heme iron — it says nothing about whether
that molecule could ever *be* a drug: soluble, permeable, non-reactive, free of a known
toxicophore. This stage adds that axis (`openafr/admet.py`, `scripts/admet_filter.py`),
a fully **local, deterministic RDKit computation** — no network — over every candidate's
SMILES.

> **Informational, never a gate.** Ro5/Veber are 30-year-old heuristics with famous
> exceptions; PAINS/BRENK/NIH flag *possible* liabilities that context can excuse; the
> hERG flag is a pharmacophore hint, not an IC50. Nothing here rejects a candidate. It
> tells a reader what to go check before spending bench time. Every marketed azole below
> bends at least one rule and is still a drug.

## The reference drugs validate the profiler in-run

The shortlist carries the real azoles as controls, and they profile exactly as a working
profiler should:

- **fluconazole, voriconazole — clean on every axis** (0 Ro5 violations, Veber OK, no
  PAINS/alerts, no hERG flag, log S −2.4 / −3.6). These are the best-tolerated, most
  soluble azoles in the clinic; a profiler that flagged them would be crying wolf.
- **ketoconazole — flags its own known liabilities.** MW 531 (one Ro5 violation), poor
  solubility (log S −5.7), a PAINS aniline substructure, **and the hERG risk-factor flag**.
  Ketoconazole is the azole withdrawn/black-boxed for hepatotoxicity and QT prolongation —
  the profiler independently recovers the shape of its real cardiac liability.

That the same rules pass the tolerated azoles and flag the dangerous one is the live proof
the layer discriminates — not just that it ran.

## What it caught among the novel hypotheses

| candidate | class | flags | note |
|---|---|---|---|
| **R-1479** | HCV nucleoside | **PAINS + azide + diazo + quaternary-N** | an azido-ribonucleoside — a stack of reactive/unstable alerts; chemically implausible as a CYP51 lead, read the alerts before anything else |
| **taranabant** | CB1 inverse agonist | **2× Ro5, insoluble** (log S −6.9, cLogP 6.3) | too greasy/insoluble to be a comfortable oral lead as-is |
| **LDN-27219** | TG2 inhibitor | acyl-hydrazine, hydrazine, N–O bond (BRENK) | reactive hydrazide warhead — a metabolic/tox liability |
| **flucloxacillin** | β-lactam antibiotic | β-lactam alert (NIH) | expected (it *is* a β-lactam); the reactive ring is the drug's mechanism, not a surprise |
| **flutrimazole** | topical azole | triphenyl-methyl (BRENK) | the trityl-like azole core; a topical drug, flagged for the bulky lipophilic centre |
| **vatalanib** | VEGFR TKI | hERG risk factors | basic center + lipophilic + aromatic bulk — go read its cardiac profile |
| verinurad, lersivirine, ravuconazole, L-838417, olprinone, nolatrexed, letrozole | — | **clean** | drug-like on every computed axis |

## How to read the flags

- **`clean`** — drug-like on every computed axis. Not a guarantee (these are estimates),
  but no computed red flag.
- **`Ro5xN` / `Veber` / `insoluble`** — physicochemical drug-likeness concerns. Soft;
  weigh against the fact that the reference azoles also bend rules.
- **`PAINS`** — a pan-assay-interference substructure. Says "this may read the assay, not
  the target" — worth extra scrutiny of any activity, not an automatic disqualification.
- **`alert`** — a BRENK/NIH reactive-or-undesirable fragment. Read the specific fragment;
  some (β-lactam here) are the molecule's intended mechanism, others (hydrazine, azide)
  are genuine liabilities.
- **`hERG?`** — the classic hERG pharmacophore *risk factors* only (a protonatable amine +
  lipophilicity + aromatic bulk). **Not a prediction** — a hint to go read the cardiac
  liability.

## What this is not

Every number is a computed estimate, not a measurement. For a candidate that is already an
approved drug, the authoritative ADMET/tox source is the **drug label**, not these
estimates — several shortlist entries are repurposed drugs whose real liabilities are
already documented (this is the useful signal #75 asked to fold in). The hERG flag ships no
model; a true off-target CYP panel (#74) and species/mutant retention (#69, #77) are
separate axes. This layer is one column in the eventual per-candidate scorecard (#88), not
a verdict on its own.

## Reproduce

```
conda activate openafr
python scripts/admet_filter.py work/candidates/shortlist_focus.smi          # table
python scripts/admet_filter.py work/candidates/shortlist_focus.smi --tsv    # flat table
python scripts/admet_filter.py work/candidates/shortlist_focus.smi --json   # full profiles
python -m pytest tests/test_admet.py -q                                     # 7 tests
```

Artifacts: `work/candidates/shortlist_admet.tsv`, `work/candidates/shortlist_admet.json`.
Module `openafr/admet.py`; CLI `scripts/admet_filter.py`; tests `tests/test_admet.py`.
