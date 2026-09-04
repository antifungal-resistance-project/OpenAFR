# Prior-art scan — have the top repurposing candidates already been tested against fungi?

Date: 2026-09-04

A companion check to [RESULTS_repurposing.md](../RESULTS_repurposing.md). The screen ranks;
it does not tell you whether someone has *already* run these drugs against a fungus. Before
any wet-lab outreach it is worth knowing which top candidates carry a "we tried that already"
risk. This is that check, run against the three largest public sources of bioactivity /
literature evidence.

> **Bottom line:** across ChEMBL, PubChem BioAssay, and PubMed, **none of the top-10
> candidates has a recorded antifungal test — active or inactive.** Every apparent hit
> resolved to a false positive (see below). This *lowers*, but does not eliminate, the
> "already tried" risk: classic antifungal MIC studies often live only in a paper's
> text/supplement and may be deposited nowhere machine-readable.

## Candidates checked (top novel hypotheses from the screen)

verinurad, vatalanib, flucloxacillin, taranabant, lersivirine, balapiravir/R-1479,
epacadostat, nolatrexed, apatinib. (Ranks per RESULTS_repurposing.md.)

## Source 1 — ChEMBL (deposited bioactivities)

Pulled every deposited activity per molecule; filtered for fungal test organisms /
antifungal assays.

| candidate | ChEMBL id | total activities | fungal-related |
|---|---|---:|---:|
| verinurad | CHEMBL3707347 | 47 | 0 |
| vatalanib | CHEMBL2142861 | 19 | 0 |
| flucloxacillin | CHEMBL1439043 | 8 | 0 |
| taranabant | CHEMBL220360 | 107 | 0 |
| lersivirine | CHEMBL571987 | 190 | 0 |
| balapiravir/R-1479 | CHEMBL577003 | 5 | 0 |
| epacadostat | CHEMBL3545369 | 427 | 0 |
| nolatrexed | CHEMBL2356692 | 30 | 0 |
| apatinib | CHEMBL3186534 | 335 | 0 |

## Source 2 — PubChem BioAssay (incl. inactive/HTS outcomes)

Pulled the full per-compound assay history (`assaysummary`), including inactive outcomes —
this is where deposited *negative* HTS data would surface.

| candidate | CID | total assay records | genuine fungal |
|---|---|---:|---:|
| verinurad | 54767229 | 32 | 0 |
| vatalanib | 151194 | 2298 | 0 |
| flucloxacillin | 21319 | 230 | 0 |
| taranabant | 11226090 | 175 | 0 |
| lersivirine | 16739244 | 285 | 0 |
| balapiravir | 11691726 | 31 | 0 |
| R-1479 | 457388 | 77 | 0 |
| epacadostat | 135564890 | 710 | 0 |
| nolatrexed | 135400184 | 144 | 0 |
| apatinib | 45139106 | 11 | 0 |

Note: a naive filter flagged 2 vatalanib + 6 epacadostat records as "fungal" — all false
positives from the substring **azole** inside chemical names (1,2,5-oxa**diazole**s,
sulfaphen**azole**), not fungi.

## Source 3 — PubMed (literature, via NCBI E-utilities)

Searched each drug (+ development-code synonyms) AND
`(antifungal OR Candida OR Aspergillus OR Cryptococcus OR "azole resistance" OR CYP51 OR
ERG11 OR mycos* OR fungus OR fungal)`. Raw hit counts are misleading; every hit was
triaged. **None reports testing the drug as an antifungal.** The co-occurrences break down as:

- **vatalanib** (15 hits): natural-product papers where compounds were *isolated from* fungi
  (Ganoderma, Aspergillus) and vatalanib is merely a VEGFR assay comparator.
- **flucloxacillin** (507 hits): it is an anti-*bacterial*; hits are Staph/endocarditis, with
  incidental "mycotic aneurysm" phrasing. No antifungal testing.
- **apatinib** (21 hits): closest brush — PMID 41476802 (2025) is a **drug–drug-interaction**
  case report (an osteosarcoma patient *on* apatinib who also received voriconazole for
  invasive aspergillosis). Clinically relevant (azole antifungals inhibit CYP3A → raise TKI
  levels), but **not** a test of apatinib against a fungus.
- **nolatrexed** (6): antifolate/thymidylate-synthase resistance papers; "antifungal" appears
  only as a folate-analog class label.
- **lersivirine / balapiravir / R-1479 / epacadostat / taranabant** (0–12): all HIV / HCV /
  dengue / IDO1 / CB1 pharmacology. No fungal efficacy work.
- **verinurad**: 0 hits.

## Interpretation & recommended next step

- The **published, indexed** "already tried" risk for these ten is low — good news for
  outreach framing.
- **Caveat:** database/abstract absence ≠ never tested. A friendly mycology lab may hold
  unpublished MIC data these sources can't see. The honest move remains: pitch the *ranking
  method* (7 known azoles recovered in the top 3.2%), not any single molecule, and explicitly
  ask collaborators "what have you already ruled out?" — turning their file-drawer negatives
  into external validation of the filter.
- One thing to carry into any apatinib/TKI discussion: these drugs are **co-administered with
  azole antifungals** in oncology (CYP3A DDI). That's a pharmacology footnote, not evidence of
  antifungal activity.

## Reproduce

Mining scripts (kept under the workspace `.context/priorart/`, gitignored per docking-artifact
convention): `mine.py` (ChEMBL activity API), `pubchem.py` (PubChem `assaysummary`),
`pubmed.py` (NCBI E-utilities esearch/esummary). All hit only public REST APIs; no keys.
Queries and organism filters are inlined in each script and summarized above.
