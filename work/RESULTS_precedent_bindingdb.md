# Broaden the precedent search to BindingDB (#78)

Date: 2026-08-29. Code: [`openafr/precedent.py`](../openafr/precedent.py),
[`scripts/check_precedent.py`](../scripts/check_precedent.py). Tests:
[`tests/test_precedent.py`](../tests/test_precedent.py) (BindingDB adapter section, 6 new
offline tests). Artifact:
[`work/candidates/shortlist_precedent_bindingdb.tsv`](candidates/shortlist_precedent_bindingdb.tsv)
(the focus shortlist run with `--source all`).

## What this adds and why

The assay-precedent check answers the one question the geometry ranker cannot: *has this molecule
already been measured against fungal CYP51, and what did it say?* Until now it asked two sources —
ChEMBL and PubChem BioAssay. Every source we do **not** check widens the "unknown unknowns": a
deposited measurement those two miss is a candidate we would propose for the bench that someone has
already characterised. #78 adds a third independent source, **BindingDB**.

BindingDB is a curated store of medicinal-chemistry binding affinities (IC50/Ki/Kd, in nM), each
record carrying a literature citation (**PMID/DOI**). It harvests SAR tables from the med-chem
literature that ChEMBL has not necessarily curated, so it is a distinct store of measured
*actives* against the enzyme — exactly the records the other two can miss.

## How it is wired in (and what it will and will not claim)

- **Target-first, not compound-first.** ChEMBL/PubChem are queried per candidate (by InChIKey).
  BindingDB is instead queried *per target*: `getLigandsByUniprots` returns every ligand ever
  measured against the CYP51 enzyme (UniProt **P10613**, *C. albicans* sterol 14α-demethylase —
  the same enzyme the ChEMBL `CHEMBL1780` target names and the docked 5TZ1 structure is). That
  index is built **once** (one call), keyed by InChIKey, then each candidate is matched into it
  locally — no extra network per molecule.
- **On-target only.** Because the query is scoped to the CYP51 UniProt, every BindingDB record
  lands in the `on_target` bucket. BindingDB is not mined for a candidate's *off-target* records
  (that would need the compound-first endpoint), so it only ever adds the direct, high-value enzyme
  verdict the issue asks for.
- **Active evidence only, never a fabricated inactive.** BindingDB affinities are nM potencies, so
  they flow through the *same* nM rule the ChEMBL path already uses (`inactives.classify_record`):
  a value ≤ 10 µM is read as **active** evidence; a large nM value on a binding endpoint is *not* a
  non-inhibitory MIC and is ignored, never counted as inactive. So BindingDB can flip a molecule
  from `no-precedent` to `tested-active`, but it can never invent a `tested-inactive` verdict. It
  merges at the tally level under the same "contrary evidence wins" rule as the other two sources.
- **Provenance carried through.** Each matched record keeps its `monomerid`, `affinity_type`,
  `affinity`, `pmid`, and `doi`. The TSV's `bindingdb` column reports the hit count plus the citing
  PMIDs, so a reader can go read the measurement.

Invoke with `--source all` (adds BindingDB to ChEMBL+PubChem) or `--source bindingdb`. The default
remains `both`, so existing reproductions are unchanged.

## The run (focus shortlist, `--source all`, live)

`python scripts/check_precedent.py work/candidates/shortlist_focus.smi --source all`. The
BindingDB index built to **34 distinct InChIKeys** across the one CYP51 UniProt target. The 16-row
focus shortlist = 10 novel candidates + 6 known-azole controls.

| candidate | verdict | n_records | BindingDB (on-target) |
|---|---|---:|---|
| flutrimazole … LDN-27219, R-1479, L-838417, olprinone, nolatrexed (10 novel) | `no-precedent` / `record-inconclusive` | — | **no BindingDB CYP51 record** |
| ravuconazole* | record-inconclusive | 676 | — |
| letrozole* | record-inconclusive | 3,103 | — |
| **fluconazole*** | tested-active | 21,284 | **8 rec** (PMID 33183866; 34731761; 34772530; 36074863; 36282007; 36439981; +2) |
| **voriconazole*** | tested-active | 6,533 | **2 rec** (PMID 29894182; 39115219) |
| ketoconazole* | tested-active | 10,595 | — |

(`*` = known-azole control. Full table:
[`shortlist_precedent_bindingdb.tsv`](candidates/shortlist_precedent_bindingdb.tsv).)

## What it tells us (honestly)

- **Acceptance criterion met:** the precedent report now cites ≥1 additional independent source
  with per-record provenance — BindingDB IC50/Ki/Kd records against the CYP51 enzyme, each with a
  PMID/DOI.
- **The controls check out.** BindingDB independently carries on-target CYP51 records for
  fluconazole (8) and voriconazole (2) — the class the tool already flags `tested-active`. This is
  the source behaving correctly on molecules we *know* were measured against the enzyme.
- **The novel candidates stay clean.** None of the 10 novel candidates gained a BindingDB CYP51
  binding record. On this shortlist BindingDB did not flip any candidate's verdict — it *confirms*
  the `no-precedent` reads rather than overturning them, and it shrinks (does not eliminate) the
  unknown-unknown gap for these specific molecules against this one enzyme.
- **A miss is still weak.** As always: no BindingDB record ≠ never tested. BindingDB covers the
  literature it has curated against P10613; a failing or unpublished measurement can still evade
  all three sources. The verdict remains `no-precedent`, never `untested` (see the module
  docstring and [`METHODS_file_drawer.md`](METHODS_file_drawer.md)).

## Extensibility

`CYP51_UNIPROTS` is frozen to P10613 (the docked *C. albicans* target) for auditability, mirroring
the frozen ChEMBL `CYP51_TARGETS`. Adding a genuine CYP51 orthologue (e.g. *C. auris* or
*A. fumigatus* ERG11) is a one-line change to that set — but it becomes a *cross-species* on-target
claim, so it should be made deliberately and noted in the run, alongside the organism-aware
phenotypic handling from [#79](RESULTS_precedent_organism.md).
