# Organism-aware phenotypic precedent (#79)

Date: 2026-08-28. Code: [`openafr/precedent.py`](../openafr/precedent.py),
[`scripts/check_precedent.py`](../scripts/check_precedent.py). Tests:
[`tests/test_precedent.py`](../tests/test_precedent.py). Regenerated artifacts:
[`work/candidates/top100_precedent.tsv`](candidates/top100_precedent.tsv),
[`work/candidates/top40_precedent.tsv`](candidates/top40_precedent.tsv),
[`work/candidates/shortlist_scorecard.json`](candidates/shortlist_scorecard.json).

## The bug this fixes

The assay-precedent check (#65) sorts a molecule's ChEMBL records into three buckets:
`on_target` (the CYP51 enzyme — the direct verdict), `phenotypic` (a whole-cell MIC), and
`off_target` (a novelty count). The `phenotypic` bucket is meant to be a *secondary antifungal*
hint — "was this molecule ever active in a fungal cell assay?". But the ChEMBL side routed
**every** non-enzyme µg/mL record into it **regardless of organism**. So a whole-cell MIC against
a *bacterium* or a *mammalian cell line* was counted as if it bore on antifungal activity, and a
`phenotypic tested-active` flag — which demotes a novel candidate in the scorecard
("phenotypic tested-active elsewhere (read assay)") — could be fired by an antibacterial result
that has nothing to say about a fungus. PubChem's side was never affected: it is matched by assay
*name* against fungal terms (`antifungal|candida|aspergillus|…`), so it was already
fungal-gated. **Only the ChEMBL side needed this fix.**

## The fix

`group_by_target` now reads each record's organism (`target_organism`, falling back to
`assay_organism`) and keeps a µg/mL record in `phenotypic` **only if that organism is a fungus**
(`is_fungal_organism`, a curated regex over the clinically-relevant fungal genera plus generic
terms — `candida|aspergillus|cryptococc|…|\bfung|\byeast|…`). Every other µg/mL record — against
a bacterium, a mammalian line, or an **unrecorded** organism — drops to `off_target`: it is still
a measured record (so the novelty count is unchanged in aggregate), but it makes **no antifungal
claim**. An unrecorded organism is treated as non-fungal on purpose: a record we cannot confirm is
fungal must not masquerade as antifungal evidence.

The verdict now **carries the organism it rests on**. `phenotypic_organism_breakdown` returns the
fungal organisms behind a phenotypic verdict and the non-fungal/unspecified organisms that were
reclassified out; `check_precedent.py` prints this as a new `phenotypic_organism` column (fungal
count + examples, or `reclassified: <organisms>` when a molecule's cell records were all
non-fungal), so a reader can see *why* a once-flagged molecule now reads `—`.

Only structured organism fields are read; the free-text `assay_description` is deliberately not
mined, to avoid a stray "vs *Candida* comparator" mention flipping a non-fungal assay to fungal.

## What changed on the shortlist

The one novel candidate this corrects is **flucloxacillin**. Its `phenotypic tested-active` flag —
counted by the scorecard as an open concern — came entirely from a **µg/mL MIC against
*Staphylococcus aureus***: flucloxacillin is an antibacterial β-lactam, and its cell-active
measurement says nothing about antifungal activity. Organism-aware, its phenotypic verdict is
now empty (`reclassified: Staphylococcus aureus`), removing a **false** demotion. (Its other
concerns — a structural ADMET alert, coordinator-fragile domain, protomer state-dependence —
stand; flucloxacillin remains a low-confidence candidate, but for real reasons, not a
misattributed antibacterial one.)

Among the controls the axis behaves as a **positive control on itself**:

- **fluconazole / ravuconazole** (real antifungal azoles) **retain** `phenotypic tested-active`,
  now annotated with the long fungal-organism list (*Candida*, *Aspergillus*, *Cryptococcus*,
  *Trichophyton*, *Fusarium*, … dozens of species) their cell activity was actually measured
  against — exactly as they should.
- **letrozole** (an aromatase inhibitor, not an antifungal) drops from `phenotypic tested-active`
  to no fungal signal (`reclassified: Rattus norvegicus`) — its "active" cell record was a rat
  assay.

## Aggregate effect across the top-100 shortlist

Re-running the full `top100.smi` against live ChEMBL + PubChem (2026-08-29), the phenotypic
column changed on **19 of 100** molecules. The count carrying *any* phenotypic flag fell from
**31 → 17**, and phenotypic `tested-active` from **22 → 8** — i.e. **14 of the 22** old
whole-cell "active" flags were non-fungal noise.

| phenotypic verdict | before | after |
|---|---:|---:|
| tested-active | 22 | 8 |
| tested-inactive | 7 | 5 |
| record-inconclusive | 2 | 4 |
| — (no fungal cell record) | 69 | 83 |

The corrections split into two honest kinds:

- **Confirmed non-fungal, de-flagged (the intended catch).** Antibacterials lose their
  `tested-active`: **cefotiam** (*E. coli*, *K. pneumoniae*, *S. aureus*, …), **tazobactam**
  (*Acinetobacter*, *Bacteroides*, *Campylobacter*, …), **flucloxacillin** (*S. aureus*). So do
  animal/human/viral cell assays: **hydralazine**/**cilostazol** (*Rattus norvegicus*),
  **resiquimod**/**tofacitinib** (*Homo sapiens*/*Mus musculus*), **sangivamycin**/**cytarabine**
  (herpesviruses, primate lines), **allopurinol** (*Pseudomonas*, human, rat). Five of these
  (**azacitidine, cilostazol, cytarabine, letrozole, resiquimod**) flip `tested-active →
  tested-inactive`/`record-inconclusive` rather than to `—`: removing a non-fungal *active*
  record unmasks the molecule's actual *fungal* cell record, which is the more honest verdict.
- **Unrecorded organism, conservatively de-flagged (a documented judgement call).** Five
  molecules (**A-769662, GDC-0879, TG-02, epacadostat, lesinurad**) changed only because their
  µg/mL records carried **no** organism in either structured field. The rule treats an
  unconfirmable organism as non-fungal, so these drop to `—` (shown `reclassified:
  (unspecified)`). This is the conservative direction — it will occasionally hide a genuine but
  unlabelled fungal MIC — chosen so the phenotypic signal never *over*-claims. None of these five
  is in the 16-molecule confidence scorecard, so the one scorecard correction (flucloxacillin)
  rests entirely on a **confirmed** non-fungal organism (*S. aureus*), not on this edge case.

## What this establishes / does not

- **Does:** stop a non-fungal whole-cell assay from firing (or suppressing) an antifungal
  phenotypic verdict, and make every phenotypic verdict carry the organism it rests on — so the
  scorecard's precedent demotions are now organism-honest. It corrects one novel-candidate false
  demotion (flucloxacillin) and validates on the controls (fungal azoles keep the flag; the
  aromatase inhibitor loses it).
- **Does not:** change the enzyme (`on_target`) verdict, which was never organism-ambiguous
  (CYP51 target id is the *C. albicans* enzyme); change any geometry, docking, or ranking; or
  upgrade the phenotypic signal's *strength* — a fungal whole-cell MIC still folds in
  permeability/efflux and remains a hint to read the assays, never an enzyme verdict. The
  file-drawer caveat (#80) still binds: `no-precedent` / `—` is **not** "never tested".

## Reproduce

    conda run -n openafr python -m pytest tests/test_precedent.py -q
    conda run -n openafr python scripts/check_precedent.py \
        work/candidates/top100.smi --near data/ligands/verified_inactives.smi \
        -o work/candidates/top100_precedent.tsv
    # top40 is the first 40 rows (top40.smi == first 40 of top100.smi):
    (head -1 work/candidates/top100_precedent.tsv; sed -n '2,41p' work/candidates/top100_precedent.tsv) \
        > work/candidates/top40_precedent.tsv
    conda run -n openafr python work/candidates/build_scorecard.py

Precedent calls hit ChEMBL (www.ebi.ac.uk) and PubChem (pubchem.ncbi.nlm.nih.gov) live on
2026-08-28.
