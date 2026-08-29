"""Tests for the check_precedent CLI's organism-aware rendering (#79).

The bucket/verdict logic is covered in test_precedent.py; this locks the one thing the CLI adds
on top: the `phenotypic_organism` column, which must (a) name the fungal organisms a phenotypic
verdict rests on, (b) show the non-fungal organisms that were reclassified out when no fungal
signal remains, (c) name a PubChem-sourced (name-matched, no structured organism) verdict as
such rather than as a bare '—', and (d) print '—' only when there is truly no whole-cell record.
"""
from scripts import check_precedent as CP


def _summary(phenotypic, fungal, nonfungal):
    verdict = {"verdict": phenotypic} if phenotypic else None
    return {"phenotypic": verdict,
            "phenotypic_organisms": {"fungal": fungal, "nonfungal": nonfungal}}


def test_fungal_organisms_are_named_with_a_count():
    s = _summary("tested-active", ["Aspergillus fumigatus", "Candida albicans"], [])
    assert CP._pheno_organism(s) == "2 fungal: Aspergillus fumigatus; Candida albicans"


def test_reclassified_nonfungal_organisms_are_flagged():
    """flucloxacillin's real case: the only whole-cell MIC was antibacterial, so no fungal signal."""
    s = _summary(None, [], ["Staphylococcus aureus"])
    assert CP._pheno_organism(s) == "reclassified: Staphylococcus aureus"


def test_pubchem_sourced_verdict_is_named_not_dashed():
    """A phenotypic verdict with an empty ChEMBL breakdown came from PubChem (fungal name-gated)."""
    s = _summary("tested-inactive", [], [])
    assert CP._pheno_organism(s) == "PubChem (fungal name-matched)"


def test_no_whole_cell_record_is_a_dash():
    assert CP._pheno_organism(_summary(None, [], [])) == "—"


def test_long_fungal_list_is_capped():
    many = [f"Candida sp{i}" for i in range(10)]
    out = CP._pheno_organism(_summary("tested-active", many, []))
    assert out.startswith("10 fungal: ") and out.endswith(f"+{10 - CP._ORG_CAP} more")
