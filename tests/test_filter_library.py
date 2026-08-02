"""Known-answer test for the applicability-domain filter (scripts/filter_library.py).

The filter enforces the envelope the gate's evidence is claimed over, before any
docking compute: <= 45 heavy atoms (pre-registered applicability domain) and >= 1
aromatic nitrogen (the N-Fe criterion is undefined otherwise). This locks both
rules and the project's non-negotiable "never silently drop" invariant — every
input molecule must come out either kept or excluded-with-a-reason.
"""
import pytest

rdkit = pytest.importorskip("rdkit")  # filter needs rdkit; skip cleanly without it

from scripts.filter_library import MAX_HEAVY_ATOMS, classify, filter_library

# fluconazole: small azole, 4 aromatic N -> in domain.
FLUCONAZOLE = "C1=CC(=C(C=C1F)F)C(CN2C=NC=N2)(CN3C=NC=N3)O"
# itraconazole: 47 heavy atoms, > 45 -> outside applicability domain.
ITRACONAZOLE = ("CCC(C)N1C(=O)N(C=N1)C2=CC=C(C=C2)N3CCN(CC3)C4=CC=C(C=C4)"
                "OC[C@H]5CO[C@](O5)(CN6C=NC=N6)C7=C(C=C(C=C7)Cl)Cl")
BENZENE = "c1ccccc1"  # aromatic but no nitrogen -> outside mechanistic scope


def test_small_azole_is_kept():
    keep, reason = classify(FLUCONAZOLE)
    assert keep and reason == ""


def test_oversized_azole_is_excluded_by_the_applicability_domain():
    keep, reason = classify(ITRACONAZOLE)
    assert not keep
    assert "applicability domain" in reason
    assert str(MAX_HEAVY_ATOMS) in reason


def test_molecule_without_aromatic_nitrogen_is_out_of_scope():
    keep, reason = classify(BENZENE)
    assert not keep
    assert "aromatic nitrogen" in reason


def test_unparseable_smiles_is_excluded_not_crashed():
    keep, reason = classify("this is not a molecule")
    assert not keep and "unparseable" in reason


def test_nothing_is_silently_dropped():
    lines = [
        f"{FLUCONAZOLE}\tfluconazole",
        f"{ITRACONAZOLE}\titraconazole",
        f"{BENZENE}\tbenzene",
        "not_smiles\tjunk",
    ]
    kept, excluded = filter_library(lines)
    # every well-formed input line is accounted for, exactly once.
    assert len(kept) + len(excluded) == len(lines)
    assert [name for _smi, name in kept] == ["fluconazole"]
    assert {name for name, _ in excluded} == {"itraconazole", "benzene", "junk"}


def test_malformed_line_is_skipped_like_prep_ligands():
    # a line without a tab-separated name is not a molecule record; skip it.
    kept, excluded = filter_library([FLUCONAZOLE, f"{FLUCONAZOLE}\tfluconazole"])
    assert len(kept) == 1 and excluded == []
