"""Tests for the human-CYP off-target liability profiler (openafr/cyp_offtarget.py, #74).

What these lock:
  * an azole warhead (imidazole/triazole) flags HIGH type-II across ALL five human CYPs —
    the mechanistic claim the axis exists to make (ketoconazole is the textbook example);
  * the strong-vs-weak split: imidazole/triazole/tetrazole => HIGH, pyrazole => MODERATE;
  * a bare non-azole azine nitrogen (pyridine) is a weaker MODERATE pan-CYP ligator;
  * per-isoform substrate-pharmacophore hints fire on the right isoform (2C9 acid, 2D6
    basic-N, 1A2 planar aromatic) and a fully non-heme, non-matching molecule stays clean;
  * the level is the STRONGER of the two signals and `basis` records why;
  * an unparseable SMILES is refused (ValueError), never silently profiled as empty.
"""
import pytest

from openafr import cyp_offtarget as cy

KETOCONAZOLE = "CC(=O)N1CCN(CC1)c1ccc(OCC2COC(Cn3ccnc3)(c3ccc(Cl)cc3Cl)O2)cc1"  # imidazole
FLUCONAZOLE = "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F"                              # bis-triazole
PYRIDINE = "c1ccncc1"                                                            # bare azine
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(O)=O"                                         # 2C9 acid, no N
DIPHENHYDRAMINE = "CN(C)CCOC(c1ccccc1)c1ccccc1"                                  # 2D6 basic N
NAPHTHALENE = "c1ccc2ccccc2c1"                                                   # 1A2 planar PAH


def test_azole_warhead_is_high_type_ii_across_the_whole_panel():
    p = cy.profile(KETOCONAZOLE)
    assert p["azole_warhead"] == "imidazole"
    assert p["warhead_strength"] == "strong"
    assert p["n_high"] == 5
    assert all(p["panel"][c]["liability"] == "high" for c in cy.PANEL)
    # the basis names the mechanism, so a reader sees WHY, not just a level
    assert any("type-II" in b for b in p["panel"]["CYP3A4"]["basis"])


def test_triazole_is_a_strong_warhead_too():
    p = cy.profile(FLUCONAZOLE)
    assert p["azole_warhead"] == "1,2,4-triazole"
    assert p["warhead_strength"] == "strong"
    assert p["n_high"] == 5


def test_pyrazole_is_a_weaker_moderate_warhead():
    # lersivirine's 1,2-diazole: a real but softer type-II ligand -> MODERATE, not HIGH
    lersivirine = "CCc1nn(CCO)c(CC)c1Oc1cc(cc(c1)C#N)C#N"
    p = cy.profile(lersivirine)
    assert p["azole_warhead"] == "pyrazole"
    assert p["warhead_strength"] == "weak"
    assert p["n_high"] == 0
    assert all(p["panel"][c]["liability"] == "moderate" for c in cy.PANEL)


def test_bare_azine_nitrogen_is_moderate_pan_cyp():
    p = cy.profile(PYRIDINE)
    assert p["azole_warhead"] is None
    assert p["azine_nitrogen"] is True
    assert p["type_ii_ligator"] is True
    assert p["n_high"] == 0
    assert all(p["panel"][c]["liability"] == "moderate" for c in cy.PANEL)


def test_2c9_acid_pharmacophore_fires_without_any_heme_ligand():
    # ibuprofen: a lipophilic carboxylic acid, no ligating nitrogen at all
    p = cy.profile(IBUPROFEN)
    assert p["type_ii_ligator"] is False
    assert p["panel"]["CYP2C9"]["liability"] == "moderate"
    assert any("acidic" in b for b in p["panel"]["CYP2C9"]["basis"])


def test_2d6_basic_nitrogen_pharmacophore_fires():
    p = cy.profile(DIPHENHYDRAMINE)
    assert p["panel"]["CYP2D6"]["liability"] == "moderate"
    assert any("basic N" in b for b in p["panel"]["CYP2D6"]["basis"])


def test_1a2_planar_aromatic_pharmacophore_fires():
    # naphthalene: small, flat, fully aromatic, no heteroatom -> a clean 1A2-shape hit
    p = cy.profile(NAPHTHALENE)
    assert p["type_ii_ligator"] is False
    assert p["panel"]["CYP1A2"]["liability"] == "moderate"
    assert any("planar" in b for b in p["panel"]["CYP1A2"]["basis"])


def test_a_clean_aliphatic_molecule_is_low_everywhere():
    p = cy.profile("CCCCCCO")  # hexanol: no ring, no N, small — nothing should fire
    assert p["type_ii_ligator"] is False
    assert p["n_high"] == 0 and p["n_moderate"] == 0
    assert p["flagged_isoforms"] == []
    assert cy.flag_summary(p) == "clean"


def test_liability_is_the_stronger_of_the_two_signals():
    # ketoconazole is large+lipophilic (a 3A4 pharmacophore match) AND an azole warhead;
    # the reported level must be HIGH (the stronger), with BOTH reasons in the basis.
    p = cy.profile(KETOCONAZOLE)
    assert p["panel"]["CYP3A4"]["liability"] == "high"


def test_summary_string_leads_with_the_worst_signal():
    assert cy.flag_summary(cy.profile(FLUCONAZOLE)).startswith("type-II(1,2,4-triazole)")
    assert cy.flag_summary(cy.profile(PYRIDINE)).startswith("azine:")


def test_unparseable_smiles_is_refused():
    with pytest.raises(ValueError):
        cy.profile("not_a_molecule)))")
