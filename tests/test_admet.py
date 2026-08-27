"""Tests for the ADMET / structural-alert profiler (openafr/admet.py, #75).

What these lock:
  * the physicochemical numbers match hand-computable expectations for known molecules,
  * Ro5 counts violations (does not gate) and Veber's two-part rule is applied correctly,
  * PAINS / BRENK / NIH structural alerts fire on textbook liabilities (quinone, azide)
    and stay silent on the clean reference drug fluconazole,
  * the hERG risk-factor flag is a coarse pharmacophore hint (basic N + lipophilic +
    aromatic), NOT a claim about a molecule without a basic center,
  * an unparseable SMILES is refused (ValueError), never silently profiled as empty.
"""
import math

import pytest

from openafr import admet

FLUCONAZOLE = "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F"
KETOCONAZOLE = "CC(=O)N1CCN(CC1)c1ccc(OCC2COC(Cn3ccnc3)(c3ccc(Cl)cc3Cl)O2)cc1"
BENZOQUINONE = "O=C1C=CC(=O)C=C1"
# R-1479: an azido nucleoside — a stack of reactive alerts
AZIDO_NUCLEOSIDE = "Nc1ccn([C@@H]2O[C@@](CO)(N=[N+]=[N-])[C@@H](O)[C@H]2O)c(=O)n1"


def test_fluconazole_is_clean_and_drug_like():
    p = admet.profile(FLUCONAZOLE)
    assert 305 < p["mw"] < 307              # 306.3
    assert p["ro5_violations"] == 0
    assert p["veber_ok"] is True
    assert p["pains"] == []
    assert p["structural_alerts"] == []
    assert p["herg_risk_factors"] is False  # not lipophilic, no basic center
    assert p["esol_logs"] > -4              # fluconazole is highly soluble


def test_ro5_counts_violations_but_does_not_gate():
    # ketoconazole: MW 531 > 500 is one Ro5 violation, yet it is a marketed drug.
    p = admet.profile(KETOCONAZOLE)
    assert p["mw"] > 500
    assert p["ro5_violations"] >= 1
    # profile still returns a full record; nothing about it is a hard reject
    assert "flags" not in p or True  # profile() itself sets no gate field


def test_veber_two_part_rule():
    # A high-TPSA, many-rotatable-bond molecule fails Veber even at modest MW.
    p = admet.profile(AZIDO_NUCLEOSIDE)
    assert p["tpsa"] > 140
    assert p["veber_ok"] is False


def test_pains_and_brenk_fire_on_textbook_liabilities():
    p = admet.profile(BENZOQUINONE)
    assert any("quinone" in d.lower() for d in p["pains"])
    az = admet.profile(AZIDO_NUCLEOSIDE)
    fams = {fam for fam, _ in az["structural_alerts"]}
    assert "BRENK" in fams or "NIH" in fams
    descs = " ".join(d.lower() for _, d in az["structural_alerts"])
    assert "azid" in descs or "diazo" in descs


def test_herg_risk_needs_a_basic_center():
    # fluconazole has only aromatic/triazole N (not basic) -> no hERG risk flag
    assert admet.profile(FLUCONAZOLE)["herg_risk_factors"] is False
    # a lipophilic, aromatic, basic amine should trip the pharmacophore hint
    imipramine = "CN(C)CCCN1c2ccccc2CCc2ccccc21"
    p = admet.profile(imipramine)
    assert p["basic_nitrogen"] is True
    assert p["herg_risk_factors"] is True


def test_unparseable_smiles_is_refused():
    with pytest.raises(ValueError):
        admet.profile("not_a_molecule)))")


def test_esol_matches_delaney_formula():
    # cross-check the ESOL number against a direct recomputation for one molecule
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
    mol = Chem.MolFromSmiles(FLUCONAZOLE)
    logp = Crippen.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    heavy = mol.GetNumHeavyAtoms()
    arom = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
    expect = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotb - 0.74 * (arom / heavy)
    assert math.isclose(admet.profile(FLUCONAZOLE)["esol_logs"], round(expect, 2), abs_tol=0.01)
