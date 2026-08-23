"""Known-answer tests for the per-molecule applicability triage (scripts/applicability_triage.py).

The triage answers "does the validated CYP51 method apply to THIS molecule, and is it novel?"
before any docking. It reuses the pre-registered envelope rules (<=45 heavy atoms, >=1 aromatic
nitrogen) and the PROVENANCE novelty rule (Morgan r=2, Tanimoto < 0.70 vs the 15 validated
actives), and must return exactly one stable verdict per molecule in a fixed priority order.

These lock each verdict branch and the "never silently drop" invariant.
"""
import pytest

rdkit = pytest.importorskip("rdkit")  # triage needs rdkit; skip cleanly without it

from scripts.applicability_triage import NOVELTY_MAX, classify, load_panel

# fluconazole: small azole, aromatic N, IS one of the validated actives -> self-match Tc 1.0.
FLUCONAZOLE = "C1=CC(=C(C=C1F)F)C(CN2C=NC=N2)(CN3C=NC=N3)O"
# posaconazole: 51 heavy atoms -> out-of-domain even though it is a defining azole.
POSACONAZOLE = ("CC[C@@H]([C@H](C)O)N1C(=O)N(C=N1)C2=CC=C(C=C2)N3CCN(CC3)C4=CC=C(C=C4)"
                "OC[C@H]5C[C@](OC5)(CN6C=NC=N6)C7=C(C=C(C=C7)F)F")
# benzene: aromatic but no nitrogen -> out of mechanistic scope.
BENZENE = "c1ccccc1"
# VT-1598: resistance-active tetrazole, 43 heavy, Tanimoto 0.312 -> in-envelope-novel.
VT_1598 = ("C1=CC(=CC=C1COC2=CC=C(C=C2)C#CC3=CN=C(C=C3)C([C@](CN4C=NN=N4)"
           "(C5=C(C=C(C=C5)F)F)O)(F)F)C#N")


@pytest.fixture(scope="module")
def panel():
    return load_panel()


def test_panel_loads_the_fifteen_validated_actives(panel):
    assert len(panel) == 15  # 8 training + 7 holdout


def test_novel_in_envelope_active_is_covered(panel):
    r = classify(VT_1598, panel)
    assert r["verdict"] == "in-envelope-novel"
    assert r["heavy"] == 43
    assert r["max_tc"] < NOVELTY_MAX


def test_member_of_the_panel_is_near_known_not_novel(panel):
    # fluconazole is IN the panel, so max Tanimoto is 1.0 -> a re-discovery, not novel chemistry.
    r = classify(FLUCONAZOLE, panel)
    assert r["verdict"] == "near-known"
    assert r["max_tc"] == pytest.approx(1.0)
    assert r["nearest"] == "fluconazole"


def test_oversized_ligand_is_out_of_domain(panel):
    r = classify(POSACONAZOLE, panel)
    assert r["verdict"] == "out-of-domain"
    assert r["heavy"] == 51


def test_no_aromatic_nitrogen_is_out_of_scope(panel):
    r = classify(BENZENE, panel)
    assert r["verdict"] == "out-of-scope"


def test_out_of_scope_takes_priority_over_domain(panel):
    # A huge molecule with no aromatic N should report out-of-scope (higher priority), not
    # out-of-domain — the criterion is undefined regardless of size.
    big_no_n = "C" * 60  # 60-carbon chain: >45 heavy, zero aromatic nitrogen
    r = classify(big_no_n, panel)
    assert r["verdict"] == "out-of-scope"


def test_unparseable_smiles_is_reported_not_dropped(panel):
    r = classify("this-is-not-a-molecule", panel)
    assert r["verdict"] == "unparseable"
    assert r["heavy"] is None
