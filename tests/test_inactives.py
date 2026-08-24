"""Tests for the verified-inactive decoy construction (openafr/inactives.py).

The benchmark this module builds is the independent dataset the preprint's stopping rule
demands, so its evidence rule has to be exactly as strict as it claims to be. What these
lock:
  * the quantitative bands -- a censored "> X" only counts as inactive once the assay
    reached a real ceiling, and a mid-band MIC is ambiguous, never "inactive by default",
  * contrary evidence always wins: ONE active or ambiguous record anywhere disqualifies a
    compound, because a false inactive silently depresses every metric downstream,
  * nM/enzyme records are read as active evidence only, never as non-inhibition,
  * the single-variable swap -- mechanism (aromatic N), property matching, novelty and the
    applicability domain are enforced identically to the presumed-inactive construction,
  * selection is deterministic and quota-bounded, so the set is reproducible from a cache.
No network; every case is a literal record.
"""
import importlib.util
import pathlib

import pytest
from rdkit import Chem

from openafr import inactives as I

_SCRIPT = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "make_verified_inactives.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("make_verified_inactives", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def rec(**kw):
    """A ChEMBL activity record with only the fields the rule reads."""
    base = {"standard_value": None, "standard_relation": None, "standard_units": None,
            "activity_comment": None, "pchembl_value": None}
    base.update(kw)
    return base


# --- the quantitative evidence bands -------------------------------------------------

@pytest.mark.parametrize("value,relation,expected", [
    ("4", "=", "active"),        # potent MIC
    ("8", "=", "active"),        # exactly the active ceiling is active
    ("16", "=", "ambiguous"),    # weakly active -- NOT usable as an inactive
    ("32", "=", "ambiguous"),    # still inside the band
    ("64", "=", "inactive"),     # a measured, non-inhibitory MIC
    ("128", "=", "inactive"),
    ("32", ">", "inactive"),     # tested to a real ceiling, nothing seen
    ("100", ">", "inactive"),
    ("16", ">", "ambiguous"),    # censored too low to mean anything
    ("4", "<", "active"),
    ("128", "<", "ambiguous"),   # "< 128" is not evidence of activity
])
def test_ugml_bands(value, relation, expected):
    assert I.classify_record(
        rec(standard_value=value, standard_relation=relation, standard_units="ug.mL-1")) == expected


def test_nm_records_are_active_evidence_only():
    """A small nM potency disqualifies; a large one is NOT read as non-inhibition, because
    a weak enzyme IC50 is not a measured non-inhibitory MIC."""
    assert I.classify_record(
        rec(standard_value="50", standard_relation="=", standard_units="nM")) == "active"
    assert I.classify_record(
        rec(standard_value="900000", standard_relation="=", standard_units="nM")) == "ignore"


def test_textual_and_pchembl_verdicts():
    assert I.classify_record(rec(activity_comment="Not Active")) == "inactive"
    assert I.classify_record(rec(activity_comment="inactive")) == "inactive"
    assert I.classify_record(rec(activity_comment="Active")) == "active"
    assert I.classify_record(rec(pchembl_value="6.2")) == "active"
    assert I.classify_record(rec(pchembl_value="4.1")) == "ignore"


def test_unusable_records_are_ignored_not_guessed():
    assert I.classify_record(rec()) == "ignore"
    assert I.classify_record(rec(standard_value="10", standard_units="%")) == "ignore"
    assert I.classify_record(rec(standard_value="not-a-number", standard_units="ug.mL-1")) == "ignore"
    # an unrecognised relation is not silently read as an equality
    assert I.classify_record(
        rec(standard_value="128", standard_relation="??", standard_units="ug.mL-1")) == "ignore"


# --- per-compound aggregation: contrary evidence wins ---------------------------------

def _ugml(v, rel="="):
    return rec(standard_value=v, standard_relation=rel, standard_units="ug.mL-1")


def test_consistent_inactivity_qualifies():
    v, tally = I.compound_verdict([_ugml("64"), _ugml("100", ">"), rec(activity_comment="Not Active")])
    assert v == "verified-inactive"
    assert tally["inactive"] == 3


def test_one_active_record_anywhere_disqualifies():
    v, _ = I.compound_verdict([_ugml("128"), _ugml("128"), _ugml("2")])
    assert v == "has-active-evidence"


def test_one_ambiguous_record_disqualifies():
    """A single weak-but-real MIC means the compound is not a clean inactive."""
    v, _ = I.compound_verdict([_ugml("128"), _ugml("16")])
    assert v == "ambiguous"


def test_no_usable_evidence_is_not_inactive():
    v, _ = I.compound_verdict([rec(), rec(standard_value="1", standard_units="%")])
    assert v == "no-usable-evidence"


def test_ignored_records_do_not_block_a_verdict():
    v, _ = I.compound_verdict([_ugml("64"), rec()])
    assert v == "verified-inactive"


# --- selection: the single-variable swap ----------------------------------------------

FLUCONAZOLE = "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F"
# Two novel-skeleton, aromatic-N, in-domain molecules verified to sit inside fluconazole's
# frozen property tolerances (Tanimoto 0.10 / 0.09) -- stand-ins for real ChEMBL hits.
MATCHED_A = "OCc1nnc(-c2ccc(O)cc2)n1CCOc1ccccc1"
MATCHED_B = "O=C(Nc1ccccc1)c1cnc(N2CCOCC2)nc1"


@pytest.fixture
def actives(tmp_path):
    p = tmp_path / "actives.smi"
    p.write_text(f"{FLUCONAZOLE}\tfluconazole\n")
    return I.load_actives(str(p))


def test_selection_enforces_every_frozen_filter(actives):
    cases = [
        ("CHEMBL1", MATCHED_A),
        ("CHEMBL2", "not a smiles"),                 # unparseable
        ("CHEMBL3", "CC(=O)Oc1ccccc1C(=O)O"),        # no aromatic N (and unmatched)
        ("CHEMBL4", FLUCONAZOLE),                    # identical to the active -> too similar
    ]
    chosen, rejects = I.select_inactives(cases, actives, [a["fp"] for a in actives],
                                         per_active=50)
    assert rejects["unparseable"] == 1
    assert rejects["no-aromatic-n"] == 1
    assert rejects["too-similar"] == 1
    assert [c["chembl_id"] for c in chosen] == ["CHEMBL1"]


def test_over_domain_molecules_are_excluded(actives):
    """The applicability domain (<= 45 heavy atoms) is enforced on decoys too -- otherwise
    the set would contain molecules the method makes no claim about."""
    big = "c1ccccc1" + "".join("C" for _ in range(60))
    chosen, rejects = I.select_inactives([("CHEMBL9", "c1ccncc1" + "C" * 60)], actives,
                                         [a["fp"] for a in actives], per_active=50)
    assert rejects["over-domain"] == 1
    assert chosen == []
    assert Chem.MolFromSmiles(big).GetNumHeavyAtoms() > I.MAX_HEAVY


def test_quota_and_determinism(actives):
    fps = [a["fp"] for a in actives]
    cands = [("CHEMBL1", MATCHED_A), ("CHEMBL2", MATCHED_B)]
    one, _ = I.select_inactives(cands, actives, fps, per_active=1)
    assert len(one) == 1 and one[0]["chembl_id"] == "CHEMBL1"
    again, _ = I.select_inactives(cands, actives, fps, per_active=1)
    assert [c["chembl_id"] for c in one] == [c["chembl_id"] for c in again]


def test_duplicate_smiles_kept_once(actives):
    smi = MATCHED_A
    chosen, rejects = I.select_inactives([("CHEMBL1", smi), ("CHEMBL2", smi)], actives,
                                         [a["fp"] for a in actives], per_active=50)
    assert len(chosen) == 1 and rejects["duplicate"] == 1


def test_property_mismatch_is_rejected(actives):
    """The anti-shortcut filter: docking score tracks heavy-atom count at r = -0.91, so an
    unmatched molecule must not enter the set even when it is measured-inactive and novel."""
    caffeine = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"     # aromatic N, novel, but far too small/polar
    chosen, rejects = I.select_inactives([("CHEMBL5", caffeine)], actives,
                                         [a["fp"] for a in actives], per_active=50)
    assert rejects["no-property-match"] == 1 and chosen == []


def test_unparseable_active_is_refused_not_skipped(tmp_path):
    """A dropped active would silently shrink the positive set and inflate every metric."""
    p = tmp_path / "bad.smi"
    p.write_text("this-is-not-a-smiles\tghost\n")
    with pytest.raises(ValueError):
        I.load_actives(str(p))


def test_malformed_active_lines_are_skipped(tmp_path):
    p = tmp_path / "mixed.smi"
    p.write_text(f"a comment line with no tab\n{FLUCONAZOLE}\tfluconazole\n")
    assert [a["name"] for a in I.load_actives(str(p))] == ["fluconazole"]


def test_matching_tolerances_match_the_presumed_inactive_construction():
    """The swap is single-variable: these are the make_decoys.py constants verbatim."""
    assert (I.MW_TOL, I.LOGP_TOL, I.ROTB_TOL, I.HBD_TOL, I.HBA_TOL) == (25.0, 1.5, 2, 1, 2)
    assert I.MAX_TANIMOTO == 0.35


# --- the CLI shell ---------------------------------------------------------------------

def test_build_refuses_a_missing_cache(tmp_path):
    mod = _load_script()
    with pytest.raises(SystemExit):
        mod.load_cache(str(tmp_path / "nope"))


def test_chembl_sort_key_orders_numerically():
    mod = _load_script()
    assert sorted(["CHEMBL100", "CHEMBL9", "CHEMBL1000"], key=mod._chembl_sort_key) == [
        "CHEMBL9", "CHEMBL100", "CHEMBL1000"]
