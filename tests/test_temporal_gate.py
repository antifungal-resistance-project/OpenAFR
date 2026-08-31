"""Known-answer + integrity tests for the temporal (publication-date) holdout gate (issue #82,
follow-up to the disjoint-population holdout #108).

What this holdout adds on top of #108 is the TEMPORAL split, so that is what has to be right:

  * the temporal predicate: a molecule is held out iff its earliest antifungal literature year
    is STRICTLY after the cutoff era (an unknown year is not temporal);
  * the selection is deterministic and genuinely never-docked: it excludes the seed-42 docked
    sample AND the #108 external holdout, and every kept molecule is post-cutoff;
  * the frozen set integrity: the actives are in-domain, all dated after the cutoff, and the
    actives/decoys are disjoint by ChEMBL id from EVERY set the criterion has used (the #108
    sets included), so nothing already seen can leak back in;
  * the pre-registration is frozen (its pinned sha matches the committed file), so the rule,
    cutoff and interpretation cannot be edited after the fact.
"""
import hashlib
import pathlib

from rdkit import Chem, DataStructs

from openafr import inactives
from scripts import make_temporal_holdout as mth
from scripts.make_external_holdout import _id_of, _read_smi
from scripts.validate_gate_external import classify_outcome  # reused verbatim by the temporal gate
from scripts.validate_gate_temporal import PREREG, PREREG_SHA

MIN_AUC = 0.70


def _ids(path):
    return {_id_of(name) for _smi, name in _read_smi(path)}


def test_temporal_predicate_requires_strictly_after_cutoff():
    assert mth.is_temporal(mth.CUTOFF_YEAR + 1) is True
    assert mth.is_temporal(mth.CUTOFF_YEAR) is False          # same year is not "after"
    assert mth.is_temporal(mth.CUTOFF_YEAR - 5) is False
    assert mth.is_temporal(None) is False                     # unknown year is not temporal


def test_cutoff_is_the_frozen_value():
    assert mth.CUTOFF_YEAR == 2023


def test_select_is_deterministic_never_docked_and_post_cutoff():
    """The selection core: same inputs -> same output; excludes docked + #108 sets; post-cutoff."""
    pool = [("C1=CC=NC=C1", f"active_CHEMBL{i}") for i in range(100, 130)]
    years = {f"CHEMBL{i}": 2024 for i in range(100, 130)}
    years["CHEMBL105"] = 2020          # pre-cutoff -> must be excluded
    years["CHEMBL110"] = 2023          # exactly cutoff -> must be excluded
    scored = {"CHEMBL101"}             # in the seed-42 docked sample -> excluded
    external = {"CHEMBL102"}           # in the #108 holdout -> excluded
    a = mth.select_temporal_actives(pool, years, scored, external, n=5)
    b = mth.select_temporal_actives(pool, years, scored, external, n=5)
    assert [n for _s, n in a] == [n for _s, n in b]           # deterministic
    picked = {_id_of(n) for _s, n in a}
    assert picked.isdisjoint({"CHEMBL101", "CHEMBL102", "CHEMBL105", "CHEMBL110"})
    for cid in picked:
        assert years[cid] > mth.CUTOFF_YEAR


def test_select_returns_all_when_fewer_than_n():
    pool = [("C1=CC=NC=C1", f"active_CHEMBL{i}") for i in range(200, 203)]
    years = {f"CHEMBL{i}": 2025 for i in range(200, 203)}
    assert len(mth.select_temporal_actives(pool, years, set(), set(), n=50)) == 3


def test_used_ids_includes_the_external_holdout_sets():
    """The temporal set must be disjoint from #108's external sets, so their ids must be counted
    among the 'already used' ids the builder excludes."""
    used = mth.used_ids()
    for path in ("data/ligands/external_actives.smi", "data/ligands/external_decoys.smi"):
        assert _ids(path) and _ids(path) <= used, path


def test_temporal_actives_are_in_domain_and_post_cutoff():
    years = mth.load_years()
    rows = _read_smi(mth.ACTIVES_OUT)
    assert len(rows) == mth.N_ACTIVES
    for smi, name in rows:
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, name
        assert mol.GetNumHeavyAtoms() <= inactives.MAX_HEAVY, name
        assert inactives.has_aromatic_nitrogen(mol), name
        assert years[_id_of(name)] > mth.CUTOFF_YEAR, name    # every active is newer-literature


def test_temporal_actives_are_never_docked_and_disjoint():
    active_ids = _ids(mth.ACTIVES_OUT)
    assert active_ids.isdisjoint(_ids("data/ligands/verified_actives_sample.smi"))  # seed-42 docked
    assert active_ids.isdisjoint(_ids("data/ligands/external_actives.smi"))          # #108 holdout
    assert active_ids.isdisjoint(_ids(mth.DECOYS_OUT))                               # actives vs decoys
    assert active_ids <= _ids("data/ligands/verified_actives.smi")                   # drawn from the corpus


def test_temporal_decoys_are_property_matched_novel_and_disjoint():
    actives = inactives.load_actives(mth.ACTIVES_OUT)
    active_fps = [a["fp"] for a in actives]
    decoy_ids = _ids(mth.DECOYS_OUT)
    assert decoy_ids.isdisjoint(mth.used_ids())               # untouched by every used set
    for smi, name in _read_smi(mth.DECOYS_OUT):
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, name
        assert inactives.has_aromatic_nitrogen(mol), name
        assert mol.GetNumHeavyAtoms() <= inactives.MAX_HEAVY, name
        fp = inactives.fingerprint(mol)
        assert max(DataStructs.BulkTanimotoSimilarity(fp, active_fps)) < inactives.MAX_TANIMOTO, name


def test_prereg_hash_matches_the_frozen_document():
    got = hashlib.sha256(pathlib.Path(PREREG).read_bytes()).hexdigest()
    assert got == PREREG_SHA


def test_three_outcome_rule_is_the_same_one_108_committed():
    """The temporal gate reuses #108's classify_outcome verbatim — pin that it still behaves."""
    assert classify_outcome(0.80, 0.71, MIN_AUC) == "PASS"
    assert classify_outcome(0.715, 0.630, MIN_AUC) == "SUPPORTED-BUT-WIDE"
    assert classify_outcome(0.69, 0.55, MIN_AUC) == "FAIL"
