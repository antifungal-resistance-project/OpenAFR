"""Known-answer + integrity tests for the external-holdout gate (issue #82).

Two things this gate adds on top of the run-2 grading math have to be right:

  * the pre-registered THREE-outcome rule (PASS / SUPPORTED-BUT-WIDE / FAIL), which turns the
    AUC and its bootstrap CI into the verdict the pre-registration committed to in advance.
  * the anti-inflation rule: a molecule that failed to dock (no pose) must be ranked LAST,
    never dropped — dropping a failed active is the easiest way to fake a passing gate.

Plus integrity checks that the frozen holdout is what the pre-registration pins: the actives
and decoys are in-domain, and are disjoint by ChEMBL id from EVERY set the criterion has used
(so nothing here can leak a molecule the criterion was developed against).
"""
import hashlib
import pathlib

from rdkit import Chem

from openafr import inactives
from scripts.make_external_holdout import USED_ID_FILES, _id_of, _read_smi
from scripts.validate_gate_external import (ACTIVES_FILE, DECOYS_FILE, PREREG, PREREG_SHA,
                                            _ordered_flags, classify_outcome)

MIN_AUC = 0.70


def test_classify_pass_when_whole_ci_clears_bar():
    assert classify_outcome(0.80, 0.71, MIN_AUC) == "PASS"


def test_classify_supported_but_wide_when_ci_dips_below():
    # the actual external result: point estimate over the bar, CI lower bound under it
    assert classify_outcome(0.715, 0.630, MIN_AUC) == "SUPPORTED-BUT-WIDE"


def test_classify_fail_when_point_estimate_below_bar():
    assert classify_outcome(0.69, 0.55, MIN_AUC) == "FAIL"


def test_classify_boundary_exactly_at_bar_is_not_a_fail():
    # AUC exactly at the bar passes the >= bar; a CI floor exactly at the bar is a clean PASS
    assert classify_outcome(0.70, 0.70, MIN_AUC) == "PASS"
    assert classify_outcome(0.70, 0.699, MIN_AUC) == "SUPPORTED-BUT-WIDE"


def test_unrankable_molecules_go_last_never_dropped():
    """A docking-failed active (key None) is placed after every ranked molecule, not removed."""
    rows = [("d1", 3.0, False), ("act_ok", 2.0, True), ("act_fail", None, True),
            ("d2", 2.5, False)]
    flags = _ordered_flags(rows)
    assert len(flags) == 4                 # nothing dropped
    assert flags == [True, False, False, True]   # act_ok, d1, d2, then the failed active last


def test_all_unrankable_actives_sink_to_the_bottom():
    rows = [("a", None, True), ("d", 1.0, False), ("b", None, True)]
    flags = _ordered_flags(rows)
    assert flags[0] is False               # the one rankable decoy leads
    assert flags[1:] == [True, True]       # both failed actives last


def test_prereg_hash_matches_the_frozen_document():
    """The grader refuses to credit a run if the pre-registration was edited after freezing;
    that guarantee is only real if the pinned hash actually matches the committed file."""
    got = hashlib.sha256(pathlib.Path(PREREG).read_bytes()).hexdigest()
    assert got == PREREG_SHA


def _ids(path):
    return {_id_of(name) for _smi, name in _read_smi(path)}


def test_external_actives_are_in_domain():
    for smi, name in _read_smi(ACTIVES_FILE):
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, name
        assert mol.GetNumHeavyAtoms() <= inactives.MAX_HEAVY, name
        assert inactives.has_aromatic_nitrogen(mol), name
    assert len(_read_smi(ACTIVES_FILE)) == 50


def test_external_sets_are_disjoint_by_id_from_every_used_set():
    active_ids = _ids(ACTIVES_FILE)
    decoy_ids = _ids(DECOYS_FILE)
    assert active_ids.isdisjoint(decoy_ids)                 # actives vs decoys
    for used in USED_ID_FILES:
        used_ids = _ids(used)
        # the decoys must touch none of the used sets; the actives are drawn FROM
        # verified_actives, so only that one overlaps them (by construction)
        assert decoy_ids.isdisjoint(used_ids), used
        if "verified_actives" not in used:
            assert active_ids.isdisjoint(used_ids), used


def test_external_actives_are_the_never_scored_subset():
    """Every external active is a verified active that was NOT in the seed-42 scored sample."""
    scored = {_id_of(n) for _s, n in _read_smi("data/ligands/verified_actives_sample.smi")}
    verified = {_id_of(n) for _s, n in _read_smi("data/ligands/verified_actives.smi")}
    for aid in _ids(ACTIVES_FILE):
        assert aid in verified
        assert aid not in scored


def test_external_decoys_are_property_matched_and_novel():
    """The frozen decoy discipline: >=1 aromatic N, in-domain, Tanimoto < 0.35 vs the actives."""
    actives = inactives.load_actives(ACTIVES_FILE)
    active_fps = [a["fp"] for a in actives]
    from rdkit import DataStructs
    for smi, name in _read_smi(DECOYS_FILE):
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, name
        assert inactives.has_aromatic_nitrogen(mol), name
        assert mol.GetNumHeavyAtoms() <= inactives.MAX_HEAVY, name
        fp = inactives.fingerprint(mol)
        assert max(DataStructs.BulkTanimotoSimilarity(fp, active_fps)) < inactives.MAX_TANIMOTO, name
