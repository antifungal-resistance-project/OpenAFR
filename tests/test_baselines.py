"""Tests for the multi-baseline comparison scaffold (issue #83).

The scorer binaries (GNINA, PLANTS) are Linux/x86 native and run cloud-gated, so — exactly as
for the GNINA parser — the load-bearing, testable pieces are the pure ones: each opponent's
saved output must map to per-pose ranking keys IN FILE ORDER (aligning 1:1 with
pdbqt.read_poses), the best key must be picked lower-is-better, an unscored pose must degrade to
None rather than a guess, and the grader must rank a scoreless molecule LAST, never drop it.
Plus the pre-registration hash the grader pins must match the frozen document.
"""
import hashlib
import pathlib

from openafr import baselines
from openafr.baselines import (REGISTRY, molecule_key, parse_plants_features, spec_by_key)
from scripts.rescore_baselines import PREREG, PREREG_SHA, _ordered_flags


# --- PLANTS features.csv parser --------------------------------------------------------------

_PLANTS_CSV = """\
LIGAND_ENTRY,TOTAL_SCORE,SCORE_RB_PEN,PLPtotal
lig_entry_00001_conf_01,-92.5,1.2,-95.0
lig_entry_00001_conf_02,-88.1,1.0,-90.0
lig_entry_00001_conf_03,-101.3,1.5,-104.0
"""


def test_plants_parses_total_score_in_pose_order():
    keys = parse_plants_features(_PLANTS_CSV)
    assert keys == [-92.5, -88.1, -101.3]        # file order preserved, TOTAL_SCORE column


def test_plants_column_located_by_header_not_position():
    # TOTAL_SCORE as the 3rd column must still be found by name, not assumed first.
    csv = "NAME,PLPtotal,TOTAL_SCORE\na,-1.0,-50.0\nb,-2.0,-60.0\n"
    assert parse_plants_features(csv) == [-50.0, -60.0]


def test_plants_non_numeric_or_missing_cell_is_none_not_guessed():
    csv = "LIGAND_ENTRY,TOTAL_SCORE\na,-50.0\nb,NA\nc\n"
    assert parse_plants_features(csv) == [-50.0, None, None]


def test_plants_without_total_score_column_returns_empty():
    assert parse_plants_features("NAME,PLPtotal\na,-1.0\n") == []
    assert parse_plants_features("") == []


# --- GNINA keys via the registry parser ------------------------------------------------------

def test_gnina_key_negates_cnnaffinity_so_higher_is_better_becomes_lower_is_better():
    spec = spec_by_key("gnina")
    text = "CNNaffinity: 6.5\nCNNaffinity: 4.0\n"
    assert spec.parse(text) == [-6.5, -4.0]      # best (max cnnaff 6.5) -> smallest key -6.5


# --- molecule_key ----------------------------------------------------------------------------

def test_molecule_key_takes_the_best_lower_is_better():
    assert molecule_key([-92.5, -88.1, -101.3]) == -101.3


def test_molecule_key_ignores_unscored_poses():
    assert molecule_key([None, -3.0, None]) == -3.0


def test_molecule_key_all_unscored_is_none():
    assert molecule_key([None, None]) is None
    assert molecule_key([]) is None


# --- registry --------------------------------------------------------------------------------

def test_registry_names_two_metal_aware_opponents():
    keys = [s.key for s in REGISTRY]
    assert keys == ["gnina", "plants"]
    assert all(s.metal_aware for s in REGISTRY)
    assert all(s.citation for s in REGISTRY)


def test_spec_by_key_rejects_unknown():
    import pytest
    with pytest.raises(KeyError):
        spec_by_key("glide")


# --- grader: scoreless molecules ranked last -------------------------------------------------

def test_ordered_flags_places_scoreless_molecules_last():
    rows = [("d1", -3.0, False), ("act", -5.0, True), ("act_noscore", None, True)]
    assert _ordered_flags(rows) == [True, False, True]   # act, d1, then the scoreless active


# --- pre-registration integrity --------------------------------------------------------------

def test_prereg_hash_matches_the_frozen_document():
    got = hashlib.sha256(pathlib.Path(PREREG).read_bytes()).hexdigest()
    assert got == PREREG_SHA
