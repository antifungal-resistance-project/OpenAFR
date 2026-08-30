"""Known-answer tests for the ensemble-consistency classifier
(scripts/candidate_ensemble.py, issue #70).

These lock the per-candidate verdict rule frozen in
work/PREREGISTRATION_candidate_ensemble.md, BEFORE any 5FSA/5V5Z distance was read. The
verdict is a pure function of how many of the considered crystallographic conformers a
molecule coordinates in (>= MIN_COORD_SEEDS/5 seeds each):
  ROBUST             all considered conformers coordinate (strongest hypothesis)
  CONSISTENT         a majority but not all
  CONFORMER-FRAGILE  exactly one (single-snapshot artifact; the informative demotion)
  NON-COORDINATING   none
Plus the conformer_consensus bridge: a conformer with no coordinating seed contributes None
(so it drops out of ensemble_min), never a spurious closest-approach distance.
"""
from scripts.candidate_ensemble import (
    MIN_COORD_SEEDS, classify, conformer_consensus,
)


def test_all_three_conformers_coordinate_is_robust():
    assert classify(n_conf_coord=3, n_conf_considered=3) == "ROBUST"


def test_majority_is_consistent():
    assert classify(n_conf_coord=2, n_conf_considered=3) == "CONSISTENT"


def test_single_conformer_is_fragile():
    # coordinates in only 1 of 3 crystal snapshots -> single-snapshot artifact
    assert classify(n_conf_coord=1, n_conf_considered=3) == "CONFORMER-FRAGILE"


def test_zero_conformers_is_non_coordinating():
    assert classify(n_conf_coord=0, n_conf_considered=3) == "NON-COORDINATING"


def test_robust_requires_all_CONSIDERED_not_a_fixed_three():
    # if the self-docking guard drops a conformer, ROBUST is measured against what remains
    assert classify(n_conf_coord=2, n_conf_considered=2) == "ROBUST"
    assert classify(n_conf_coord=1, n_conf_considered=2) == "CONFORMER-FRAGILE"


def test_no_considered_conformers_is_unrankable():
    assert classify(n_conf_coord=0, n_conf_considered=0) == "UNRANKABLE"


def test_conformer_with_coordinating_seeds_returns_min_distance():
    # 3/5 seeds reach the iron -> value is the MIN over those, n_coord counts them
    byseed = {"42": 2.4, "1": 2.6, "2": 2.5, "3": 5.0, "4": None}
    value, n_coord, n_seeds = conformer_consensus(byseed)
    assert value == 2.4
    assert n_coord == 3
    assert n_seeds == 5


def test_conformer_with_no_coordinating_seed_contributes_none():
    # no seed within FE_CUTOFF -> None so it drops out of ensemble_min, never a stray distance
    byseed = {"42": 4.0, "1": 5.0, "2": None}
    value, n_coord, _ = conformer_consensus(byseed)
    assert value is None
    assert n_coord == 0


def test_coordination_bar_is_the_frozen_min_coord_seeds():
    # a conformer coordinating in exactly MIN_COORD_SEEDS-1 seeds is not enough to count
    assert MIN_COORD_SEEDS == 3
