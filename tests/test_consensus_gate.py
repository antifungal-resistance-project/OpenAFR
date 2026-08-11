"""Known-answer tests for the consensus criterion's aggregation
(scripts/validate_gate_consensus.py).

The consensus gate re-docks each molecule under several seeds and ranks by the MINIMUM
N-Fe across seeds (a less-biased estimate of the achievable coordination), reporting the
run-to-run spread as confidence. These lock that aggregation and its non-negotiable rules:
the minimum is taken, not the mean (the mean would penalise flexible actives); seeds that
never coordinate are ignored, never counted as coordinating; and a molecule that coordinates
in no seed is unrankable (None), so the gate ranks it last rather than dropping it.
"""
from scripts.validate_gate_consensus import consensus_stats


def test_consensus_takes_the_minimum_not_the_mean():
    # a flexible active: mostly mediocre, occasionally excellent -> consensus is the best
    byseed = {"42": 2.72, "1": 2.58, "2": 3.70, "3": 2.94, "4": 2.69}
    cons, spread = consensus_stats(byseed)
    assert cons == 2.58                      # minimum, not the (worse) mean ~2.93
    assert round(spread, 2) == 1.12          # max - min


def test_steady_molecule_has_small_spread():
    byseed = {"42": 2.71, "1": 2.71, "2": 2.70, "3": 2.72, "4": 2.69}
    cons, spread = consensus_stats(byseed)
    assert cons == 2.69
    assert round(spread, 2) == 0.03


def test_non_coordinating_seeds_are_ignored_not_counted():
    # None means that search found no coordinating pose; it must not lower the estimate
    byseed = {"42": None, "1": 3.10, "2": None, "3": 3.20}
    cons, spread = consensus_stats(byseed)
    assert cons == 3.10
    assert round(spread, 2) == 0.10


def test_molecule_that_never_coordinates_is_unrankable():
    cons, spread = consensus_stats({"42": None, "1": None})
    assert cons is None                      # -> gate ranks it LAST, never dropped
    assert spread == 0.0


def test_single_seed_has_zero_spread():
    cons, spread = consensus_stats({"42": 2.80})
    assert cons == 2.80
    assert spread == 0.0
