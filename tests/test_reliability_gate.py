"""Known-answer tests for the reliability criterion's aggregation
(scripts/validate_gate_reliability.py).

The reliability gate ranks by the MEDIAN of a molecule's per-seed best N-Fe distances (a
pool-size-neutral estimate of how reliably it reaches the iron), reporting the consensus
minimum for context. These lock that the median is taken (not the min or mean), that
non-coordinating seeds are ignored rather than counted, and that a molecule coordinating in
no seed is unrankable (None) so the gate ranks it last.
"""
from scripts.validate_gate_reliability import reliability_stats


def test_median_not_min_is_the_reliability_value():
    # erratic molecule: one lucky close pose, mostly mediocre -> median stays honest
    byseed = {"42": 2.72, "1": 2.58, "2": 3.70, "3": 2.94, "4": 2.69}
    median, minimum = reliability_stats(byseed)
    assert median == 2.72       # 3rd of sorted [2.58,2.69,2.72,2.94,3.70]
    assert minimum == 2.58      # consensus value, for context


def test_reliable_coordinator_has_low_median():
    byseed = {"42": 2.71, "1": 2.71, "2": 2.70, "3": 2.72, "4": 2.69}
    median, minimum = reliability_stats(byseed)
    assert median == 2.71
    assert minimum == 2.69


def test_non_coordinating_seeds_are_ignored():
    byseed = {"42": None, "1": 3.10, "2": None, "3": 3.20, "4": 3.00}
    median, minimum = reliability_stats(byseed)
    assert median == 3.10       # median of [3.00, 3.10, 3.20], Nones dropped
    assert minimum == 3.00


def test_molecule_that_never_coordinates_is_unrankable():
    median, minimum = reliability_stats({"42": None, "1": None})
    assert median is None       # -> gate ranks it LAST, never dropped
    assert minimum is None


def test_even_count_median_averages_middle_two():
    median, _ = reliability_stats({"42": 2.6, "1": 2.8, "2": 3.0, "3": 3.4})
    assert median == 2.9        # (2.8 + 3.0) / 2
