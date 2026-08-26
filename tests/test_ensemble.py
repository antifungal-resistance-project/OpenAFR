"""Locks openafr.ensemble — the aggregation the within-azole ensemble gate rests on."""
import pathlib
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr.ensemble import aggregate, ensemble_rows


def test_min_picks_smallest_across_conformers():
    assert aggregate([2.7, 2.4, 2.9], "min") == 2.4


def test_mean_averages_present_values():
    assert aggregate([2.0, 3.0], "mean") == pytest.approx(2.5)


def test_none_values_are_ignored_not_treated_as_zero():
    # A conformer with no coordinating pose must not pull the aggregate toward 0.
    assert aggregate([None, 2.6, None], "min") == 2.6
    assert aggregate([None, 2.0, 4.0], "mean") == pytest.approx(3.0)


def test_all_none_is_unrankable():
    assert aggregate([None, None], "min") is None
    assert aggregate([None, None], "mean") is None


def test_unknown_mode_raises():
    with pytest.raises(ValueError):
        aggregate([1.0], "median")


def test_ensemble_rows_min_over_conformers():
    per = {"5TZ1": {"a": 2.8, "b": 3.0}, "5FSA": {"a": 2.5, "b": None}}
    rows = dict(ensemble_rows(["a", "b"], per, "min"))
    assert rows["a"] == 2.5      # min(2.8, 2.5)
    assert rows["b"] == 3.0      # min(3.0, ignore None)


def test_self_pair_exclusion_drops_named_conformer():
    # itraconazole must not be scored against 5V5Z (its own co-crystal).
    per = {"5TZ1": {"itra": 2.9}, "5V5Z": {"itra": 2.1}}
    excl = {"itra": {"5V5Z"}}
    rows = dict(ensemble_rows(["itra"], per, "min", exclude=excl))
    assert rows["itra"] == 2.9   # the leaky 2.1 from 5V5Z is excluded


def test_missing_molecule_in_a_conformer_is_treated_as_none():
    per = {"5TZ1": {"a": 2.7}, "5FSA": {}}   # 'a' absent from 5FSA
    rows = dict(ensemble_rows(["a"], per, "min"))
    assert rows["a"] == 2.7
