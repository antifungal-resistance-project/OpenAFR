"""Tests for the continuous-potency correlation core (openafr/potency.py).

The #81 claim — geometry tracks *how potent* a molecule is — rests on two rules a reviewer
must be able to trust without rereading the grader: how several ChEMBL records collapse to
one continuous potency per compound, and how the rank statistics are computed in a
scipy-free environment. What these lock:
  * median_pchembl uses ONLY exact-relation, allowed-type, numeric pchembl records; a
    censored ('>') record, a wrong standard_type, or a missing value is dropped, not zeroed,
    and a compound with no usable record is placed nowhere (None), never at potency 0;
  * spearman equals scipy's rank correlation, ties included, and is exactly -1 / +1 on
    monotone data;
  * the permutation p-value is one-sided in the pre-registered (negative-rho) direction and
    is small for a clean negative relationship, ~0.5 for noise;
  * bootstrap_ci brackets the point estimate.
No network, no docking; every case is literal.
"""
import numpy as np
import pytest

from openafr import potency


def _rec(**kw):
    base = dict(standard_type="IC50", standard_relation="=", pchembl_value="7.0")
    base.update(kw)
    return base


def test_median_pchembl_basic_median():
    val, n, types = potency.median_pchembl(
        [_rec(pchembl_value="6.0"), _rec(pchembl_value="7.0"), _rec(pchembl_value="9.0")])
    assert val == 7.0 and n == 3 and types == ("IC50",)


def test_median_pchembl_averages_even_count():
    val, n, _ = potency.median_pchembl(
        [_rec(pchembl_value="6.0"), _rec(pchembl_value="8.0")])
    assert val == 7.0 and n == 2


def test_median_pchembl_drops_censored_and_wrong_type_and_missing():
    recs = [
        _rec(pchembl_value="8.0"),                       # kept
        _rec(standard_relation=">", pchembl_value="9.0"),  # censored -> dropped
        _rec(standard_type="MIC", pchembl_value="9.0"),    # wrong type -> dropped
        _rec(pchembl_value=None),                          # missing -> dropped
        _rec(pchembl_value=""),                            # blank -> dropped
    ]
    val, n, types = potency.median_pchembl(recs)
    assert val == 8.0 and n == 1 and types == ("IC50",)


def test_median_pchembl_none_when_no_usable_record():
    val, n, types = potency.median_pchembl([_rec(standard_type="MIC")])
    assert (val, n, types) == (None, 0, ())


def test_median_pchembl_drops_non_numeric_value():
    # a garbage pchembl string must be skipped, not crash, and not counted
    val, n, _ = potency.median_pchembl(
        [_rec(pchembl_value="n/a"), _rec(pchembl_value="6.0")])
    assert val == 6.0 and n == 1


def test_median_pchembl_mixes_allowed_types():
    _, n, types = potency.median_pchembl(
        [_rec(standard_type="Ki", pchembl_value="7.0"),
         _rec(standard_type="Kd", pchembl_value="8.0")])
    assert n == 2 and types == ("Kd", "Ki")


def test_rankdata_average_ties():
    # values 10,10,20 -> ranks 1.5,1.5,3
    assert list(potency.rankdata([10, 10, 20])) == [1.5, 1.5, 3.0]


def test_spearman_perfect_monotone():
    x = [1, 2, 3, 4, 5]
    assert potency.spearman(x, [2, 4, 6, 8, 10]) == pytest.approx(1.0)
    assert potency.spearman(x, [10, 8, 6, 4, 2]) == pytest.approx(-1.0)


def test_spearman_matches_pearson_of_ranks_with_ties():
    x = [1, 2, 2, 3, 5, 8]
    y = [5, 3, 3, 2, 1, 9]
    # independent reference: Pearson on average ranks
    rx, ry = potency.rankdata(x), potency.rankdata(y)
    ref = np.corrcoef(rx, ry)[0, 1]
    assert potency.spearman(x, y) == pytest.approx(ref)


def test_permutation_p_small_for_negative_relationship():
    # tighter N-Fe (small x) with higher potency (large y): strong negative rho
    x = np.linspace(2.0, 4.0, 12)
    y = -x + np.array([0.01 * ((-1) ** i) for i in range(12)])
    p = potency.permutation_p(x, y, n_perm=2000, seed=1, alternative="less")
    assert p < 0.01


def test_permutation_p_midrange_for_noise():
    rng = np.random.default_rng(0)
    x = rng.normal(size=40)
    y = rng.normal(size=40)
    p = potency.permutation_p(x, y, n_perm=2000, seed=7, alternative="less")
    assert 0.2 < p < 0.8


def test_permutation_p_greater_and_two_sided_branches():
    x = np.linspace(2.0, 4.0, 12)
    y = x + np.array([0.01 * ((-1) ** i) for i in range(12)])  # strong POSITIVE rho
    assert potency.permutation_p(x, y, n_perm=2000, seed=2, alternative="greater") < 0.01
    assert potency.permutation_p(x, y, n_perm=2000, seed=2, alternative="two-sided") < 0.02


def test_bootstrap_ci_survives_degenerate_resamples():
    # constant y -> spearman is NaN for any resample touching only ties; the CI code must
    # drop those and still return finite bounds from the non-degenerate resamples.
    x = np.arange(20, dtype=float)
    y = np.concatenate([np.zeros(19), [1.0]])
    lo, hi = potency.bootstrap_ci(x, y, n_boot=500, seed=1)
    assert np.isfinite(lo) and np.isfinite(hi) and lo <= hi


def test_bootstrap_ci_brackets_point_estimate():
    x = np.linspace(2.0, 4.0, 30)
    y = -x + np.random.default_rng(3).normal(scale=0.3, size=30)
    rho = potency.spearman(x, y)
    lo, hi = potency.bootstrap_ci(x, y, n_boot=1000, seed=5)
    assert lo <= rho <= hi
    assert lo < 0  # a real negative relationship's CI sits below zero


def test_spearman_length_mismatch_is_refused():
    # mismatched vectors cannot be ranked against each other -> explicit refusal
    with pytest.raises(ValueError):
        potency.spearman([1.0, 2.0, 3.0], [1.0, 2.0])
