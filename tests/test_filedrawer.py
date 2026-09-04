"""Tests for the file-drawer bound (openafr/filedrawer.py, issue #80).

What these lock:
  * the leak is 1 - deposit_rate, and out-of-range probabilities are refused,
  * the residual is a BAND (higher deposit rate -> lower leak, so the band flips), and there is
    NO point estimate exposed -- false precision is designed out,
  * adding independent sources compounds the miss geometrically (the issue-#78 payoff),
  * the inline caveat states the leak honestly: it names the note, calls itself illustrative,
    and never claims to be a probability the compound failed.
"""
import pytest

from openafr import filedrawer as fd


# --- the single-source leak -----------------------------------------------------------

def test_leak_is_one_minus_deposit_rate():
    assert fd.leak(0.0) == 1.0
    assert fd.leak(1.0) == 0.0
    assert fd.leak(0.6) == pytest.approx(0.4)


@pytest.mark.parametrize("bad", [-0.01, 1.01, 2.0, -1.0])
def test_leak_refuses_non_probabilities(bad):
    with pytest.raises(ValueError):
        fd.leak(bad)


# --- the residual is a band, not a point ----------------------------------------------

def test_residual_band_flips_deposit_to_leak():
    # default deposit band [0.40, 0.70] -> leak band [0.30, 0.60]
    lo, hi = fd.residual_hidden_band()
    assert lo == pytest.approx(1 - fd.DEPOSIT_RATE_HI)
    assert hi == pytest.approx(1 - fd.DEPOSIT_RATE_LO)
    assert lo < hi  # a genuine band


def test_residual_band_rejects_disordered_or_out_of_range():
    with pytest.raises(ValueError):
        fd.residual_hidden_band(lo=0.7, hi=0.4)   # lo > hi
    with pytest.raises(ValueError):
        fd.residual_hidden_band(lo=-0.1, hi=0.5)  # out of range


def test_no_point_estimate_is_exposed():
    # false precision is designed out: there is deliberately no scalar residual accessor.
    assert not hasattr(fd, "residual_hidden_point")


# --- adding sources compounds the miss (issue #78 payoff) -----------------------------

def test_adding_one_source_halves_leak_at_half_deposit():
    base = fd.residual_hidden_band()
    after = fd.leak_after_adding_sources([0.5])
    assert after[0] == pytest.approx(base[0] * 0.5)
    assert after[1] == pytest.approx(base[1] * 0.5)


def test_more_sources_only_shrink_the_leak():
    base = fd.residual_hidden_band()
    one = fd.leak_after_adding_sources([0.5])
    two = fd.leak_after_adding_sources([0.5, 0.5])
    assert two[1] < one[1] < base[1]
    assert two[0] < one[0] < base[0]


def test_adding_no_sources_returns_base_band():
    base = fd.residual_hidden_band()
    assert fd.leak_after_adding_sources([]) == pytest.approx(base)


# --- the inline caveat is honest ------------------------------------------------------

def test_caveat_line_is_honest_and_cited():
    line = fd.caveat_line()
    assert fd.METHODS_NOTE in line              # cites the note
    assert "illustrative" in line.lower()       # flags itself as assumption-driven
    assert "not 'untested'" in line             # never implies novelty
    # the default 30-60% band shows up as percentages
    assert "30" in line and "60" in line


# --- adding sources refuses a malformed base band --------------------------------------

@pytest.mark.parametrize("bad", [(-0.1, 0.5), (0.6, 0.4), (0.5, 1.2), (1.2, 0.5)])
def test_leak_after_adding_sources_refuses_bad_base_band(bad):
    with pytest.raises(ValueError):
        fd.leak_after_adding_sources([0.3], base_band=bad)
