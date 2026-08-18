"""Tests for the early-warning backtest (openafr/backtest.py, issue #27).

These lock the guarantees the honest-validation write-up leans on:
  * deposit lag is measured only over isolates with BOTH dates as a real day; a bare
    year or a missing date is undatable and excluded, never guessed,
  * the pre-registered H1-arrival criterion (median <= 365d AND >=60% within 365d) is
    applied exactly, and returns UNDECIDED (not a default pass) with no datable data,
  * the walk-forward flags a target only when the detector does, and reports the first
    such as-of date,
  * on records with empty calls the walk-forward is an honest null (no callable step,
    no flag) -- the pipeline never invents a warning from absent calls,
  * lead time is reference - first_flag, and is None (never silently zero) if a date
    is missing.
No network: everything runs against synthetic snapshot records.
"""
import datetime as _dt

from openafr import backtest as bt


def _rec(key, created, collected, call="", country="Pakistan"):
    return {
        "isolate_key": key, "target_acc": key + ".1",
        "target_creation_date": created, "collection_date": collected,
        "country": country, "erg11_call": call,
        "run_acc": "SRR" + key[-1], "resistance_source": "recalled",
    }


# ---- parse_day ---------------------------------------------------------------

def test_parse_day_iso_and_year_month():
    assert bt.parse_day("2024-03-01") == _dt.date(2024, 3, 1)
    assert bt.parse_day("2024-03") == _dt.date(2024, 3, 1)  # year-month -> first


def test_parse_day_bare_year_and_junk_are_none():
    # A bare year is not a day: pinning it to Jan 1 would invent up to a year.
    assert bt.parse_day("2024") is None
    assert bt.parse_day("") is None
    assert bt.parse_day(None) is None
    assert bt.parse_day("not-a-date") is None
    assert bt.parse_day("2024-13-01") is None  # invalid month


# ---- deposit_lags / lag_summary ----------------------------------------------

def test_deposit_lags_only_fully_datable_contribute():
    recs = [
        _rec("A1", "2024-04-10", "2024-01-01"),  # 100 days
        _rec("A2", "2024-06-01", "2024"),         # collected bare year -> undatable
        _rec("A3", "", "2024-01-01"),             # no creation -> undatable
    ]
    lags, undatable = bt.deposit_lags(recs)
    assert lags == [(_dt.date(2024, 4, 10) - _dt.date(2024, 1, 1)).days]
    assert undatable == 2


def test_lag_summary_quantiles_and_frac_within():
    # lags: 0, 100, 200, 400 days (one beyond the 365 window)
    recs = [
        _rec("B0", "2024-01-01", "2024-01-01"),
        _rec("B1", "2024-04-10", "2024-01-01"),
        _rec("B2", "2024-07-19", "2024-01-01"),
        _rec("B3", "2025-02-04", "2024-01-01"),
    ]
    s = bt.lag_summary(recs, within_days=365)
    assert s["n_datable"] == 4 and s["n_undatable"] == 0
    assert s["min"] == 0 and s["max"] == 400
    assert s["median"] == 150.0  # interp between 100 and 200
    assert s["n_within"] == 3 and abs(s["frac_within"] - 0.75) < 1e-9


def test_lag_summary_surfaces_negative_lag():
    # Collected AFTER made visible: a real data quirk we surface, never clamp.
    recs = [_rec("N1", "2024-01-01", "2024-02-01")]
    s = bt.lag_summary(recs)
    assert s["min"] == -31 and s["n_negative"] == 1


# ---- arrival_verdict (the pre-registered H1 criterion) -----------------------

def test_arrival_verdict_pass():
    recs = [_rec(f"P{i}", "2024-02-01", "2024-01-01") for i in range(10)]  # 31d each
    passed, reason = bt.arrival_verdict(bt.lag_summary(recs))
    assert passed is True and "median lag 31d" in reason


def test_arrival_verdict_fail_on_long_median():
    recs = [_rec(f"L{i}", "2026-01-01", "2024-01-01") for i in range(10)]  # ~2y each
    passed, _ = bt.arrival_verdict(bt.lag_summary(recs))
    assert passed is False


def test_arrival_verdict_undecided_without_datable():
    recs = [_rec("U1", "2024", "2024")]  # both bare years -> no day
    passed, reason = bt.arrival_verdict(bt.lag_summary(recs))
    assert passed is None and "cannot decide" in reason


# ---- month_steps / walk_forward ----------------------------------------------

def test_month_steps_inclusive():
    steps = bt.month_steps(_dt.date(2023, 11, 20), _dt.date(2024, 2, 3))
    assert steps == [_dt.date(2023, 11, 1), _dt.date(2023, 12, 1),
                     _dt.date(2024, 1, 1), _dt.date(2024, 2, 1)]


def test_walk_forward_flags_first_appearance_with_lead_time():
    # Baseline (2022) has no Y132F; three Y132F carriers become visible late 2023.
    recs = [_rec(f"BG{i}", f"2022-0{3 + i}-01", "2021", "K143R") for i in range(3)]
    recs += [
        _rec("Y1", "2023-11-05", "2023-06", "Y132F"),
        _rec("Y2", "2023-12-10", "2023-07", "Y132F"),
        _rec("Y3", "2024-01-15", "2023-08", "Y132F"),
    ]
    start, stop = bt.creation_date_range(recs)
    wf = bt.walk_forward(recs, bt.month_steps(start, stop), target="Y132F")
    assert wf["first_flag"] is not None
    assert wf["any_callable"] is True
    lead = bt.lead_time_days(wf["first_flag"], "2024-09-01")
    assert lead is not None and lead > 0  # flagged before the reference date


def test_walk_forward_is_honest_null_on_empty_calls():
    # Real-data shape: dated, callable-by-reads, but no erg11_call anywhere.
    recs = [_rec(f"E{i}", f"2024-0{1 + i}-01", "2023", "") for i in range(6)]
    start, stop = bt.creation_date_range(recs)
    wf = bt.walk_forward(recs, bt.month_steps(start, stop), target="Y132F")
    assert wf["any_callable"] is False       # no step had a called isolate
    assert wf["first_flag"] is None          # nothing flagged
    assert all(st["n_signals"] == 0 for st in wf["steps"])


# ---- lead_time_days ----------------------------------------------------------

def test_lead_time_positive_and_none():
    assert bt.lead_time_days("2023-08-01", "2024-09-01") == 397
    assert bt.lead_time_days(None, "2024-09-01") is None      # no flag -> not zero
    assert bt.lead_time_days("2023-08-01", "bad") is None


# ---- panel_prevalence: the azole event frequency (the deliverable, TODOS step 2) ----

def _pv(key, source, call=""):
    """A minimal record carrying only what panel_prevalence reads: source + call."""
    return {"isolate_key": key, "resistance_source": source, "erg11_call": call}


def test_prevalence_undefined_when_nothing_resolved():
    # A snapshot the re-caller has not touched: every row pending -> 'cannot measure',
    # never a zero-frequency all-clear.
    recs = [_pv("K1", "pending:sra-recaller"), _pv("K2", "pending:sra-recaller")]
    p = bt.panel_prevalence(recs)
    assert p["n_resolved"] == 0
    assert p["event_frequency"] is None
    assert p["excluded"]["pending"] == 2


def test_prevalence_denominator_is_resolved_only():
    # Only called/wild-type count; partial/failed/pending are excluded, NOT azole-negative.
    recs = [
        _pv("POS1", "sra-recaller:called", "Y132F"),
        _pv("POS2", "sra-recaller:called", "K143R,Y132F"),
        _pv("WT1", "sra-recaller:wild-type", ""),
        _pv("PART", "sra-recaller:partial(panel-uncalled:132)", ""),   # excluded
        _pv("FAIL", "sra-recaller:failed(fetch)", ""),                 # excluded
        _pv("PEND", "pending:sra-recaller", ""),                       # excluded
    ]
    p = bt.panel_prevalence(recs)
    assert p["n_resolved"] == 3                       # POS1, POS2, WT1 only
    assert p["n_panel_positive"] == 2                 # POS1, POS2
    assert abs(p["event_frequency"] - 2 / 3) < 1e-9
    assert p["excluded"] == {"partial": 1, "failed": 1, "pending": 1}
    # per-mutation counts each panel token across resolved carriers, sorted by residue
    assert list(p["per_mutation"].items()) == [("Y132F", 2), ("K143R", 1)]


def test_prevalence_called_nonpanel_is_resolved_but_azole_negative():
    # A 'called' isolate whose only substitution is a non-panel clade SNP is resolved and
    # counts in the denominator, but is NOT a panel hit -- positivity is decided by the
    # token, not the 'called' source.
    recs = [
        _pv("SNP", "sra-recaller:called", "A61T"),      # non-panel substitution
        _pv("POS", "sra-recaller:called", "Y132F"),
    ]
    p = bt.panel_prevalence(recs)
    assert p["n_resolved"] == 2
    assert p["n_panel_positive"] == 1                 # only Y132F
    assert p["n_called_nonpanel"] == 1                # A61T isolate
    assert p["per_mutation"] == {"Y132F": 1}


def test_prevalence_non_panel_mutant_letter_at_panel_position_is_not_a_hit():
    # Y132H sits at a panel RESIDUE but is not the panel MUTANT (Y132F): not azole-positive.
    recs = [_pv("H", "sra-recaller:called", "Y132H")]
    p = bt.panel_prevalence(recs)
    assert p["n_panel_positive"] == 0
    assert p["n_called_nonpanel"] == 1
