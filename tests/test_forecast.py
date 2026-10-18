"""Tests for openafr.forecast -- the resistance-prevalence forecast layer.

Covers the three contracts the page/grader rely on: (1) the series is resolved-denominator,
panel-gated, and date-binned exactly like the trend section; (2) the model degrades honestly
(<2 months -> None, no-residual-d.o.f. -> point without interval, intervals ordered); (3) the
walk-forward backtest never scores a month it trained on and only claims skill when the model
beats BOTH naive baselines.
"""
import math

from openafr import forecast


def _rec(call, source="sra-recaller:called", date="2026-06-15", key="k"):
    return {"isolate_key": key, "fks1_call": call,
            "fks1_resistance_source": source, "target_creation_date": date}


# --- monthly_prevalence_series --------------------------------------------------------

def test_series_empty():
    assert forecast.monthly_prevalence_series([]) == []


def test_series_panel_gating_and_resolved_denominator():
    recs = [
        _rec("S639F", date="2026-06-01", key="a"),          # panel-positive, resolved
        _rec("A100V", date="2026-06-02", key="b"),          # non-panel VUS: resolved, NOT resistant
        _rec("", date="2026-06-03", key="c"),               # re-called wild-type: resolved, NOT resistant
        _rec("S639Y", source="pending:", date="2026-06-04", key="d"),  # NOT resolved -> excluded entirely
    ]
    series = forecast.monthly_prevalence_series(recs)
    assert len(series) == 1
    e = series[0]
    assert e["month"] == "2026-06"
    assert e["n_resolved"] == 3          # a, b, c (d is pending)
    assert e["n_resistant"] == 1         # only the panel token S639F
    assert math.isclose(e["prevalence"], 1 / 3)


def test_series_drops_undated_and_sorts():
    recs = [
        _rec("S639F", date="2026-08-01", key="a"),
        _rec("", date="2026-05-01", key="b"),
        _rec("S639F", date="", key="c"),          # undated -> excluded
        _rec("", date="not-a-date", key="d"),     # unparseable -> excluded
    ]
    series = forecast.monthly_prevalence_series(recs)
    assert [e["month"] for e in series] == ["2026-05", "2026-08"]


# --- forecast_next --------------------------------------------------------------------

def test_forecast_none_below_two_months():
    assert forecast.forecast_next([]) is None
    one = [{"month": "2026-06", "n_resolved": 10, "n_resistant": 1, "prevalence": 0.1}]
    assert forecast.forecast_next(one) is None


def test_forecast_two_months_point_no_interval():
    series = [
        {"month": "2026-06", "n_resolved": 100, "n_resistant": 5, "prevalence": 0.05},
        {"month": "2026-07", "n_resolved": 100, "n_resistant": 7, "prevalence": 0.07},
    ]
    fc = forecast.forecast_next(series)
    assert fc is not None
    assert fc["has_interval"] is False      # 2 points -> no residual d.o.f.
    assert fc["lo"] is None and fc["hi"] is None
    assert 0.0 < fc["point"] < 1.0
    assert fc["next_month"] == "2026-08"
    assert fc["slope_logit"] > 0            # rising


def test_forecast_interval_ordered_and_brackets_point():
    series = [
        {"month": f"2026-0{m}", "n_resolved": 200,
         "n_resistant": n, "prevalence": n / 200}
        for m, n in zip(range(4, 9), (10, 14, 20, 28, 40))
    ]
    fc = forecast.forecast_next(series)
    assert fc["has_interval"] is True
    assert fc["lo"] <= fc["point"] <= fc["hi"]
    assert 0.0 < fc["lo"] < fc["hi"] < 1.0
    assert fc["next_month"] == "2026-09"


def test_forecast_reports_naive_baselines():
    series = [
        {"month": "2026-06", "n_resolved": 100, "n_resistant": 2, "prevalence": 0.02},
        {"month": "2026-07", "n_resolved": 100, "n_resistant": 6, "prevalence": 0.06},
    ]
    fc = forecast.forecast_next(series)
    assert math.isclose(fc["naive_last"], 0.06)
    assert math.isclose(fc["naive_pooled"], 8 / 200)


# --- backtest / has_skill -------------------------------------------------------------

def _clean_rising_series():
    # A deterministically rising, low-noise trend the logit model should fit better than
    # either flat baseline.
    return [
        {"month": f"2026-{m:02d}", "n_resolved": 500,
         "n_resistant": n, "prevalence": n / 500}
        for m, n in zip(range(1, 8), (5, 9, 16, 28, 48, 82, 135))
    ]


def test_backtest_too_short_zero_folds():
    series = [
        {"month": "2026-06", "n_resolved": 100, "n_resistant": 3, "prevalence": 0.03},
        {"month": "2026-07", "n_resolved": 100, "n_resistant": 4, "prevalence": 0.04},
    ]
    bt = forecast.backtest(series, min_train=3)
    assert bt["n_folds"] == 0
    assert forecast.has_skill(bt) is False


def test_backtest_folds_are_out_of_sample_and_ordered():
    bt = forecast.backtest(_clean_rising_series(), min_train=3, horizon=1)
    assert bt["n_folds"] == 4
    for f in bt["folds"]:
        # the target is strictly later than the last training month
        assert f["target_month"] > f["train_end"]


def test_has_skill_true_on_clean_trend():
    bt = forecast.backtest(_clean_rising_series(), min_train=3, horizon=1)
    assert bt["model_mae"] < bt["naive_last_mae"]
    assert bt["model_mae"] < bt["naive_pooled_mae"]
    assert forecast.has_skill(bt) is True


def test_has_skill_false_on_stationary_noise():
    # Prevalence oscillates around a flat mean: the pooled baseline is hard to beat, so the
    # model must NOT claim skill (the real release-710 outcome).
    series = [
        {"month": f"2026-{m:02d}", "n_resolved": 400,
         "n_resistant": n, "prevalence": n / 400}
        for m, n in zip(range(1, 8), (20, 12, 24, 10, 22, 14, 18))
    ]
    bt = forecast.backtest(series, min_train=3, horizon=1)
    assert forecast.has_skill(bt) is False
