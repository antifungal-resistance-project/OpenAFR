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


# --- weighted fit (v2) ----------------------------------------------------------------

def test_wls_equal_weights_matches_ols():
    # With equal weights the weighted fit must reproduce OLS exactly (plus the mx/sxx fields).
    xs = [0, 1, 2, 3]
    ys = [0.1, 0.4, 0.35, 0.8]
    o_int, o_slope, o_s2, o_n = forecast._ols(xs, ys)
    w_int, w_slope, w_s2, w_n, mx_w, sxx_w = forecast._wls(xs, ys, [1.0] * 4)
    assert math.isclose(w_int, o_int) and math.isclose(w_slope, o_slope)
    assert math.isclose(w_s2, o_s2) and w_n == o_n
    assert math.isclose(mx_w, 1.5) and math.isclose(sxx_w, 5.0)


def test_wls_weights_pull_fit_toward_heavy_points():
    # A high-weight point should dominate the fit. Three points collinear at slope 1 plus one
    # off-trend point: up-weighting the off point must move the slope toward it.
    xs = [0, 1, 2, 3]
    ys = [0.0, 1.0, 2.0, 0.0]            # last point breaks the slope-1 line
    _, light, *_ = forecast._wls(xs, ys, [1, 1, 1, 1])
    _, heavy, *_ = forecast._wls(xs, ys, [1, 1, 1, 50])
    assert heavy < light                 # the heavy off-trend point drags the slope down


def test_binomial_weight_rewards_sample_size():
    # Same prevalence, more isolates -> more weight (more information).
    assert forecast._binomial_weight(10, 200) > forecast._binomial_weight(5, 100)
    # Equal N: a near-0 month carries less weight than a mid-prevalence one.
    assert forecast._binomial_weight(1, 400) < forecast._binomial_weight(200, 400)


def test_forecast_weighted_default_differs_from_unweighted_on_uneven_N():
    # Uneven monthly denominators: the weighted (default) and unweighted fits must diverge,
    # proving the default path is actually the weighted one.
    series = [
        {"month": "2026-04", "n_resolved": 100, "n_resistant": 2, "prevalence": 0.02},
        {"month": "2026-05", "n_resolved": 800, "n_resistant": 40, "prevalence": 0.05},
        {"month": "2026-06", "n_resolved": 120, "n_resistant": 12, "prevalence": 0.10},
    ]
    w = forecast.forecast_next(series)                    # weighted default
    u = forecast.forecast_next(series, weighted=False)    # v1 model
    assert not math.isclose(w["slope_logit"], u["slope_logit"])
    assert w["has_interval"] and u["has_interval"]
    assert 0.0 < w["lo"] <= w["point"] <= w["hi"] < 1.0


def test_forecast_unweighted_path_reproduces_v1_ols():
    # v1 reproducibility guard: forecast_next(weighted=False) must match a from-scratch OLS
    # forecast on the Anscombe-corrected logit -- the exact v1 model, untouched by the v2 code.
    series = [
        {"month": f"2026-0{m}", "n_resolved": 300, "n_resistant": n, "prevalence": n / 300}
        for m, n in zip(range(4, 9), (9, 12, 18, 21, 30))
    ]
    u = forecast.forecast_next(series, weighted=False)
    xs = list(range(len(series)))
    ys = [forecast._corrected_logit(e["n_resistant"], e["n_resolved"]) for e in series]
    a, b, _, _ = forecast._ols(xs, ys)
    assert math.isclose(u["point"], forecast._expit(a + b * (xs[-1] + 1)))
    assert math.isclose(u["slope_logit"], b)


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
