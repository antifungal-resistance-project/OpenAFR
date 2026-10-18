"""Near-term resistance-prevalence forecast for the C. auris echinocandin weather page.

Pipeline stage: a forecasting LAYER over the same committed FKS1 calls the trend section
reads -- it adds the one thing a "weather" service owes that a hindcast does not, a
*projection* of where resistance prevalence is heading, with an honestly-bounded error
(docs/MOONSHOT.md item 3). Like openafr.backtest, it adds no tunable that could be
reverse-fit to a flattering number, and it refuses to claim a forecast it cannot defend:
the walk-forward backtest (scripts/backtest_forecast.py) must show the model beats a
naive persistence baseline on held-out months, or the page shows the trajectory alone.

What this forecasts, and what it honestly cannot
------------------------------------------------
The time axis is each isolate's NCBI *deposit* date (`target_creation_date`), the only
dated field the committed call table carries. So the series is the prevalence of
echinocandin-resistance among *deposited* genomes per month -- NOT true epidemiological
incidence (collection lags deposit by a ~97-day median, and deposit volume is driven by
sequencing programmes, not transmission). The forecast projects the deposited series
forward; the page states that scope rather than dressing it up as incidence.

Resistance is PANEL-GATED: only carriers of a known-resistance FKS1 substitution
(openafr.fks1_caller.is_panel_token) count, never a hot-spot variant of uncertain
significance -- the same gate the trend/geo sections use, and the distinction the
real-data backtest forced (work/RESULTS_fks1_backtest.md). The denominator is the
*resolved* set (openafr.emergence.resolved), not mutation-carriers -- the coverage fix.

Stdlib only; pure functions over records (list of dicts) or a prevalence series; no I/O,
no network, no wall-clock. The model is a least-squares trend on the logit of monthly
prevalence (with an Anscombe continuity correction so a zero-count month is finite),
projected forward with a Student-t prediction interval -- deterministic, no RNG.
"""
import math

from openafr import emergence, fks1_caller


def _is_resistant(rec, call_field, is_resistant_token):
    """True iff this resolved record carries a counted-resistance token for the gene.

    Panel-gated by `is_resistant_token` (default fks1_caller.is_panel_token): a hot-spot
    variant of uncertain significance is NOT counted, matching the page's resistance headline.
    """
    toks = emergence.parse_call(rec.get(call_field))
    return any(is_resistant_token(t) for t in toks)


def _month(rec, date_field):
    """The YYYY-MM bin for a record, or None when the date is missing/unparseable.

    Deliberately strict: a row we cannot place in time is excluded from the series, never
    bucketed by guess (same discipline as emergence.split_by_creation's `undated`)."""
    raw = (rec.get(date_field) or "")[:7]
    if len(raw) == 7 and raw[:4].isdigit() and raw[4] == "-" and raw[5:].isdigit():
        return raw
    return None


def monthly_prevalence_series(records, *, call_field="fks1_call",
                              date_field="target_creation_date",
                              is_resistant_token=fks1_caller.is_panel_token):
    """Resolved-denominator, panel-gated prevalence per deposit-month, ascending by month.

    Each entry: {month, n_resolved, n_resistant, prevalence}. The denominator is the
    *resolved* subset (openafr.emergence.resolved -- a re-called wild-type isolate counts),
    binned by `date_field`; a month with no placeable resolved isolate does not appear.
    `prevalence` is None only when a month somehow has n_resolved == 0 (it never does here,
    since a month only appears once it has a resolved isolate). Reproducible from the
    committed call table alone.
    """
    buckets = {}
    for rec in emergence.resolved(records, call_field=call_field):
        month = _month(rec, date_field)
        if month is None:
            continue
        n_res, n_pos = buckets.get(month, (0, 0))
        buckets[month] = (n_res + 1,
                          n_pos + (1 if _is_resistant(rec, call_field, is_resistant_token) else 0))
    series = []
    for month in sorted(buckets):
        n_resolved, n_resistant = buckets[month]
        series.append({
            "month": month,
            "n_resolved": n_resolved,
            "n_resistant": n_resistant,
            "prevalence": (n_resistant / n_resolved) if n_resolved else None,
        })
    return series


# --- the model ------------------------------------------------------------------------

def _logit(p):
    return math.log(p / (1.0 - p))


def _expit(x):
    # Guard against overflow for large |x|.
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def _corrected_logit(n_resistant, n_resolved):
    """logit of an Anscombe-corrected prevalence (n+0.5)/(N+1), finite even at n=0 or n=N.

    The 0.5/1 continuity correction keeps a zero-count month on the logit scale instead of
    dropping it, so the trend fit sees the real early-months-were-zero shape."""
    p = (n_resistant + 0.5) / (n_resolved + 1.0)
    return _logit(p)


def _ols(xs, ys):
    """Ordinary least-squares slope/intercept + residual variance for y ~ a + b*x.

    Returns (intercept, slope, s2, n) where s2 is the residual variance with (n-2) d.o.f.
    (0.0 when n <= 2, i.e. no residual d.o.f.). Pure stdlib."""
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx if sxx else 0.0
    intercept = my - slope * mx
    if n > 2:
        resid = [y - (intercept + slope * x) for x, y in zip(xs, ys)]
        s2 = sum(r * r for r in resid) / (n - 2)
    else:
        s2 = 0.0
    return intercept, slope, s2, n


# Student-t 95% two-sided critical values by residual d.o.f. (n-2). A small committed table
# keeps us honest without pulling in scipy; d.o.f. beyond the table use the normal limit 1.96.
_T95 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365,
        8: 2.306, 9: 2.262, 10: 2.228, 15: 2.131, 20: 2.086, 30: 2.042}


def _t_crit(df):
    if df <= 0:
        return None
    if df in _T95:
        return _T95[df]
    for k in sorted(_T95):
        if df <= k:
            return _T95[k]
    return 1.96


def forecast_next(series, *, horizon=1):
    """Project prevalence `horizon` deposit-month(s) past the end of `series`.

    Model: least-squares trend on the Anscombe-corrected logit of monthly prevalence vs a
    0-based month index, inverted to a probability; a Student-t prediction interval (on the
    logit scale, inverted) carries the regression uncertainty. Returned alongside are two
    naive baselines the backtest must be beaten against -- `naive_last` (persistence: the
    last observed prevalence) and `naive_pooled` (the pooled prevalence over all months) --
    so a caller can never dress a no-skill model as a forecast.

    Returns None when fewer than two months are available (a point is not a trajectory, the
    same <2 -> nothing degradation the page's trend section uses). Otherwise a dict:
      {horizon, next_month, point, lo, hi, slope_logit, naive_last, naive_pooled,
       n_months, has_interval}. `lo`/`hi` are None when there are no residual d.o.f.
      (<=2 months): a point with no honest interval, flagged by has_interval=False.
    """
    pts = [e for e in series if e.get("prevalence") is not None]
    if len(pts) < 2:
        return None
    xs = list(range(len(pts)))
    ys = [_corrected_logit(e["n_resistant"], e["n_resolved"]) for e in pts]
    intercept, slope, s2, n = _ols(xs, ys)
    x_star = xs[-1] + horizon
    mu = intercept + slope * x_star

    mx = sum(xs) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    lo = hi = None
    has_interval = False
    tcrit = _t_crit(n - 2)
    if tcrit is not None and s2 > 0 and sxx > 0:
        # Prediction (not mean) interval: adds the 1.0 for a new observation's own scatter.
        se = math.sqrt(s2 * (1.0 + 1.0 / n + (x_star - mx) ** 2 / sxx))
        lo = _expit(mu - tcrit * se)
        hi = _expit(mu + tcrit * se)
        has_interval = True

    total_res = sum(e["n_resolved"] for e in pts)
    total_pos = sum(e["n_resistant"] for e in pts)
    return {
        "horizon": horizon,
        "next_month": _add_months(pts[-1]["month"], horizon),
        "point": _expit(mu),
        "lo": lo,
        "hi": hi,
        "slope_logit": slope,
        "naive_last": pts[-1]["prevalence"],
        "naive_pooled": (total_pos / total_res) if total_res else None,
        "n_months": len(pts),
        "has_interval": has_interval,
    }


def _add_months(month, k):
    """YYYY-MM plus k months, as YYYY-MM."""
    y, m = (int(x) for x in month.split("-"))
    idx = (y * 12 + (m - 1)) + k
    return f"{idx // 12:04d}-{idx % 12 + 1:02d}"


# --- the walk-forward backtest -------------------------------------------------------

def backtest(series, *, min_train=3, horizon=1):
    """Walk-forward honest validation: does the trend model beat naive persistence?

    For each split that leaves at least `min_train` training months and one held-out month
    `horizon` steps ahead, fit on the prefix, forecast the held-out month, and record the
    absolute error of the model AND of both naive baselines against that month's ACTUAL
    prevalence. No month is ever used to both fit and score itself.

    Returns {n_folds, model_mae, naive_last_mae, naive_pooled_mae, interval_coverage,
    folds:[{train_end, target_month, actual, model, naive_last, naive_pooled, in_interval}]}
    -- or {n_folds: 0, ...} when the series is too short to leave one honest fold. The
    PASS/skill decision is the grader's (scripts/backtest_forecast.py), read off the
    pre-registered bar; this function only measures, it never grades.
    """
    pts = [e for e in series if e.get("prevalence") is not None]
    folds = []
    for i in range(min_train, len(pts) - horizon + 1):
        train = pts[:i]
        target = pts[i + horizon - 1]
        fc = forecast_next(train, horizon=horizon)
        if fc is None:
            continue
        actual = target["prevalence"]
        in_interval = (fc["lo"] is not None and fc["hi"] is not None
                       and fc["lo"] <= actual <= fc["hi"])
        folds.append({
            "train_end": train[-1]["month"],
            "target_month": target["month"],
            "actual": actual,
            "model": fc["point"],
            "naive_last": fc["naive_last"],
            "naive_pooled": fc["naive_pooled"],
            "in_interval": in_interval,
        })

    def _mae(key):
        errs = [abs(f["actual"] - f[key]) for f in folds if f[key] is not None]
        return (sum(errs) / len(errs)) if errs else None

    covered = [f for f in folds if f["model"] is not None]
    interval_coverage = (
        sum(1 for f in covered if f["in_interval"]) / len(covered)) if covered else None
    return {
        "n_folds": len(folds),
        "min_train": min_train,
        "horizon": horizon,
        "model_mae": _mae("model"),
        "naive_last_mae": _mae("naive_last"),
        "naive_pooled_mae": _mae("naive_pooled"),
        "interval_coverage": interval_coverage,
        "folds": folds,
    }


def has_skill(bt):
    """True iff the model's held-out MAE beats BOTH naive baselines (the pre-registered skill
    definition). Requires at least one fold; a tie does not count as skill."""
    if not bt or bt.get("n_folds", 0) < 1 or bt.get("model_mae") is None:
        return False
    rivals = [bt.get("naive_last_mae"), bt.get("naive_pooled_mae")]
    rivals = [r for r in rivals if r is not None]
    return bool(rivals) and all(bt["model_mae"] < r for r in rivals)
