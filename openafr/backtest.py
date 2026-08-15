"""Backtest / honest validation for the C. auris azole early-warning track (issue #27).

Pipeline stage: the credibility test for the whole track (#20-#26).

What this proves, and what it honestly cannot
---------------------------------------------
The track promises an *early warning*: replayed over history, it would flag an
emerging ERG11 azole-resistance mutation (canonically Y132F) before the published
record characterises it, with a positive lead time. Testing that on real data needs
one thing the track never built: an **ERG11 re-caller** that turns SRA reads into an
`erg11_call`. NCBI exposes no such call (spike #20), so on every real snapshot
`erg11_call` is empty and the detector honestly returns "cannot detect yet."

So this module refuses to fake the missing input, and instead measures the three
things that ARE honest (pre-registered in work/PREREGISTRATION_backtest.md):

  * H1-arrival  -- the lead-time *budget*: the deposit lag collection_date ->
                   target_creation_date, over REAL isolates. Bounds any achievable
                   lead time; even a perfect re-caller cannot warn before the data is
                   visible. This is the arrival number spike #20 deferred to #27.
  * H2-null     -- a walk-forward over the REAL snapshot: every step returns the honest
                   empty state, because calls are absent. The correct output, not a
                   pass -- the pipeline refusing to emit a false all-clear.
  * H3-mechanism-- the SAME walk-forward over a labelled (synthetic) replay, to prove
                   the detector converts a first-appearance into a dated, positive-lead-
                   time flag once calls exist -- isolating H2 to the missing re-caller.

Stdlib only, operating on snapshot records (list of dicts) as produced by
openafr.earlywarning. Detection is delegated to openafr.emergence with the frozen
track defaults; this module adds no tunable that could be reverse-fit to a lead time.
"""
import datetime as _dt
import re

from openafr import emergence


# ----------------------------- date helpers ----------------------------------
# NCBI dates come in three shapes: full ISO (2024-03-01), year-month (2024-03),
# and bare year (2024). We parse the first two to a real day; a bare year is NOT a
# day and is treated as undatable *for the lag*, because pinning it to Jan 1 would
# invent up to a year of precision the source never had.

def parse_day(date):
    """A datetime.date from a full ISO ('YYYY-MM-DD') or year-month ('YYYY-MM').

    Returns None for a bare year, an empty value, or anything unparseable -- the
    caller counts these as undatable rather than guessing a day.
    """
    s = (date or "").strip()
    m = re.match(r"^(\d{4})-(\d{2})(?:-(\d{2}))?", s)
    if not m:
        return None
    y, mo, d = int(m.group(1)), int(m.group(2)), int(m.group(3) or 1)
    try:
        return _dt.date(y, mo, d)
    except ValueError:
        return None


# ----------------------- H1: the deposit-lag budget ---------------------------

def deposit_lags(records):
    """Per-isolate deposit lag in days: target_creation_date - collection_date.

    Only isolates with BOTH dates parseable to a day contribute a lag. Returns
    (lags, undatable) where `lags` is a list of ints (may be negative if a sample
    is recorded as collected after it was made visible -- a data quirk we surface
    rather than clamp) and `undatable` counts isolates missing either usable date.
    """
    lags, undatable = [], 0
    for r in records:
        collected = parse_day(r.get("collection_date"))
        visible = parse_day(r.get("target_creation_date"))
        if collected is None or visible is None:
            undatable += 1
            continue
        lags.append((visible - collected).days)
    return lags, undatable


def _quantile(sorted_vals, q):
    """Linear-interpolated quantile of a pre-sorted non-empty list."""
    if not sorted_vals:
        return None
    if len(sorted_vals) == 1:
        return float(sorted_vals[0])
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    frac = pos - lo
    hi = min(lo + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * frac


def lag_summary(records, within_days=365):
    """Summarise the deposit-lag distribution over `records` (H1-arrival).

    Returns a dict:
      n_total, n_datable, n_undatable   -- coverage of the lag measurement
      median, p25, p75, p90             -- lag quantiles in days (None if no data)
      min, max
      within_days, frac_within          -- fraction of datable isolates visible
                                           within `within_days` of collection
      n_negative                        -- lags < 0 (collected-after-visible quirk)
    """
    lags, undatable = deposit_lags(records)
    n_datable = len(lags)
    summary = {
        "n_total": len(records),
        "n_datable": n_datable,
        "n_undatable": undatable,
        "within_days": within_days,
        "median": None, "p25": None, "p75": None, "p90": None,
        "min": None, "max": None,
        "frac_within": None, "n_within": 0, "n_negative": 0,
    }
    if not lags:
        return summary
    s = sorted(lags)
    n_within = sum(1 for x in lags if x <= within_days)
    summary.update({
        "median": _quantile(s, 0.50),
        "p25": _quantile(s, 0.25),
        "p75": _quantile(s, 0.75),
        "p90": _quantile(s, 0.90),
        "min": s[0], "max": s[-1],
        "n_within": n_within,
        "frac_within": n_within / n_datable,
        "n_negative": sum(1 for x in lags if x < 0),
    })
    return summary


# Pre-registered H1 criterion (work/PREREGISTRATION_backtest.md), fixed in advance:
# PASS if median lag <= 365 days AND fraction visible within 365 days >= 0.60.
ARRIVAL_MEDIAN_MAX_DAYS = 365
ARRIVAL_MIN_FRAC_WITHIN = 0.60


def arrival_verdict(summary):
    """Apply the pre-registered H1-arrival criterion to a lag_summary.

    Returns (passed: bool | None, reason: str). `None` when there is no datable
    isolate to judge -- an honest "cannot decide", never a default pass.
    """
    if summary["n_datable"] == 0 or summary["median"] is None:
        return None, "no isolate has both a collection and a visibility day; cannot decide"
    median_ok = summary["median"] <= ARRIVAL_MEDIAN_MAX_DAYS
    frac_ok = summary["frac_within"] >= ARRIVAL_MIN_FRAC_WITHIN
    passed = median_ok and frac_ok
    reason = (f"median lag {summary['median']:.0f}d "
              f"({'<=' if median_ok else '>'} {ARRIVAL_MEDIAN_MAX_DAYS}d) and "
              f"{summary['frac_within']:.0%} visible within {summary['within_days']}d "
              f"({'>=' if frac_ok else '<'} {ARRIVAL_MIN_FRAC_WITHIN:.0%})")
    return passed, reason


# --------------------- H2/H3: the walk-forward replay -------------------------

def creation_date_range(records):
    """(earliest, latest) parseable target_creation_date across records, or (None, None)."""
    days = [parse_day(r.get("target_creation_date")) for r in records]
    days = [d for d in days if d is not None]
    return (min(days), max(days)) if days else (None, None)


def month_steps(start, stop):
    """First-of-month dates from `start`'s month through `stop`'s month, inclusive."""
    if start is None or stop is None:
        return []
    steps = []
    y, m = start.year, start.month
    while (y, m) <= (stop.year, stop.month):
        steps.append(_dt.date(y, m, 1))
        m += 1
        if m > 12:
            y, m = y + 1, 1
    return steps


def walk_forward(records, as_of_dates, target=None, window_days=180, min_count=3,
                 min_delta=0.05, backlog_frac=0.5):
    """Replay `records` as of each date in `as_of_dates`, running the detector.

    At each as-of date the detector sees the full record set but splits it by
    NCBI visibility date exactly as in production, so this is a faithful "what would
    we have said on that day" replay. Records the per-step outcome and, if `target`
    is given, the first as-of date at which `target` is flagged.

    Returns a dict:
      steps        -- list of {as_of, callable(bool), n_signals, target_flagged(bool)}
      first_flag   -- earliest as_of where `target` was flagged, or None
      any_callable -- True if any step had >=1 called isolate in its recent window
    """
    target = target.upper() if target else None
    steps, first_flag, any_callable = [], None, False
    for as_of in sorted(as_of_dates):
        as_of_iso = as_of.isoformat() if hasattr(as_of, "isoformat") else str(as_of)
        result = emergence.detect_emergence(
            records, as_of_iso, window_days=window_days, min_count=min_count,
            min_delta=min_delta, backlog_frac=backlog_frac)
        callable_step = "note" not in result
        any_callable = any_callable or callable_step
        muts = {s["mutation"] for s in result["signals"]}
        hit = target in muts if target else False
        if hit and first_flag is None:
            first_flag = as_of_iso
        steps.append({
            "as_of": as_of_iso,
            "callable": callable_step,
            "n_signals": len(result["signals"]),
            "target_flagged": hit,
        })
    return {"steps": steps, "first_flag": first_flag, "any_callable": any_callable}


def lead_time_days(first_flag, reference_date):
    """reference_date - first_flag, in days (positive => flagged before the reference).

    Returns None if either date is missing or unparseable, so a missing flag never
    silently reads as zero lead time.
    """
    flag = parse_day(first_flag)
    ref = parse_day(reference_date)
    if flag is None or ref is None:
        return None
    return (ref - flag).days
