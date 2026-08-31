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
import math
import re

from openafr import emergence
from openafr import fks1_caller
from openafr import recaller


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


# --------------------------- azole event frequency ----------------------------
# The number the whole re-caller run exists to produce (TODOS step 2): once `fill` has
# populated erg11_call, what fraction of the isolates the caller could actually place at
# the ERG11 panel carry a KNOWN azole-resistance substitution. It is also the gate on the
# FKS1/v2 decision -- a measured near-zero (not a measured-zero-because-unmeasured) is what
# justifies redirecting to echinocandins.

# resistance_source values that mean the panel was RESOLVED (assessable) for an isolate.
# Only these enter the denominator; partial / failed / pending are excluded, never counted
# as azole-negative -- the same "uncalled is not wild type" honesty the re-caller enforces.
_RESOLVED_SOURCES = ("sra-recaller:called", "sra-recaller:wild-type")


def _resolution_bucket(source):
    """Classify an isolate's resistance_source into resolved / partial / failed / pending."""
    src = (source or "").strip()
    if src in _RESOLVED_SOURCES:
        return "resolved"
    if src.startswith("sra-recaller:partial"):
        return "partial"
    if src.startswith("sra-recaller:failed") or src.startswith("sra-recaller:refused"):
        return "failed"
    return "pending"          # pending:sra-recaller, empty, or any other uncalled state


def wilson_interval(k, n, z=1.96):
    """Wilson score confidence interval for a binomial proportion k/n.

    Why Wilson and not the textbook Wald (p +/- z*sqrt(p(1-p)/n)): the prevalence
    number this run produces will often sit near an extreme (C. auris is largely
    azole-resistant at baseline, so p may be ~0.9) and the first samples are small.
    Wald misbehaves exactly there -- it can emit bounds below 0 or above 1 and its
    coverage collapses when p is near 0/1 -- whereas Wilson stays inside [0,1] and
    keeps close to nominal coverage at small n and extreme p. z=1.96 is 95%.

    Returns (lo, hi); for n == 0 returns None (nothing resolved -> no interval).
    """
    if n == 0:
        return None
    z2 = z * z
    center = (k + z2 / 2) / (n + z2)
    half = (z / (n + z2)) * math.sqrt(k * (n - k) / n + z2 / 4)
    return (max(0.0, center - half), min(1.0, center + half))


def wilson_halfwidth_worst_case(n, z=1.96):
    """Half-width of the Wilson 95% CI at the worst case p=0.5, for run planning.

    'How many isolates before the interval is tight enough?' -- answered *before*
    looking at the data (the same decide-n-before-you-peek discipline as the
    pre-registrations). p=0.5 maximizes the width, so this is the conservative
    planning number; a real p near 0.9 yields a tighter interval than this.
    """
    lo, hi = wilson_interval(n // 2, n, z)
    return (hi - lo) / 2


def panel_prevalence(records):
    """Azole event frequency on a (re-called) snapshot.

    Of the isolates whose ERG11 panel the caller could RESOLVE (resistance_source
    `sra-recaller:called` or `:wild-type`), how many carry a known panel substitution.
    Note `called` != panel-positive: an isolate called only for a non-panel (clade) SNP is
    resolved-but-azole-negative, so positivity is decided from the tokens via
    recaller.is_panel_token, not from the source string.

    Returns a dict:
      n_total            all isolates in the snapshot
      n_resolved         isolates with a resolved panel (the denominator)
      n_panel_positive   resolved isolates carrying >=1 panel substitution (the numerator)
      event_frequency    n_panel_positive / n_resolved, or None if nothing is resolved yet
      ci95               Wilson 95% CI (lo, hi) on event_frequency, or None if nothing resolved
      per_mutation       {panel token -> # resolved isolates carrying it}, sorted by residue
      n_called_nonpanel  resolved isolates called ONLY for non-panel substitution(s)
      excluded           {partial, failed, pending} counts (reported, never counted negative)
    """
    n_resolved = n_positive = n_called_nonpanel = 0
    per_mutation = {}
    excluded = {"partial": 0, "failed": 0, "pending": 0}
    for r in records:
        bucket = _resolution_bucket(r.get("resistance_source"))
        if bucket != "resolved":
            excluded[bucket] += 1
            continue
        n_resolved += 1
        panel_toks = [t for t in emergence.parse_call(r.get("erg11_call"))
                      if recaller.is_panel_token(t)]
        if panel_toks:
            n_positive += 1
            for t in panel_toks:
                per_mutation[t] = per_mutation.get(t, 0) + 1
        elif (r.get("erg11_call") or "").strip():
            n_called_nonpanel += 1
    per_mutation = dict(sorted(per_mutation.items(),
                               key=lambda kv: int(re.search(r"\d+", kv[0]).group())))
    return {
        "n_total": len(records),
        "n_resolved": n_resolved,
        "n_panel_positive": n_positive,
        "event_frequency": (n_positive / n_resolved) if n_resolved else None,
        "ci95": wilson_interval(n_positive, n_resolved),   # None until anything is resolved
        "per_mutation": per_mutation,
        "n_called_nonpanel": n_called_nonpanel,
        "excluded": excluded,
    }


# ------------------------ echinocandin (FKS1) event frequency (v2) ----------------------
# The FKS1 analog of panel_prevalence, and the number the FKS1/v2 track exists to produce:
# once an FKS1 `fill` has populated fks1_call, what fraction of the isolates the caller could
# RESOLVE at the FKS1 panel carry a known echinocandin-resistance substitution (S639F/P/Y).
# Where ERG11 carries ONE status per isolate, FKS1's source is per-window
# (sra-fks1-recaller:HS1=...,HS2=...) because the two hot-spots are called independently -- so
# resolution is judged over the panel-bearing window(s) only.

_FKS1_SOURCE_PREFIX = "sra-fks1-recaller:"
# The windows that carry a known-resistance panel position -- the ones whose coverage decides
# whether the panel is assessable. Derived from the caller so it can never drift (today: HS1;
# HS2 becomes panel-bearing the moment its mutant set is pinned, with no change here).
_FKS1_PANEL_WINDOWS = tuple(name for name, w in fks1_caller.FKS1_WINDOWS.items() if w.panel)


def _fks1_window_statuses(source):
    """Parse a `sra-fks1-recaller:HS1=<st>,HS2=<st>` source into {window: base_status}.

    The base status strips any parenthetical detail (`partial(uncalled:639)` -> `partial`).
    Returns {} for a non-FKS1-recaller source (pending / empty / other) OR a whole-isolate
    failure source (`sra-fks1-recaller:failed(...)`, which carries no per-window part)."""
    src = (source or "").strip()
    if not src.startswith(_FKS1_SOURCE_PREFIX):
        return {}
    out = {}
    for part in src[len(_FKS1_SOURCE_PREFIX):].split(","):
        name, sep, st = part.partition("=")
        if sep:
            out[name.strip()] = st.split("(")[0].strip()
    return out


def _fks1_resolution_bucket(source):
    """Classify an FKS1 source into resolved / partial / failed / pending over the
    panel-bearing windows only.

      resolved -- every panel window was called or wild-type (the panel is assessable)
      failed   -- the whole isolate failed to fetch/align, OR a panel window was refused
                  (an in-window indel)
      partial  -- a panel window was partial (a hot-spot residue uncalled) or missing
      pending  -- not re-called yet (pending:sra-fks1-recaller / :no-reads / empty)

    Only 'resolved' enters the prevalence denominator; the rest are excluded, never counted
    echinocandin-negative -- the same 'uncalled is not wild type' honesty as ERG11."""
    src = (source or "").strip()
    if not src.startswith(_FKS1_SOURCE_PREFIX):
        return "pending"          # pending:sra-fks1-recaller / :no-reads / empty / other
    statuses = _fks1_window_statuses(src)
    if not statuses:
        return "failed"           # sra-fks1-recaller:failed(...) -- no per-window part
    panel_sts = [statuses.get(w) for w in _FKS1_PANEL_WINDOWS]
    if any(s == "refused" for s in panel_sts):
        return "failed"
    if all(s in ("called", "wild-type") for s in panel_sts):
        return "resolved"
    return "partial"              # any panel window missing/partial/unknown


def fks1_panel_prevalence(records):
    """Echinocandin event frequency on an (FKS1-re-called) snapshot -- the FKS1 analog of
    panel_prevalence.

    Of the isolates whose FKS1 panel the caller could RESOLVE (every panel window called or
    wild-type), how many carry a known echinocandin-resistance panel substitution (S639F/P/Y).
    As with ERG11, `resolved` != panel-positive: positivity is decided from the tokens via
    fks1_caller.is_panel_token, so an isolate called only for a non-panel/clade SNP is
    resolved-but-negative. Returns the same shape as panel_prevalence (so a caller/report can
    treat the two genes uniformly). event_frequency/ci95 are None until anything is resolved --
    the honest 'NOT MEASURED YET' state on today's un-filled snapshots."""
    n_resolved = n_positive = n_called_nonpanel = 0
    per_mutation = {}
    excluded = {"partial": 0, "failed": 0, "pending": 0}
    for r in records:
        bucket = _fks1_resolution_bucket(r.get("fks1_resistance_source"))
        if bucket != "resolved":
            excluded[bucket] += 1
            continue
        n_resolved += 1
        panel_toks = [t for t in emergence.parse_call(r.get("fks1_call"))
                      if fks1_caller.is_panel_token(t)]
        if panel_toks:
            n_positive += 1
            for t in panel_toks:
                per_mutation[t] = per_mutation.get(t, 0) + 1
        elif (r.get("fks1_call") or "").strip():
            n_called_nonpanel += 1
    per_mutation = dict(sorted(per_mutation.items(),
                               key=lambda kv: int(re.search(r"\d+", kv[0]).group())))
    return {
        "n_total": len(records),
        "n_resolved": n_resolved,
        "n_panel_positive": n_positive,
        "event_frequency": (n_positive / n_resolved) if n_resolved else None,
        "ci95": wilson_interval(n_positive, n_resolved),   # None until anything is resolved
        "per_mutation": per_mutation,
        "n_called_nonpanel": n_called_nonpanel,
        "excluded": excluded,
    }
