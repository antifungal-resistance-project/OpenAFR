"""Is the FKS1 fill due? — the non-billable freshness check behind the weather cadence.

The resistance-weather page is a thin view over a live NCBI *C. auris* pull with the durable,
git-tracked FKS1 call table (`data/earlywarning/fks1_calls/`) overlaid on top. The calls are
manufactured out-of-band by the billable GCP re-caller (issue #171,
`work/RUNBOOK_fks1_fill.md`); nothing re-runs that fill on a cadence, so as NCBI deposits new
in-window genomes past the last fill's coverage the page quietly degrades toward "awaiting
re-call". The only signal that a (paid) refill is worth running is that degradation — and nobody
is watching the page.

This module answers "is a refill due, and how big is the delta?" from a *fresh, non-billable*
record set (a live pull or a snapshot) overlaid with the newest committed table. It does NO
network or tool work of its own and never invokes the billable fill; it just measures the gap so
a scheduled watchdog (`scripts/fill_due_check.py`) can ring the bell precisely when the paid step
pays off, and stay quiet otherwise.

The incremental-delta subtlety: a fresh NCBI pull carries NO echinocandin call (NCBI runs no AMR
pipeline), so every row arrives `pending:`. Counting pending rows on the raw pull would say "the
whole window is uncalled" on every run. So we first overlay the newest committed call table
(replaying the billable calls onto the matching isolates), and only THEN count the in-window
isolates still pending with a run accession — the isolates a refill would genuinely add.
"""
import datetime
import os

from openafr import earlywarning as ew
from openafr import emergence

# The resistance-source sentinel a not-yet-re-called FKS1 isolate carries (mirrors
# scripts/recall_fks1.FKS1_PENDING; kept local so this layer takes no dependency on the CLI).
FKS1_PENDING = "pending:sra-fks1-recaller"

# A re-called delta worth paying the billable GCP fill for. Below this the weekly publish keeps
# showing the last fill honestly (with its freshness clause) rather than re-paying for a window
# that has barely moved; the runbook's "monthly is enough" cadence assumes a signal like this.
DEFAULT_THRESHOLD = 25

# How much runway to demand before the page would go dark. The watchdog runs monthly, so 90 days
# is ~3 checks of warning -- enough to schedule and run the billable GCP fill before the committed
# window rolls fully off the page. This is the "never a surprise" lead time: the bell rings while
# there is still time to act, not on the month the page actually goes blank.
DEFAULT_DARK_LEAD_DAYS = 90


def _parse_iso_date(value):
    """Parse an ISO date (YYYY-MM-DD, extra time suffix tolerated) -> date, or None."""
    s = str(value or "")[:10]
    try:
        return datetime.date.fromisoformat(s)
    except ValueError:
        return None


def _window_aged_out(committed_as_of, committed_window_days, current_as_of):
    """True when the last fill's entire window has rolled off the displayed window.

    The committed fill covered the window (committed_as_of, committed_as_of + window]; the page
    now displays (current_as_of, current_as_of + window]. Once the committed window's END is at
    or before the current window's START the two no longer overlap, so none of the committed
    calls land on a displayed isolate and the page shows "awaiting re-call" regardless of how
    many new isolates exist. Unknown committed dates (e.g. a --no-window fill) -> None (can't
    tell); the caller treats that as "rely on the delta count instead"."""
    committed = _parse_iso_date(committed_as_of)
    start = _parse_iso_date(current_as_of)
    try:
        window_days = int(committed_window_days)
    except (TypeError, ValueError):
        return None
    if committed is None or start is None:
        return None
    committed_end = committed + datetime.timedelta(days=window_days)
    return committed_end <= start


def projected_dark_date(committed_as_of, committed_window_days):
    """The calendar date the page goes dark if no newer fill lands -- computed in ADVANCE.

    The page displays the trailing window `(start, start + window]` with `start = today - window`.
    The committed fill's calls live in `(committed_as_of, committed_as_of + window]`; once the
    display window's START reaches the committed window's END the two stop overlapping and the
    page shows "awaiting re-call". So `start = committed_as_of + window`, i.e.
    `today = committed_as_of + 2*window`. Returns that `date`, or None when the committed date or
    window is unknown (e.g. a `--no-window` fill) -- same "can't tell" contract as
    `_window_aged_out`."""
    committed = _parse_iso_date(committed_as_of)
    try:
        window_days = int(committed_window_days)
    except (TypeError, ValueError):
        return None
    if committed is None:
        return None
    return committed + datetime.timedelta(days=2 * window_days)


def display_window(calls_dir, fallback_as_of, fallback_window):
    """The (as_of, window_days) the PAGE should display, anchored to the latest committed fill.

    The page shows a frozen, git-tracked fill whose called isolates sit in a fixed date range.
    A trailing-`today` window marches forward and slowly ages those isolates out, so the page
    silently decays toward "awaiting re-call" even though the fill is the freshest data we have.
    Anchoring the DISPLAY window to the fill's own covered window `(fill_as_of, fill_as_of +
    window_days]` keeps the page showing the fill at full coverage, with an honest "calls current
    through <fill date>" clause, until a newer fill is committed. Returns (as_of, window_days) as
    (str, int); falls back to (fallback_as_of, fallback_window) when no committed fill exists.

    This is DISPLAY-only. `fill_due` deliberately stays `today`-relative: its job is to measure the
    real new-deposit gap against now and decide when the billable refill pays off."""
    idx = ew.latest_fks1_calls_index_entry(calls_dir)
    if not idx:
        return fallback_as_of, fallback_window
    as_of = idx.get("as_of") or None
    try:
        window_days = int(idx.get("window_days"))
    except (TypeError, ValueError):
        window_days = None
    if as_of is None or window_days is None:
        return fallback_as_of, fallback_window
    return as_of, window_days


def fill_due(records, calls_dir, as_of, window_days=180, threshold=DEFAULT_THRESHOLD,
             dark_lead_days=DEFAULT_DARK_LEAD_DAYS):
    """Measure the FKS1 fill gap for `records` (a fresh pull/snapshot) against the committed
    calls in `calls_dir`, scoped to the trailing-`window_days` window ending after `as_of`.

    Overlays the newest committed call table onto `records` IN PLACE (so counts reflect the real
    echinocandin picture, not an all-pending raw pull), then reports the in-window coverage and
    the incremental re-call delta. The dark-date forecast is derived from `as_of` (the watchdog
    runs with `as_of = today - window`, so `as_of` tracks the clock 1:1 and the forecast needs no
    separate wall-clock); `dark_lead_days` is how much runway to demand before the page goes dark.
    Returns a dict:

      due              -- bool: is a (billable) refill worth running now?
      reason           -- short human string naming why (or why not)
      uncalled_new     -- in-window isolates still pending WITH a run accession (the refill's
                          actual workload == new SRA downloads)
      total_in_window  -- isolates NCBI-visible in the displayed window
      resolved_in_window / coverage_pct -- how much of the window the committed calls cover now
      window_aged_out  -- bool|None: the last fill's whole window has rolled off the page
      dark_date        -- ISO str|None: the date the page goes dark if no newer fill lands
      days_until_dark  -- int|None: calendar days from `today` to `dark_date` (<=0 == already dark)
      dark_lead_days   -- the runway bar applied to the dark-date trigger
      committed_as_of / committed_window_days / newest_release -- the last fill's provenance
      threshold        -- the uncalled_new bar that was applied
    """
    table = ew.latest_fks1_call_table(calls_dir)
    if table and os.path.exists(table):
        ew.overlay_fks1_calls(records, table)

    _baseline, recent, _undated = emergence.split_by_creation(records, as_of, window_days)
    total_in_window = len(recent)
    resolved_in_window = len(emergence.resolved(recent, "fks1_call"))
    uncalled_new = sum(
        1 for r in recent
        if (r.get("fks1_resistance_source") or "") == FKS1_PENDING
        and (r.get("run_acc") or "").strip()
    )
    coverage_pct = (100.0 * resolved_in_window / total_in_window) if total_in_window else 0.0

    idx = ew.latest_fks1_calls_index_entry(calls_dir)
    committed_as_of = (idx or {}).get("as_of") or None
    committed_window_days = (idx or {}).get("window_days") or None
    newest_release = (idx or {}).get("release_tag") or None

    if idx is None or table is None:
        return {
            "due": True, "reason": "no committed FKS1 fill yet — the page is at day-0",
            "uncalled_new": uncalled_new, "total_in_window": total_in_window,
            "resolved_in_window": resolved_in_window, "coverage_pct": coverage_pct,
            "window_aged_out": None, "dark_date": None, "days_until_dark": None,
            "dark_lead_days": dark_lead_days, "committed_as_of": None,
            "committed_window_days": None, "newest_release": None, "threshold": threshold,
        }

    # Forecast, derived from as_of (the displayed window's start) so it stays consistent with
    # _window_aged_out: the page goes dark when as_of reaches the committed window's end. today is
    # as_of + window_days (the watchdog's own relation), so the calendar dark date resolves to the
    # canonical committed_as_of + 2*window.
    dark = projected_dark_date(committed_as_of, committed_window_days)
    start = _parse_iso_date(as_of)
    dark_date = dark.isoformat() if dark else None
    if dark and start is not None:
        today_eff = start + datetime.timedelta(days=window_days)
        days_until_dark = (dark - today_eff).days
    else:
        days_until_dark = None

    aged_out = _window_aged_out(committed_as_of, committed_window_days, as_of)
    if aged_out:
        reason = (f"last fill ({newest_release}, through {committed_as_of}) has fully rolled "
                  f"off the {window_days}d window — the page is showing no real calls")
        due = True
    elif uncalled_new >= threshold:
        reason = (f"{uncalled_new} new in-window isolate(s) uncalled "
                  f"(>= threshold {threshold}) since {newest_release}")
        due = True
    elif days_until_dark is not None and days_until_dark <= dark_lead_days:
        reason = (f"page goes dark {dark_date} ({days_until_dark}d) — within the {dark_lead_days}d "
                  f"lead time; run the billable refill before then to keep {newest_release} live")
        due = True
    else:
        reason = (f"only {uncalled_new} new in-window isolate(s) uncalled "
                  f"(< threshold {threshold}); {newest_release} still covers the window "
                  f"(dark {dark_date}, {days_until_dark}d out)")
        due = False

    return {
        "due": due, "reason": reason,
        "uncalled_new": uncalled_new, "total_in_window": total_in_window,
        "resolved_in_window": resolved_in_window, "coverage_pct": coverage_pct,
        "window_aged_out": aged_out, "dark_date": dark_date,
        "days_until_dark": days_until_dark, "dark_lead_days": dark_lead_days,
        "committed_as_of": committed_as_of,
        "committed_window_days": committed_window_days, "newest_release": newest_release,
        "threshold": threshold,
    }
