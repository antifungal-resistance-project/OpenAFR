"""Emergence detection for the C. auris azole early-warning track (issue #22).

Pipeline stage: resistance early-warning / detection method.

Why this exists
---------------
This is the core method NCBI does not provide. NCBI's Pathogen Detection feed is
a metadata catalogue (spike #20); it tells you an isolate *exists*, never that a
resistance mutation is *newly appearing* or *rising*. Given the versioned snapshot
store (#21), this module answers the question that matters for early warning:

  Compared with what we had seen before, which ERG11 (CYP51) azole-resistance
  substitutions are new, and which are climbing?

Two signals, defined honestly
-----------------------------
  * FIRST-APPEARANCE -- a substitution present among recently-visible isolates
    that was never seen in the baseline window.
  * RISING           -- a substitution whose share of *called* isolates in the
    recent window exceeds its baseline share by at least `min_delta`.

Both require a minimum absolute isolate count (`min_count`) so a single sequenced
isolate cannot raise an alarm.

What "called" means, and why the denominator is honest
------------------------------------------------------
The store leaves `erg11_call` empty until the ERG11 re-caller populates it, and
records *why* it is empty (`resistance_source`) rather than defaulting an
uncalled isolate to wild-type (#21). This detector inherits that discipline: an
isolate with an empty call is NOT counted as susceptible. Proportions are taken
over *called* isolates only, and the calling coverage is reported alongside every
verdict, so "6% carry Y132F" is never silently read off 2%-called data. On today's
real snapshots coverage is 0% (calls pending), so `detect_emergence` returns no
signals and says so -- an honest "cannot detect yet", not a false "all clear".

The backlog artifact, and the guard
------------------------------------
NCBI visibility (`target_creation_date`) is not sample age (`collection_date`). A
lab depositing a years-old backlog makes old isolates appear *now*: a spike in the
recent window that is really archival. Each flagged mutation is therefore checked
for backlog -- if most of its recent carriers were *collected* on or before the
baseline cut, the signal is annotated `possible_backlog` with the reason, rather
than being dropped or trusted blindly.

Stdlib only, operating on snapshot records (list of dicts) as produced by
openafr.earlywarning.read_snapshot / normalize.
"""
import collections
import re

# A canonical ERG11 point substitution looks like Y132F (wt / position / mut).
# We normalize to upper-case and split multi-substitution calls on the usual
# separators. Tokens that do not look like substitutions are kept verbatim so we
# never silently discard something the re-caller emitted -- but see parse_call.
_SUB_RE = re.compile(r"^[A-Z]\d+[A-Z*]$")
_SPLIT_RE = re.compile(r"[,;/|\s]+")


def parse_call(call):
    """Turn an erg11_call cell into a set of substitution tokens.

    "Y132F,K143R" -> {"Y132F", "K143R"}. Empty / whitespace -> empty set.
    Tokens are upper-cased and de-duplicated. A token that does not match the
    wt-pos-mut shape is still returned (upper-cased); the caller decides whether
    to trust it. We never invent a call and never drop a nonempty one silently.
    """
    call = (call or "").strip()
    if not call:
        return set()
    return {t.upper() for t in _SPLIT_RE.split(call) if t}


def is_substitution(token):
    """True if `token` looks like a canonical point substitution (Y132F)."""
    return bool(_SUB_RE.match(token))


def _year(date):
    """Leading 4-digit year from an ISO-ish or bare-year date, else None."""
    m = re.match(r"\s*(\d{4})", date or "")
    return int(m.group(1)) if m else None


def _iso_day(date):
    """First 10 chars if they are a valid-looking ISO date, else ''."""
    d = (date or "")[:10]
    return d if len(d) == 10 and d[:4].isdigit() and d[4] == "-" else ""


def split_by_creation(records, as_of, window_days=180):
    """Partition records into (baseline, recent) by NCBI visibility date.

    baseline = isolates created on or before `as_of`.
    recent   = isolates created after `as_of` and within `window_days` of it.
    Isolates with no usable creation date cannot be placed in time and are
    returned separately so the gap is visible, never bucketed by guess.
    Returns (baseline, recent, undated).
    """
    as_of = str(as_of)[:10]
    cutoff = _add_days(as_of, window_days)
    baseline, recent, undated = [], [], []
    for rec in records:
        created = _iso_day(rec.get("target_creation_date", ""))
        if not created:
            undated.append(rec)
        elif created <= as_of:
            baseline.append(rec)
        elif created <= cutoff:
            recent.append(rec)
        # created beyond the window is intentionally ignored: it is newer than the
        # recent window this call is scoped to.
    return baseline, recent, undated


def _add_days(iso_day, days):
    import datetime as _dt
    y, m, d = (int(x) for x in iso_day.split("-"))
    return (_dt.date(y, m, d) + _dt.timedelta(days=days)).isoformat()


def called(records, call_field="erg11_call"):
    """The subset of records that carry a nonempty call in `call_field`.

    `call_field` selects which re-caller's column to read -- `erg11_call` (azole, the
    default) or `fks1_call` (echinocandin, v2). The honesty is identical for both: an empty
    call is NOT counted, so the denominator is always *called* isolates for that gene."""
    return [r for r in records if (r.get(call_field) or "").strip()]


def mutation_carriers(records, call_field="erg11_call"):
    """Map each substitution token -> list of records carrying it (called only).

    Only iterates over called isolates, so the returned counts are over the
    called denominator, never over all isolates. `call_field` selects the gene's call
    column (erg11_call or fks1_call), so the same detector serves both tracks.
    """
    carriers = collections.defaultdict(list)
    for r in called(records, call_field):
        for tok in parse_call(r.get(call_field)):
            carriers[tok].append(r)
    return carriers


def assess_backlog(recent_carriers, backlog_frac=0.5, max_lag_years=2):
    """Is a recent-window mutation spike really an archival backlog?

    The tell is deposit latency, not sample age per se: a carrier is "archival"
    if it became visible in NCBI (`target_creation_date`) at least
    `max_lag_years` after it was *collected* (`collection_date`) -- an old sample
    freshly deposited, not real-time emergence. Normal sequencing turnaround
    (collected, then deposited within a year or so) is not archival. If at least
    `backlog_frac` of the datable carriers are archival, the spike is flagged.
    Carriers missing either date are excluded from the ratio but reported, since
    we cannot judge their lag.

    Returns (possible_backlog: bool, reason: str).
    """
    archival = datable = undated = 0
    for r in recent_carriers:
        cy = _year(r.get("collection_date", ""))
        vy = _year(r.get("target_creation_date", ""))
        if cy is None or vy is None:
            undated += 1
            continue
        datable += 1
        if vy - cy >= max_lag_years:
            archival += 1
    if datable == 0:
        return False, f"undatable ({undated} carriers lack a usable date pair)"
    frac = archival / datable
    if frac >= backlog_frac:
        return True, (f"{archival}/{datable} carriers deposited >={max_lag_years}y "
                      f"after collection (archival {frac:.0%}); likely a backlog")
    return False, (f"{archival}/{datable} carriers show a long deposit lag "
                   f"(archival {frac:.0%}); reads as genuinely recent")


def detect_emergence(records, as_of, window_days=180, min_count=3,
                     min_delta=0.05, backlog_frac=0.5, call_field="erg11_call"):
    """Flag substitutions that are newly-appearing or rising, over one gene's calls.

    Splits `records` into a baseline (created <= as_of) and a recent window
    (as_of < created <= as_of + window_days) by NCBI visibility date, then, over
    *called* isolates only, compares each substitution's recent share with its
    baseline share.

    A substitution is reported if it clears `min_count` recent carriers AND is
    either FIRST-APPEARANCE (unseen in baseline) or RISING (recent share exceeds
    baseline share by >= `min_delta`). Each report carries a backlog verdict.

    Returns a dict:
      as_of, window_days, thresholds     -- echoed so a verdict is self-describing
      coverage                            -- baseline/recent totals and called counts
      signals                             -- list of flagged substitutions, each:
        {mutation, baseline_count, recent_count, baseline_called, recent_called,
         baseline_prop, recent_prop, delta, first_appearance, rising,
         possible_backlog, backlog_reason, carrier_keys}
      note                                -- set when nothing is callable yet
    """
    baseline, recent, undated = split_by_creation(records, as_of, window_days)
    b_called, r_called = called(baseline, call_field), called(recent, call_field)
    nb, nr = len(b_called), len(r_called)

    coverage = {
        "baseline_total": len(baseline), "baseline_called": nb,
        "recent_total": len(recent), "recent_called": nr,
        "undated": len(undated),
    }
    thresholds = {"min_count": min_count, "min_delta": min_delta,
                  "backlog_frac": backlog_frac}
    result = {"as_of": str(as_of)[:10], "window_days": window_days,
              "thresholds": thresholds, "coverage": coverage, "signals": []}

    result["call_field"] = call_field
    if nr == 0:
        gene = "FKS1" if call_field == "fks1_call" else "ERG11"
        result["note"] = (
            "no called isolates in the recent window -- cannot detect mutation "
            f"emergence yet. Calls are pending the {gene} re-caller; an empty call "
            "is honestly uncalled, never assumed susceptible.")
        return result

    base_carriers = mutation_carriers(baseline, call_field)
    rec_carriers = mutation_carriers(recent, call_field)

    signals = []
    for mut, carriers in rec_carriers.items():
        recent_count = len(carriers)
        if recent_count < min_count:
            continue
        baseline_count = len(base_carriers.get(mut, ()))
        baseline_prop = baseline_count / nb if nb else 0.0
        recent_prop = recent_count / nr
        delta = recent_prop - baseline_prop
        first_appearance = baseline_count == 0
        rising = delta >= min_delta
        if not (first_appearance or rising):
            continue
        possible_backlog, reason = assess_backlog(carriers, backlog_frac)
        signals.append({
            "mutation": mut,
            "canonical": is_substitution(mut),
            "baseline_count": baseline_count,
            "recent_count": recent_count,
            "baseline_called": nb,
            "recent_called": nr,
            "baseline_prop": baseline_prop,
            "recent_prop": recent_prop,
            "delta": delta,
            "first_appearance": first_appearance,
            "rising": rising,
            "possible_backlog": possible_backlog,
            "backlog_reason": reason,
            "carrier_keys": sorted(r["isolate_key"] for r in carriers),
        })

    # Loudest first: genuine (non-backlog) before backlog, then by recent share.
    signals.sort(key=lambda s: (s["possible_backlog"], -s["delta"],
                                -s["recent_count"], s["mutation"]))
    result["signals"] = signals
    return result
