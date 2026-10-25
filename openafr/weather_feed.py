"""Machine-readable feeds for the resistance weather page (data.json + alerts.atom).

The HTML page is for humans; a surveillance team wants to pull the same facts into their own
tooling. These are pure functions over exactly the inputs `weather.render_page` already takes
-- no I/O, no network, no wall-clock unless `now` is None -- so the feeds can never disagree
with the page. Detection-only and research-use, like the page.

  build_data(...)   -> JSON-serialisable dict (the `data.json` document)
  render_atom(...)  -> Atom 1.0 XML string (the `alerts.atom` feed; one entry per alert)
"""
import datetime
import hashlib
from xml.sax.saxutils import escape

from openafr.weather import SCOPE_LIMITS, _fill_age_tier

# v2 (2026-10): added `fill_age_days` + `fill_stale` so a consumer keying freshness off
# `generated` (refreshed on every weekly publish) is not misled by an anchored, frozen fill.
# v3 (2026-10): added `panel_census` -- the per-marker breakdown of the known-resistance panel
# (which substitution each resistant carrier has), whose carrier counts sum to the headline.
# v4 (2026-10): added `geography_prevalence` -- panel-gated prevalence (resistant / resolved,
# with a Wilson CI) per place plus the national headline rate, so a consumer can see WHERE
# resistance is proportionally elevated, not just where raw carriers are.
# v5 (2026-10): added `scope_limits` -- the claims this service deliberately does NOT make yet
# (chiefly a per-isolate calibrated probability) with the reason + the unlock condition, from the
# shared openafr.weather.SCOPE_LIMITS policy, so a partner sees what is withheld, not a silent gap.
SCHEMA_VERSION = 5
_TITLE = "C. auris Resistance Weather Report"
_NOTICE = "Research use only -- a triage aid, not a clinical claim. Detection-only."


def _now_iso(now):
    if now is None:
        now = datetime.datetime.now(datetime.timezone.utc)
    if isinstance(now, str):
        now = datetime.datetime.fromisoformat(now.replace("Z", "+00:00"))
    return now.astimezone(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _alert_record(a):
    e = a["emergence"]
    rec = {
        "mutation": a["mutation"],
        "priority": a["priority"],
        "headline": a.get("headline") or "",
        "first_seen": a.get("first_seen"),
        "regions": list(a.get("regions") or []),
        "n_carriers_total": a.get("n_carriers_total", 0),
        "emergence": {k: e.get(k) for k in (
            "first_appearance", "rising", "possible_backlog", "recent_count",
            "recent_called", "baseline_count", "baseline_called", "recent_prop",
            "baseline_prop", "delta")},
    }
    if a.get("position_context"):
        rec["position_context"] = a["position_context"]
        rec["position_note"] = a.get("position_note") or ""
    return rec


def build_data(alert_result, *, now=None, calls_provenance=None, certification=None,
               release_trend=None, geo_breakdown=None, geo_history=None, geo_prevalence=None,
               panel_census=None, vus_census=None, forecast=None, page_url=None):
    """The `data.json` document: status, alerts, coverage, trend, forecast, geography,
    panel breakdown, geography_prevalence, VUS, provenance. `geography_prevalence` is the
    per-place panel-gated prevalence bundle (places with Wilson CIs + the national headline
    rate), or None. `forecast` is the openafr.forecast bundle (monthly series
    + walk-forward skill + one-step projection); a consumer reads `forecast.skill` before
    trusting `projection`. `panel_census` is the per-marker breakdown of the known-resistance
    panel, each `{token, carriers, low_ppv}`, whose carrier counts sum to the headline.
    `scope_limits` is the fixed policy list of claims the service deliberately does NOT make yet
    (each `{claim, status, reason, unlock, source}`), from openafr.weather.SCOPE_LIMITS so the
    feed and the page state the same deferrals."""
    r = alert_result or {}
    now_date = datetime.datetime.fromisoformat(
        _now_iso(now).replace("Z", "+00:00")).date()
    fill_tier, fill_age_days = _fill_age_tier(
        calls_provenance, r.get("window_days"), now_date)
    return {
        "schema_version": SCHEMA_VERSION,
        "generated": _now_iso(now),
        "notice": _NOTICE,
        "gene": r.get("gene"),
        "page_url": page_url,
        "as_of": r.get("as_of"),
        "window_days": r.get("window_days"),
        # Age of the overlaid fill at build time, and whether it has aged past its sampling
        # window -- the freshness signal `generated` alone cannot carry (see weather._age_clause).
        "fill_age_days": fill_age_days,
        "fill_stale": fill_tier in ("stale", "dark"),
        "thresholds": r.get("thresholds") or {},
        "coverage": r.get("coverage") or {},
        "alerts": [_alert_record(a) for a in (r.get("alerts") or [])],
        "release_trend": list(release_trend or []),
        "forecast": dict(forecast) if forecast else None,
        "geography": dict(sorted((geo_breakdown or {}).items())),
        "geography_history": [
            {"release_tag": e.get("release_tag"), "as_of": e.get("as_of"),
             "geography": dict(sorted((e.get("geography") or {}).items()))}
            for e in (geo_history or [])],
        "panel_census": list(panel_census or []),
        "geography_prevalence": dict(geo_prevalence) if geo_prevalence else None,
        "vus_census": list(vus_census or []),
        "scope_limits": [dict(lim) for lim in SCOPE_LIMITS],
        "calls_provenance": dict(calls_provenance or {}),
        "certification": dict(certification) if certification else None,
    }


def _entry_id(a, as_of):
    """Stable per (mutation, tier) so a reader marks an alert read once, but a tier change
    (e.g. WATCH -> ACT NOW) surfaces as a fresh entry."""
    digest = hashlib.sha256(f"{a['mutation']}|{a['priority']}".encode()).hexdigest()[:16]
    return f"urn:openafr:weather:fks1:{digest}"


def render_atom(alert_result, *, now=None, page_url=None):
    """Atom 1.0 feed of the current alerts. An empty feed (no alerts) is valid and honest."""
    r = alert_result or {}
    updated = _now_iso(now)
    alerts = r.get("alerts") or []
    link = escape(str(page_url)) if page_url else ""
    out = ['<?xml version="1.0" encoding="utf-8"?>',
           '<feed xmlns="http://www.w3.org/2005/Atom">',
           f"  <title>{escape(_TITLE)}</title>",
           f"  <subtitle>{escape(_NOTICE)}</subtitle>",
           "  <id>urn:openafr:weather:fks1</id>",
           f"  <updated>{updated}</updated>"]
    if link:
        out.append(f'  <link rel="alternate" href="{link}"/>')
    for a in alerts:
        e = a["emergence"]
        summary = (f"{a.get('headline') or ''} Recent {e['recent_count']}/"
                   f"{e['recent_called']}, regions: "
                   f"{', '.join(a.get('regions') or []) or 'undisclosed'}.").strip()
        out += ["  <entry>",
                f"    <title>{escape(a['priority'].upper())}: {escape(str(a['mutation']))}</title>",
                f"    <id>{_entry_id(a, r.get('as_of'))}</id>",
                f"    <updated>{updated}</updated>",
                f"    <category term=\"{escape(a['priority'])}\"/>"]
        if link:
            out.append(f'    <link rel="alternate" href="{link}"/>')
        out += [f"    <summary>{escape(summary)}</summary>", "  </entry>"]
    out.append("</feed>")
    return "\n".join(out) + "\n"
