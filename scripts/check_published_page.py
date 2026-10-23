"""Publish gate for the resistance-weather page: fail-closed on a silent regression (#170).

Why this exists
---------------
`scripts/render_weather.py` is deliberately *fail-open*: on any data hiccup (a failed live
NCBI pull, or a committed FKS1 call table that stops matching by isolate_key) it swallows the
error and renders the honest day-0 "watching" page instead of crashing. That is the right call
for a genuinely quiet week -- a page that 500s on a quiet week is worse than one that says
"quiet week".

But the publish workflow (`.github/workflows/weather.yml`) renders and then *immediately
deploys*, so the same fail-open silently ships a REGRESSION once the page has been lit: the
real 3375-isolate echinocandin picture quietly reverts to "watching since ..." day-0 and
nobody is told. This script is the gate that turns that silent fail-open into a loud
fail-closed AT PUBLISH TIME, run between "render" and `deploy-pages`:

  - render still always produces a valid file (unchanged);
  - this gate refuses to let an *unexpected* regression ship.

What it checks (reads the already-rendered `site/` artifact -- it renders nothing)
----------------------------------------------------------------------------------
1. MUST-STAY-LIT (and WELL-COVERED). If a committed FKS1 call table exists (a fill has landed,
   so the page is supposed to show real data), the published `data.json` must be lit -- coverage
   reports a re-called fraction at or above COVERAGE_FLOOR, not merely > 0. This catches both a
   total day-0 revert AND a *partial* overlay collapse (e.g. isolate_key match breaks and only a
   handful overlay, which the old `recent_called > 0` check waved through). A lit page must also
   carry a non-empty geography breakdown. When no table is committed (genuine day-0), the watching
   feed is expected and these checks pass -- preserving fail-open for quiet weeks.
2. FEED SCHEMA. `data.json` conforms to the contract in docs/FEED_SPEC.md -- the exact
   `schema_version`, top-level keys, and per-alert / emergence keys that
   tests/test_weather_feed.py pins. A partner builds against these; a malformed feed never ships.
3. CONTENT SANITY. Nested producers keep their pinned keys -- `release_trend[]` entries carry the
   #192 panel-vs-VUS split (`n_resistant` AND `n_vus`, never silently re-conflated) and the
   `forecast` bundle carries its documented keys; every prevalence is a 0..1 fraction with no
   implausible release-to-release swing (MAX_PREVALENCE_SWING) that would mark a bad fill.
4. ATOM FEED. `alerts.atom` is present and well-formed XML (a subscribed consumer depends on it).

Exit status 0 = healthy (publish), 1 = regression (block). Stdlib + openafr only; reads only.
"""
import argparse
import json
import pathlib
import sys
import xml.etree.ElementTree as ET

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from openafr import earlywarning as ew  # noqa: E402
from openafr import weather_feed  # noqa: E402

DEFAULT_CALLS_DIR = ROOT / "data" / "earlywarning" / "fks1_calls"

# The schema contract from docs/FEED_SPEC.md, mirroring tests/test_weather_feed.py so the gate
# and the unit test enforce exactly the same keys (changing either is a schema_version bump).
TOP_KEYS = {"schema_version", "generated", "notice", "gene", "page_url", "as_of", "window_days",
            "fill_age_days", "fill_stale",
            "thresholds", "coverage", "alerts", "release_trend", "forecast", "geography",
            "geography_history", "vus_census", "calls_provenance", "certification"}
ALERT_KEYS = {"mutation", "priority", "headline", "first_seen", "regions", "n_carriers_total",
              "emergence"}
EMERGENCE_KEYS = {"first_appearance", "rising", "possible_backlog", "recent_count",
                  "recent_called", "baseline_count", "baseline_called", "recent_prop",
                  "baseline_prop", "delta"}
# Nested producers a consumer indexes into (scripts/render_weather.py _release_trend / _forecast).
# Pinning these catches a silent producer drift the top-level key check cannot see -- e.g. the
# #192 panel-vs-VUS split (n_resistant vs n_vus) silently reverting and re-conflating "resistant".
RELEASE_TREND_KEYS = {"release_tag", "as_of", "n_called", "n_resistant", "n_vus", "prevalence"}
FORECAST_KEYS = {"series", "skill", "n_folds", "model_mae", "naive_last_mae", "naive_pooled_mae",
                 "projection"}

# A lit page re-calls most of its displayed window; a partial-overlay collapse (e.g. isolate_key
# match breaks and only a handful overlay) leaves coverage > 0 but far below this floor. 0.5 sits
# well under the real page's ~100% anchored-window coverage yet far above a collapse.
COVERAGE_FLOOR = 0.5
# A plausible release-to-release move in panel-gated prevalence; a larger jump signals a bad fill
# (miscount, denominator slip) rather than real epidemiology. Absolute, on the 0..1 prevalence.
MAX_PREVALENCE_SWING = 0.25


def _check(cond, label, failures):
    print(f"  [{'PASS' if cond else 'FAIL'}] {label}")
    if not cond:
        failures.append(label)
    return bool(cond)


def _coverage_ratio(data):
    """Best re-called fraction across the recent/baseline axes (whichever carries isolates), or
    None when the coverage dict is empty/day-0. A *partial* overlay collapse shows up here as a
    low-but-nonzero ratio -- which the old `recent_called > 0` liveness check could not see."""
    cov = data.get("coverage") or {}
    ratios = []
    for total_k, called_k in (("recent_total", "recent_called"),
                              ("baseline_total", "baseline_called")):
        total = int(cov.get(total_k) or 0)
        if total > 0:
            ratios.append(int(cov.get(called_k) or 0) / total)
    return max(ratios) if ratios else None


def _feed_is_lit(data):
    """True if the feed shows a well-covered re-called window. A nonzero-but-collapsed overlay
    (ratio below COVERAGE_FLOOR) reads as NOT lit -- it is a regression, not a healthy page."""
    ratio = _coverage_ratio(data)
    return ratio is not None and ratio >= COVERAGE_FLOOR


def _validate_content(data, failures):
    """Content sanity beyond key-shape: a lit page must carry geography, nested producers must
    keep their pinned keys, and every prevalence must be a sane fraction with no implausible
    release-to-release swing. Catches bad-fill regressions the schema-key check sails past."""
    # Nested producer keys (release_trend entries + the forecast bundle).
    for i, e in enumerate(data.get("release_trend") or []):
        if not _check(isinstance(e, dict) and RELEASE_TREND_KEYS <= set(e),
                      f"release_trend[{i}] has all pinned keys", failures):
            if isinstance(e, dict):
                print(f"        missing: {sorted(RELEASE_TREND_KEYS - set(e))}")
    fc = data.get("forecast")
    if fc is not None:
        if not _check(isinstance(fc, dict) and FORECAST_KEYS <= set(fc),
                      "forecast carries all pinned keys", failures):
            if isinstance(fc, dict):
                print(f"        missing: {sorted(FORECAST_KEYS - set(fc))}")
    # Prevalence sanity: every release's prevalence is a 0..1 fraction, and consecutive releases
    # do not jump more than MAX_PREVALENCE_SWING (a bad fill, not real epidemiology).
    prevs = []
    for i, e in enumerate(data.get("release_trend") or []):
        p = (e or {}).get("prevalence")
        if p is None:
            continue
        _check(isinstance(p, (int, float)) and 0.0 <= p <= 1.0,
               f"release_trend[{i}] prevalence is a 0..1 fraction", failures)
        prevs.append(p)
    for a, b in zip(prevs, prevs[1:]):
        _check(abs(b - a) <= MAX_PREVALENCE_SWING,
               f"prevalence swing {a:.3f}->{b:.3f} within {MAX_PREVALENCE_SWING}", failures)


def _validate_schema(data, failures):
    """Enforce the docs/FEED_SPEC.md contract on an already-parsed data.json."""
    _check(data.get("schema_version") == weather_feed.SCHEMA_VERSION,
           f"schema_version == {weather_feed.SCHEMA_VERSION}", failures)
    keys = set(data)
    ok_top = _check(keys == TOP_KEYS, "data.json top-level keys match the pinned schema",
                    failures)
    if not ok_top:
        missing, extra = TOP_KEYS - keys, keys - TOP_KEYS
        if missing:
            print(f"        missing: {sorted(missing)}")
        if extra:
            print(f"        unexpected: {sorted(extra)}")
    # types for the fields a consumer indexes into
    _check(isinstance(data.get("alerts"), list), "alerts is a list", failures)
    _check(isinstance(data.get("coverage"), dict), "coverage is an object", failures)
    _check(isinstance(data.get("release_trend"), list), "release_trend is a list", failures)
    _check(isinstance(data.get("geography"), dict), "geography is an object", failures)
    gh = data.get("geography_history")
    if _check(isinstance(gh, list), "geography_history is a list", failures):
        _check(all(isinstance(e, dict) and isinstance(e.get("geography"), dict) for e in gh),
               "geography_history entries carry a geography object", failures)
    _check(isinstance(data.get("vus_census"), list), "vus_census is a list", failures)
    for i, a in enumerate(data.get("alerts") or []):
        if not _check(ALERT_KEYS <= set(a), f"alert[{i}] has all required keys", failures):
            print(f"        missing: {sorted(ALERT_KEYS - set(a))}")
        em = a.get("emergence")
        if not _check(isinstance(em, dict) and set(em) == EMERGENCE_KEYS,
                      f"alert[{i}].emergence keys match the schema", failures):
            if isinstance(em, dict):
                print(f"        got: {sorted(em)}")


def _calls_table(calls_dir):
    """The newest committed FKS1 call table under `calls_dir`, or None (genuine day-0)."""
    return ew.latest_fks1_call_table(str(calls_dir))


def run(site_dir, calls_dir, allow_watching=False):
    failures = []
    site = pathlib.Path(site_dir)

    _check((site / "index.html").is_file(), "index.html was written", failures)

    data_path = site / "data.json"
    if not _check(data_path.is_file(), "data.json was written", failures):
        print("\nFAIL — no data.json to validate; render did not complete.")
        return 1
    try:
        data = json.loads(data_path.read_text())
        parsed = _check(True, "data.json is valid JSON", failures)
    except ValueError as exc:
        _check(False, f"data.json is valid JSON ({exc})", failures)
        parsed = False

    if parsed:
        _validate_schema(data, failures)
        _validate_content(data, failures)

        table = _calls_table(calls_dir)
        if table:
            lit = _feed_is_lit(data) if parsed else False
            if allow_watching and not lit:
                print(f"  [SKIP] a fill is committed ({pathlib.Path(table).name}) but "
                      "--allow-watching lets the watching page through")
            else:
                if _check(lit,
                          f"page is lit (a fill is committed: {pathlib.Path(table).name}) — "
                          f"coverage ≥ {COVERAGE_FLOOR:.0%} re-called, not a day-0 revert "
                          "or a partial-overlay collapse",
                          failures):
                    # A lit page always has geography from the fill; an empty map is a regression.
                    _check(bool(data.get("geography")),
                           "lit page carries a non-empty geography breakdown", failures)
        else:
            print("  [SKIP] no committed FKS1 call table — genuine day-0, watching page is "
                  "healthy (fail-open preserved)")

    # The Atom feed is published alongside data.json; a consumer subscribes to it, so it must at
    # least be present and well-formed XML (render writes it from the same composed inputs).
    atom_path = site / "alerts.atom"
    if _check(atom_path.is_file(), "alerts.atom was written", failures):
        try:
            ET.fromstring(atom_path.read_text())
            _check(True, "alerts.atom is well-formed XML", failures)
        except ET.ParseError as exc:
            _check(False, f"alerts.atom is well-formed XML ({exc})", failures)

    print()
    if failures:
        print(f"FAIL — {len(failures)} check(s) failed: {'; '.join(failures)}")
        print("Refusing to publish: the previously-good page stays live.")
        return 1
    print("PASS — rendered artifact is healthy; safe to publish.")
    return 0


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site", default="site",
                    help="rendered artifact dir holding index.html + data.json (default: site)")
    ap.add_argument("--calls-dir", default=str(DEFAULT_CALLS_DIR),
                    help="committed FKS1 calls dir; a table here means the page must be lit")
    ap.add_argument("--allow-watching", action="store_true",
                    help="permit the watching page even when a fill is committed "
                         "(an intentional day-0 republish); off by default")
    args = ap.parse_args()
    return run(args.site, args.calls_dir, allow_watching=args.allow_watching)


if __name__ == "__main__":
    sys.exit(main())
