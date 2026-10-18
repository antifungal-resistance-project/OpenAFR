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
1. MUST-STAY-LIT. If a committed FKS1 call table exists (a fill has landed, so the page is
   supposed to show real data), the published `data.json` must be lit -- coverage reports
   re-called isolates (`recent_called > 0` or `baseline_called > 0`). If a table exists but
   the feed came out empty, the overlay/pull regressed: FAIL, do not deploy. When no table is
   committed (genuine day-0), the watching feed is expected and this check passes -- preserving
   fail-open for quiet weeks.
2. FEED SCHEMA. `data.json` conforms to the contract in docs/FEED_SPEC.md -- the exact
   `schema_version`, top-level keys, and per-alert / emergence keys that
   tests/test_weather_feed.py pins. A partner builds against these; a malformed feed never
   ships.

Exit status 0 = healthy (publish), 1 = regression (block). Stdlib + openafr only; reads only.
"""
import argparse
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from openafr import earlywarning as ew  # noqa: E402
from openafr import weather_feed  # noqa: E402

DEFAULT_CALLS_DIR = ROOT / "data" / "earlywarning" / "fks1_calls"

# The schema contract from docs/FEED_SPEC.md, mirroring tests/test_weather_feed.py so the gate
# and the unit test enforce exactly the same keys (changing either is a schema_version bump).
TOP_KEYS = {"schema_version", "generated", "notice", "gene", "page_url", "as_of", "window_days",
            "thresholds", "coverage", "alerts", "release_trend", "forecast", "geography",
            "vus_census", "calls_provenance", "certification"}
ALERT_KEYS = {"mutation", "priority", "headline", "first_seen", "regions", "n_carriers_total",
              "emergence"}
EMERGENCE_KEYS = {"first_appearance", "rising", "possible_backlog", "recent_count",
                  "recent_called", "baseline_count", "baseline_called", "recent_prop",
                  "baseline_prop", "delta"}


def _check(cond, label, failures):
    print(f"  [{'PASS' if cond else 'FAIL'}] {label}")
    if not cond:
        failures.append(label)
    return bool(cond)


def _feed_is_lit(data):
    """True if the feed shows a re-called window (coverage has resolved isolates). The empty
    `{}` coverage render_weather emits on the watching-page fallback reads as not lit."""
    cov = data.get("coverage") or {}
    return int(cov.get("recent_called") or 0) > 0 or int(cov.get("baseline_called") or 0) > 0


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

        table = _calls_table(calls_dir)
        if table:
            lit = _feed_is_lit(data) if parsed else False
            if allow_watching and not lit:
                print(f"  [SKIP] a fill is committed ({pathlib.Path(table).name}) but "
                      "--allow-watching lets the watching page through")
            else:
                _check(lit,
                       f"page is lit (a fill is committed: {pathlib.Path(table).name}) — "
                       "coverage reports re-called isolates, not a silent day-0 revert",
                       failures)
        else:
            print("  [SKIP] no committed FKS1 call table — genuine day-0, watching page is "
                  "healthy (fail-open preserved)")

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
