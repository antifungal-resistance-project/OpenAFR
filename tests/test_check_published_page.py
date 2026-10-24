"""Guard the publish gate (scripts/check_published_page.py).

The gate turns render_weather's deliberate fail-open into a loud fail-closed at PUBLISH time:
once a fill is committed, a page that silently reverted to day-0 "watching" (a broken pull or
overlay) must not deploy, and data.json must match its pinned schema. These tests lock the
three paths that matter:

  * lit feed + committed table            -> pass (publish)
  * watching feed + committed table       -> FAIL (the silent regression this closes)
  * watching feed + no committed table    -> pass (genuine quiet week, fail-open preserved)

plus the schema checks (version bump / dropped key / malformed alert). End-to-end fixtures
drive the REAL render_weather CLI through the prestage harness so they can't drift from
production; the schema cases build data.json directly.
"""
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
import check_published_page as gate  # noqa: E402
import prestage_weather_dryrun as dry  # noqa: E402
from openafr import weather_feed  # noqa: E402

AS_OF = "2026-01-01"
WINDOW = 180
RELEASE = 710


def _commit_table(calls_dir, n=50, resistant_frac=0.06, release=RELEASE):
    """Write a committed FKS1 call table + INDEX via the producer's own serializers."""
    calls_dir.mkdir(parents=True, exist_ok=True)
    tsv = calls_dir.parent / "pull.tsv"
    table = calls_dir / f"{dry.ew.ORG_PREFIX}.{release}.tsv"
    raw_rows, call_records = dry.synthesize(n, resistant_frac, AS_OF, WINDOW)
    dry._write_raw_tsv(tsv, raw_rows)
    sha, n_called = dry.ew.write_fks1_call_table(str(table), call_records)
    dry.ew.update_fks1_calls_index(str(calls_dir / "INDEX.tsv"), {
        "release_tag": f"{dry.ew.ORG_PREFIX}.{release}", "as_of": AS_OF,
        "window_days": WINDOW, "n_called": n_called, "sha256": sha,
        "written_at_utc": "2026-01-02T00:00:00Z",
    })
    return tsv, table


def _render_site(tmp_path, *, lit, n=50):
    """Render a real site/ dir. `lit=True` overlays a committed table (page lights up);
    `lit=False` renders --watching-only (the day-0 fallback)."""
    calls_dir = tmp_path / "fks1_calls"
    site = tmp_path / "site"
    site.mkdir()
    tsv, table = _commit_table(calls_dir, n=n)
    if lit:
        dry._run_render(tsv, table, AS_OF, WINDOW, site / "index.html")
    else:
        # The watching-page fallback, but with the table still committed: exactly the shape
        # of a broken live pull after a fill has landed.
        dry.subprocess.run(
            [sys.executable, str(dry.ROOT / "scripts" / "render_weather.py"),
             "--gene", "fks1", "--watching-only", "--out", str(site / "index.html")],
            capture_output=True, text=True, check=True)
    return site, calls_dir


# ---- the three live-loop paths -----------------------------------------------------------

def test_lit_feed_with_committed_table_passes(tmp_path):
    site, calls_dir = _render_site(tmp_path, lit=True)
    assert gate.run(str(site), str(calls_dir)) == 0


def test_watching_feed_with_committed_table_fails(tmp_path):
    # THE regression this plan closes: a fill is committed but the page came out day-0.
    site, calls_dir = _render_site(tmp_path, lit=False)
    assert gate.run(str(site), str(calls_dir)) == 1


def test_watching_feed_with_no_committed_table_passes(tmp_path):
    # Genuine quiet week: no fill committed, so the watching page is healthy (fail-open kept).
    site, _ = _render_site(tmp_path, lit=False)
    empty_calls = tmp_path / "empty_calls"
    empty_calls.mkdir()
    assert gate.run(str(site), str(empty_calls)) == 0


def test_allow_watching_overrides_the_lit_requirement(tmp_path):
    # An intentional day-0 republish while a fill is committed is permitted only on the flag.
    site, calls_dir = _render_site(tmp_path, lit=False)
    assert gate.run(str(site), str(calls_dir), allow_watching=True) == 0


# ---- schema contract ---------------------------------------------------------------------

def _write_site(tmp_path, data):
    site = tmp_path / "site"
    site.mkdir()
    (site / "index.html").write_text("<html></html>")
    (site / "data.json").write_text(json.dumps(data))
    (site / "alerts.atom").write_text(
        '<?xml version="1.0" encoding="utf-8"?>\n<feed xmlns="http://www.w3.org/2005/Atom"/>\n')
    return site


def _valid_data():
    # A schema-valid, lit data.json with one well-formed alert.
    return {k: v for k, v in {
        "schema_version": weather_feed.SCHEMA_VERSION,
        "generated": "2026-01-01T00:00:00Z",
        "notice": "Research use only.",
        "gene": "FKS1", "page_url": None, "as_of": AS_OF, "window_days": WINDOW,
        "fill_age_days": 20, "fill_stale": False,
        "thresholds": {},
        "coverage": {"recent_total": 10, "recent_called": 10,
                     "baseline_total": 20, "baseline_called": 20},
        "alerts": [{
            "mutation": "S639F", "priority": "watch", "headline": "", "first_seen": None,
            "regions": [], "n_carriers_total": 3,
            "emergence": {k: None for k in gate.EMERGENCE_KEYS},
        }],
        "release_trend": [], "forecast": None, "geography": {}, "geography_history": [],
        "panel_census": [], "geography_prevalence": None,
        "vus_census": [], "calls_provenance": {}, "certification": None,
    }.items()}


def test_valid_schema_passes_without_a_table(tmp_path):
    site = _write_site(tmp_path, _valid_data())
    empty = tmp_path / "empty"
    empty.mkdir()
    assert gate.run(str(site), str(empty)) == 0


def test_bumped_schema_version_fails(tmp_path):
    data = _valid_data()
    data["schema_version"] = weather_feed.SCHEMA_VERSION + 1
    site = _write_site(tmp_path, data)
    empty = tmp_path / "empty"
    empty.mkdir()
    assert gate.run(str(site), str(empty)) == 1


def test_malformed_geography_history_fails(tmp_path):
    # an entry missing its geography object must not slip through the gate
    data = _valid_data()
    data["geography_history"] = [{"release_tag": "PDG000000067.710", "as_of": "2026-04-14"}]
    site = _write_site(tmp_path, data)
    empty = tmp_path / "empty"
    empty.mkdir()
    assert gate.run(str(site), str(empty)) == 1


def test_dropped_top_level_key_fails(tmp_path):
    data = _valid_data()
    del data["coverage"]
    site = _write_site(tmp_path, data)
    empty = tmp_path / "empty"
    empty.mkdir()
    assert gate.run(str(site), str(empty)) == 1


def test_alert_missing_emergence_field_fails(tmp_path):
    data = _valid_data()
    data["alerts"][0]["emergence"].pop("delta")
    site = _write_site(tmp_path, data)
    empty = tmp_path / "empty"
    empty.mkdir()
    assert gate.run(str(site), str(empty)) == 1


# ---- content sanity beyond key-shape (strengthened gate) ---------------------------------

def _lit_site_with_table(tmp_path, data):
    """A site/ carrying arbitrary data.json PLUS a committed table, so the lit/coverage/geography
    content checks (which only run when a fill is committed) are exercised."""
    calls_dir = tmp_path / "fks1_calls"
    _commit_table(calls_dir)
    return _write_site(tmp_path, data), calls_dir


def test_partial_overlay_collapse_fails_even_though_coverage_is_nonzero(tmp_path):
    # 10 of 3375 overlaid: recent_called > 0 (old check passed) but ratio far below the floor.
    data = _valid_data()
    data["coverage"] = {"recent_total": 3375, "recent_called": 10,
                        "baseline_total": 0, "baseline_called": 0}
    data["geography"] = {"USA": 5}
    site, calls_dir = _lit_site_with_table(tmp_path, data)
    assert gate.run(str(site), str(calls_dir)) == 1


def test_lit_page_with_empty_geography_fails(tmp_path):
    data = _valid_data()  # coverage is 100% -> lit, but geography is {}
    site, calls_dir = _lit_site_with_table(tmp_path, data)
    assert gate.run(str(site), str(calls_dir)) == 1


def test_release_trend_missing_panel_vs_vus_key_fails(tmp_path):
    # #192 split must not silently revert: an entry without n_vus is a regression.
    data = _valid_data()
    data["release_trend"] = [{"release_tag": "PDG000000067.710", "as_of": "2026-04-14",
                              "n_called": 3375, "n_resistant": 128, "prevalence": 0.038}]
    site = _write_site(tmp_path, data)
    empty = tmp_path / "empty"
    empty.mkdir()
    assert gate.run(str(site), str(empty)) == 1


def test_implausible_prevalence_swing_fails(tmp_path):
    data = _valid_data()
    data["release_trend"] = [
        {"release_tag": "A", "as_of": "2026-01-01", "n_called": 100, "n_resistant": 4,
         "n_vus": 2, "prevalence": 0.04},
        {"release_tag": "B", "as_of": "2026-04-14", "n_called": 100, "n_resistant": 60,
         "n_vus": 2, "prevalence": 0.60},  # +0.56 jump: a bad fill, not epidemiology
    ]
    site = _write_site(tmp_path, data)
    empty = tmp_path / "empty"
    empty.mkdir()
    assert gate.run(str(site), str(empty)) == 1


def test_out_of_range_prevalence_fails(tmp_path):
    data = _valid_data()
    data["release_trend"] = [{"release_tag": "A", "as_of": "2026-01-01", "n_called": 100,
                              "n_resistant": 4, "n_vus": 2, "prevalence": 1.5}]
    site = _write_site(tmp_path, data)
    empty = tmp_path / "empty"
    empty.mkdir()
    assert gate.run(str(site), str(empty)) == 1


def test_malformed_atom_fails(tmp_path):
    site = _write_site(tmp_path, _valid_data())
    (site / "alerts.atom").write_text("<feed><unclosed>")
    empty = tmp_path / "empty"
    empty.mkdir()
    assert gate.run(str(site), str(empty)) == 1


def test_missing_atom_fails(tmp_path):
    site = _write_site(tmp_path, _valid_data())
    (site / "alerts.atom").unlink()
    empty = tmp_path / "empty"
    empty.mkdir()
    assert gate.run(str(site), str(empty)) == 1


def test_malformed_json_fails(tmp_path):
    site = tmp_path / "site"
    site.mkdir()
    (site / "index.html").write_text("<html></html>")
    (site / "data.json").write_text("{not valid json")
    empty = tmp_path / "empty"
    empty.mkdir()
    assert gate.run(str(site), str(empty)) == 1


def test_missing_data_json_fails(tmp_path):
    site = tmp_path / "site"
    site.mkdir()
    (site / "index.html").write_text("<html></html>")
    empty = tmp_path / "empty"
    empty.mkdir()
    assert gate.run(str(site), str(empty)) == 1
