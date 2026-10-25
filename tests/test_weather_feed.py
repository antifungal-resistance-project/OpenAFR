"""Tests for the machine-readable weather feeds (openafr/weather_feed.py): data.json + Atom.

Same inputs as the page, so these lock that the feeds carry the page's facts, stay
deterministic, are valid XML/JSON, and escape page-controlled text.
"""
import json
import pathlib
import xml.etree.ElementTree as ET

from openafr import weather, weather_feed

from test_weather import NOW, _empty_result, _fks1_hit_result  # noqa: E402

ATOM = "{http://www.w3.org/2005/Atom}"
URL = "https://example.org/weather/"


def test_data_matches_alerts_and_is_json_serialisable():
    res = _fks1_hit_result()
    data = weather_feed.build_data(res, now=NOW, geo_breakdown={"USA": 3, "Kenya": 1},
                                   page_url=URL)
    assert data["schema_version"] == 5 and data["generated"] == "2026-09-17T08:30:00Z"
    assert [a["mutation"] for a in data["alerts"]] == [a["mutation"] for a in res["alerts"]]
    assert list(data["geography"]) == ["Kenya", "USA"]  # sorted, stable
    json.dumps(data)  # must not raise


def test_empty_result_is_calm_not_broken():
    data = weather_feed.build_data(None, now=NOW)
    assert data["alerts"] == [] and data["certification"] is None
    root = ET.fromstring(weather_feed.render_atom(_empty_result(), now=NOW))
    assert root.findall(f"{ATOM}entry") == []


def test_atom_is_valid_one_entry_per_alert_with_stable_ids():
    res = _fks1_hit_result()
    xml = weather_feed.render_atom(res, now=NOW, page_url=URL)
    root = ET.fromstring(xml)
    entries = root.findall(f"{ATOM}entry")
    assert len(entries) == len(res["alerts"]) > 0
    ids = [e.find(f"{ATOM}id").text for e in entries]
    assert len(set(ids)) == len(ids)
    assert ids == [e.find(f"{ATOM}id").text for e in ET.fromstring(
        weather_feed.render_atom(res, now="2027-01-01T00:00:00Z")).findall(f"{ATOM}entry")]


def test_atom_escapes_text():
    res = _fks1_hit_result()
    res["alerts"][0]["mutation"] = "<b>&x"
    ET.fromstring(weather_feed.render_atom(res, now=NOW))  # still well-formed


def test_deterministic():
    res = _fks1_hit_result()
    assert weather_feed.render_atom(res, now=NOW) == weather_feed.render_atom(res, now=NOW)


def test_page_links_feeds_only_when_asked():
    assert "alerts.atom" not in weather.render_page(None, now=NOW)
    page = weather.render_page(None, now=NOW, feeds=True)
    assert 'href="data.json"' in page and 'rel="alternate" type="application/atom+xml"' in page


# ---- schema contract (docs/FEED_SPEC.md) -------------------------------------------------
# A partner builds against these keys. Changing them is a schema_version bump, not a refactor.

TOP_KEYS = {"schema_version", "generated", "notice", "gene", "page_url", "as_of", "window_days",
            "fill_age_days", "fill_stale",
            "thresholds", "coverage", "alerts", "release_trend", "forecast", "geography",
            "geography_history", "panel_census", "geography_prevalence", "vus_census",
            "scope_limits", "calls_provenance", "certification"}
ALERT_KEYS = {"mutation", "priority", "headline", "first_seen", "regions", "n_carriers_total",
              "emergence"}
EMERGENCE_KEYS = {"first_appearance", "rising", "possible_backlog", "recent_count",
                  "recent_called", "baseline_count", "baseline_called", "recent_prop",
                  "baseline_prop", "delta"}


def test_data_json_schema_is_pinned():
    data = weather_feed.build_data(_fks1_hit_result(), now=NOW)
    assert set(data) == TOP_KEYS
    for a in data["alerts"]:
        assert ALERT_KEYS <= set(a)
        assert set(a["emergence"]) == EMERGENCE_KEYS


def test_scope_limits_are_published_and_name_the_deferred_probability():
    # The feed states what the service deliberately does NOT claim yet (chiefly a calibrated
    # per-isolate probability), from the same policy the page renders, so a partner sees the
    # deferral, not a silent gap.
    data = weather_feed.build_data(_fks1_hit_result(), now=NOW)
    sl = data["scope_limits"]
    assert isinstance(sl, list) and len(sl) > 0
    first = sl[0]
    assert {"claim", "status", "reason"} <= set(first)
    assert "probabilit" in first["claim"].lower()      # calibrated-probability deferral leads
    assert sl == [dict(lim) for lim in weather.SCOPE_LIMITS]  # feed == shared policy


def test_feed_spec_doc_tracks_the_emitted_schema():
    # docs/FEED_SPEC.md is the partner-facing contract, but nothing asserted it tracks the code,
    # so it silently drifted to schema 3 while the feed emitted 4 and geography_prevalence went
    # undocumented. A partner building against the spec would never learn the key exists.
    # Regression: FEED_SPEC-drift — found by /qa on 2026-10-25.
    # Report: .gstack/qa-reports/qa-report-weather-2026-10-25.md
    spec = (pathlib.Path(__file__).resolve().parents[1] / "docs" / "FEED_SPEC.md").read_text()
    assert f"currently `{weather_feed.SCHEMA_VERSION}`" in spec, \
        "FEED_SPEC schema_version is stale vs weather_feed.SCHEMA_VERSION"
    for key in TOP_KEYS:
        assert f"`{key}`" in spec, f"FEED_SPEC does not document top-level key {key!r}"


def test_fill_age_and_staleness_reported_for_consumers():
    # A consumer keying freshness off `generated` (bumped every weekly publish) must also see the
    # overlaid fill's real age, and whether it has aged past its sampling window.
    res = _fks1_hit_result()
    res["window_days"] = 180
    fresh = weather_feed.build_data(res, now="2026-05-01T00:00:00Z",
                                    calls_provenance={"as_of": "2026-04-14"})
    assert fresh["fill_age_days"] == 17 and fresh["fill_stale"] is False
    stale = weather_feed.build_data(res, now="2026-10-22T00:00:00Z",
                                    calls_provenance={"as_of": "2026-04-14"})
    assert stale["fill_age_days"] == 191 and stale["fill_stale"] is True


def test_release_trend_carries_panel_vs_vus_split_unchanged():
    # #192: n_resistant is the panel-gated count and n_vus is reported separately; the feed
    # must pass both through untouched so a consumer never has to re-derive "resistant".
    trend = [{"release_tag": "PDG000000067.710", "as_of": "2026-04-14", "n_called": 3375,
              "n_resistant": 128, "n_vus": 69, "prevalence": 128 / 3375}]
    data = weather_feed.build_data(_fks1_hit_result(), now=NOW, release_trend=trend)
    assert data["release_trend"] == trend


def test_geography_history_passes_through_sorted_and_defaults_empty():
    # Absent by default (backward-compatible addition), and each release's geography is sorted
    # by country so the feed is byte-stable regardless of input order.
    assert weather_feed.build_data(_fks1_hit_result(), now=NOW)["geography_history"] == []
    hist = [{"release_tag": "PDG000000067.709", "as_of": "2026-04-11",
             "geography": {"USA": 8, "China": 1}}]
    data = weather_feed.build_data(_fks1_hit_result(), now=NOW, geo_history=hist)
    assert data["geography_history"] == [
        {"release_tag": "PDG000000067.709", "as_of": "2026-04-11",
         "geography": {"China": 1, "USA": 8}}]


def test_forecast_passes_through_and_defaults_null():
    # Absent by default (backward-compatible addition), and passed through verbatim so a
    # consumer reads `forecast.skill` before trusting `forecast.projection`.
    assert weather_feed.build_data(_fks1_hit_result(), now=NOW)["forecast"] is None
    fc = {"series": [{"month": "2026-04", "n_resolved": 369, "n_resistant": 19,
                      "prevalence": 19 / 369}],
          "skill": False, "projection": {"naive_pooled": 0.038}}
    data = weather_feed.build_data(_fks1_hit_result(), now=NOW, forecast=fc)
    assert data["forecast"] == fc
    json.dumps(data)  # must not raise


def test_vus_census_passes_through_with_position_context():
    # The per-variant breakdown behind n_vus: token + carriers + position context (Rungs 1–3),
    # passed through untouched so page and feed can't disagree. Absent -> empty list, not null.
    assert weather_feed.build_data(_fks1_hit_result(), now=NOW)["vus_census"] == []
    census = [{"token": "F635I", "carriers": 6,
               "position_context": {"tier": "novel_substitution_at_panel_residue"},
               "position_note": "It sits in hot-spot HS1, at the F635 anchor — ..."}]
    data = weather_feed.build_data(_fks1_hit_result(), now=NOW, vus_census=census)
    assert data["vus_census"] == census
    json.dumps(data)


def test_panel_census_passes_through():
    # The per-marker breakdown behind n_resistant: token + carriers + low_ppv flag, passed
    # through untouched so page and feed can't disagree. Absent -> empty list, not null.
    assert weather_feed.build_data(_fks1_hit_result(), now=NOW)["panel_census"] == []
    census = [{"token": "S639F", "carriers": 20, "low_ppv": False},
              {"token": "D642Y", "carriers": 50, "low_ppv": True}]
    data = weather_feed.build_data(_fks1_hit_result(), now=NOW, panel_census=census)
    assert data["panel_census"] == census
    json.dumps(data)


def test_geography_prevalence_passes_through():
    # Per-place prevalence bundle: null when absent (backward-compatible), passed through
    # untouched otherwise so page and feed can't disagree on WHERE resistance is elevated.
    assert weather_feed.build_data(_fks1_hit_result(), now=NOW)["geography_prevalence"] is None
    geo_prev = {
        "headline_prevalence": 0.1,
        "places": [{"place": "USA", "n_resolved": 40, "n_panel_positive": 12,
                    "event_frequency": 0.3, "ci95": (0.18, 0.45)}],
        "pooled_small": {"n_places": 1, "n_resolved": 1, "n_panel_positive": 0},
        "n_unplaced": 2,
    }
    data = weather_feed.build_data(_fks1_hit_result(), now=NOW, geo_prevalence=geo_prev)
    assert data["geography_prevalence"]["places"][0]["place"] == "USA"
    assert data["geography_prevalence"]["headline_prevalence"] == 0.1
    json.dumps(data)


def test_alert_position_context_with_structural_verdict_passes_through():
    # Rung 3: when an alert carries a position_context dict that includes a structural_verdict
    # sub-field (proximity_label), _alert_record must pass the whole dict through intact so the
    # feed consumer can read evidence_class, label, and note without parsing the HTML page.
    sv = {
        "evidence_class": "proximity_label",
        "proximity_tier": "drug_contact",
        "min_dist_to_drug_A": 3.2,
        "label": "at drug-contact residue",
        "note": "Residue 637 in HS1 sits within direct contact distance of the enfumafungin binding site.",
    }
    ctx = {
        "tier": "in_hotspot_adjacent_to_panel",
        "window": "HS1",
        "anchor": "S639",
        "anchor_distance": -2,
        "at_panel_residue": False,
        "panel_letters_at_position": [],
        "nearest_core_marker": "639",
        "distance_to_nearest": 2,
        "proximity": {"proximity_tier": "drug_contact", "min_dist_to_drug_A": 3.2},
        "structural_verdict": sv,
    }
    res = _fks1_hit_result()
    res["alerts"][0]["position_context"] = ctx
    res["alerts"][0]["position_note"] = "It sits adjacent to S639 in HS1."
    data = weather_feed.build_data(res, now=NOW)
    rec = data["alerts"][0]
    assert rec["position_context"]["structural_verdict"] == sv
    assert rec["position_context"]["structural_verdict"]["evidence_class"] == "proximity_label"
    json.dumps(data)  # must not raise
