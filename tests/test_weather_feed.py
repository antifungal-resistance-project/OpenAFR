"""Tests for the machine-readable weather feeds (openafr/weather_feed.py): data.json + Atom.

Same inputs as the page, so these lock that the feeds carry the page's facts, stay
deterministic, are valid XML/JSON, and escape page-controlled text.
"""
import json
import xml.etree.ElementTree as ET

from openafr import weather, weather_feed

from test_weather import NOW, _empty_result, _fks1_hit_result  # noqa: E402

ATOM = "{http://www.w3.org/2005/Atom}"
URL = "https://example.org/weather/"


def test_data_matches_alerts_and_is_json_serialisable():
    res = _fks1_hit_result()
    data = weather_feed.build_data(res, now=NOW, geo_breakdown={"USA": 3, "Kenya": 1},
                                   page_url=URL)
    assert data["schema_version"] == 1 and data["generated"] == "2026-09-17T08:30:00Z"
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
            "thresholds", "coverage", "alerts", "release_trend", "geography", "vus_census",
            "calls_provenance", "certification"}
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


def test_release_trend_carries_panel_vs_vus_split_unchanged():
    # #192: n_resistant is the panel-gated count and n_vus is reported separately; the feed
    # must pass both through untouched so a consumer never has to re-derive "resistant".
    trend = [{"release_tag": "PDG000000067.710", "as_of": "2026-04-14", "n_called": 3375,
              "n_resistant": 128, "n_vus": 69, "prevalence": 128 / 3375}]
    data = weather_feed.build_data(_fks1_hit_result(), now=NOW, release_trend=trend)
    assert data["release_trend"] == trend
    assert data["release_trend"][0]["prevalence"] == 128 / 3375
    json.dumps(data)


def test_vus_census_passes_through_with_position_context():
    # The per-variant breakdown behind n_vus: token + carriers + geometry-free position context,
    # passed through untouched so page and feed can't disagree. Absent -> empty list, not null.
    assert weather_feed.build_data(_fks1_hit_result(), now=NOW)["vus_census"] == []
    census = [{"token": "F635I", "carriers": 6,
               "position_context": {"tier": "novel_substitution_at_panel_residue"},
               "position_note": "It sits in hot-spot HS1, at the F635 anchor — ..."}]
    data = weather_feed.build_data(_fks1_hit_result(), now=NOW, vus_census=census)
    assert data["vus_census"] == census
    json.dumps(data)
