"""Cross-artifact agreement for the resistance-weather surface.

The page (openafr.weather.render_page), the JSON feed (weather_feed.build_data ->
data.json) and the Atom feed (weather_feed.render_atom -> alerts.atom) are three VIEWS over
one composed alert result. Each is already unit-tested in isolation
(test_weather.py / test_weather_feed.py); this module locks the contract *between* them:
given ONE result and ONE `now`, the three can never disagree on the alert facts or the
freshness framing, across the three honest states the product renders.

  * day-0 / WATCHING : no data -> all three calm, no alert, no badge, no false all-clear,
  * QUIET / stale fill: a fill aged out of the window -> all three carry the SAME
    "last batch covered through <date>" provenance and still show zero alerts,
  * fresh HIT        : real alerts -> every mutation + tier + Rung-1/2/3 position context
    appears in the page and the JSON record, and every alert has exactly one Atom entry
    titled "<PRIORITY>: <mutation>"; the counts match.

No network: reuses the shared #22 / FKS1 fixtures like test_weather.py.
"""
import json
import xml.etree.ElementTree as ET

from openafr import alert, weather, weather_feed

from test_weather import (  # noqa: E402  (shared fixtures / constants)
    NOW, _cov_result, _empty_result, _fks1_hit_result)

ATOM = "{http://www.w3.org/2005/Atom}"
_PROV = {"as_of": "2026-04-11", "n_called": 300}


def _render_all(result, *, now=NOW, **kw):
    """Render all three artifacts from the SAME result + now + kwargs (one source of truth)."""
    page = weather.render_page(result, now=now, **kw)
    data = weather_feed.build_data(result, now=now, **{
        k: v for k, v in kw.items() if k not in ("watching_since", "feeds")})
    root = ET.fromstring(weather_feed.render_atom(result, now=now,
                                                  page_url=kw.get("page_url")))
    return page, data, root


# ---- day-0 / WATCHING --------------------------------------------------------

def test_day0_all_three_are_calm_and_empty():
    page, data, root = _render_all(None)
    # page: watching copy, no tier badge
    assert "Watching" in page
    assert "No echinocandin-resistance variant is over threshold" in page
    assert 'class="badge"' not in page
    # json + atom: no alerts, no false all-clear artefact (certification stays null)
    assert data["alerts"] == [] and data["certification"] is None
    assert root.findall(f"{ATOM}entry") == []


def test_day0_erg11_arm_keeps_azole_copy_and_still_empty():
    r = _empty_result()
    r["gene"] = "ERG11"
    page, data, root = _render_all(r)
    assert "No azole-resistance variant is over threshold" in page
    assert data["gene"] == "ERG11" and data["alerts"] == []
    assert root.findall(f"{ATOM}entry") == []


# ---- QUIET / stale fill ------------------------------------------------------

def test_stale_fill_carries_one_provenance_date_and_no_alerts():
    # recent_total > 0 but recent_called == 0: the fill has aged out of the window. All three
    # must agree on the same "last batch covered through <date>" and still claim nothing.
    result = _cov_result(200, 0)
    page, data, root = _render_all(result, calls_provenance=_PROV)
    assert "awaiting re-call" in page
    assert f"the last batch covered through {_PROV['as_of']}" in page
    assert "window isolates re-called" not in page  # never a confident coverage %
    # the feed exposes the identical provenance date, and agrees there are no alerts
    assert data["calls_provenance"]["as_of"] == _PROV["as_of"]
    assert data["alerts"] == []
    assert root.findall(f"{ATOM}entry") == []


# ---- fresh HIT ---------------------------------------------------------------

def _rung3_ctx():
    sv = {
        "evidence_class": "proximity_label",
        "proximity_tier": "drug_contact",
        "min_dist_to_drug_A": 2.9,
        "label": "at drug-contact residue",
        "note": "Residue 635 sits within direct contact distance of the enfumafungin site.",
    }
    return {
        "tier": "novel_substitution_at_panel_residue",
        "window": "HS1", "anchor": "S639", "anchor_distance": -4,
        "at_panel_residue": True, "panel_letters_at_position": ["C", "Y"],
        "proximity": {"proximity_tier": "drug_contact", "min_dist_to_drug_A": 2.9},
        "structural_verdict": sv,
    }, sv


def test_hit_mutations_tiers_and_counts_agree_across_all_three():
    result = _fks1_hit_result()
    assert result["alerts"], "FKS1 fixture must produce alerts"
    page, data, root = _render_all(result)

    muts = [a["mutation"] for a in result["alerts"]]
    # json: same mutations, same order
    assert [a["mutation"] for a in data["alerts"]] == muts
    # atom: one entry per alert, titled "<PRIORITY>: <mutation>"
    entries = root.findall(f"{ATOM}entry")
    assert len(entries) == len(muts)
    titles = [e.find(f"{ATOM}title").text for e in entries]
    for a in result["alerts"]:
        badge = alert._BADGE[a["priority"]]
        assert a["mutation"] in page and badge in page          # page carries mutation + tier
        rec = next(r for r in data["alerts"] if r["mutation"] == a["mutation"])
        assert rec["priority"] == a["priority"]                  # json carries the same tier
        assert f"{a['priority'].upper()}: {a['mutation']}" in titles  # atom carries both


def test_hit_rung3_structural_verdict_agrees_between_page_and_feed():
    # When an alert carries a Rung-3 structural_verdict, the page and the JSON record must
    # surface the SAME label (atom is title/summary only, so it is not asserted here).
    ctx, sv = _rung3_ctx()
    result = _fks1_hit_result()
    result["alerts"][0]["position_context"] = ctx
    result["alerts"][0]["position_note"] = "It sits at a panel residue in HS1."
    page, data, _ = _render_all(result)
    # page: the certifiable contact label + the "spatial co-location only" honesty caveat
    assert "Structural label" in page
    assert sv["label"] in page and "spatial co-location only" in page
    # json: the whole verdict dict passes through intact, evidence_class readable without HTML
    rec = data["alerts"][0]
    assert rec["position_context"]["structural_verdict"] == sv
    assert rec["position_context"]["structural_verdict"]["evidence_class"] == "proximity_label"
    # a proximity_label must never have escalated the tier to ACT NOW
    assert result["alerts"][0]["priority"] != "act_now"


# ---- determinism + safety (integration-level) --------------------------------

def test_all_three_are_deterministic_given_now():
    result = _fks1_hit_result()
    p1, d1, r1 = _render_all(result)
    p2, d2, r2 = _render_all(result)
    assert p1 == p2                                   # page byte-stable
    assert json.dumps(d1, sort_keys=True) == json.dumps(d2, sort_keys=True)
    assert ET.tostring(r1) == ET.tostring(r2)         # atom byte-stable


def test_page_controlled_text_cannot_break_any_artifact():
    result = _fks1_hit_result()
    result["alerts"][0]["mutation"] = "<b>&x"         # hostile, page-controlled
    page, data, root = _render_all(result)
    assert "<b>&x" not in page                         # page escapes it
    assert data["alerts"][0]["mutation"] == "<b>&x"    # json keeps the raw value (JSON-safe)
    json.dumps(data)                                   # must not raise
    assert root.find(f"{ATOM}entry") is not None       # atom stayed well-formed (parsed above)
