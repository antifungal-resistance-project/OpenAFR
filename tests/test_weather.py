"""Tests for the resistance weather page (openafr/weather.py).

The page is a VIEW over alert.compose_alerts output, so these lock the three states it
must render honestly and the guarantee that it never invents or leaks:

  * WATCHING / day-0: no alerts -> a calm, explicit "nothing over threshold" page, never
    a badge and never a false all-clear,
  * HIT: real fixture alerts -> one card each, carrying the mutation, its tier badge, and
    its structural fit verdict pulled straight from the composed dict,
  * QUIET framing: a prior watch-start date -> a "N day(s)" line, not the day-0 line,
  * determinism: given `now`, the page is byte-stable (no wall-clock leak),
  * safety: text from the result is HTML-escaped, and the pocket asset path is
    filename-safe so a compound token can't escape the pockets/ dir.
No network: the hit case runs against the shared #22 fixture like test_alert.py.
"""
import pathlib
import sys

from openafr import alert, weather

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
from detect_emergence import _demo_records  # noqa: E402  (the shared #22 fixture)
from compose_alert import _demo_fks1_records  # noqa: E402  (the FKS1 analog)

AS_OF = "2024-12-31"
NOW = "2026-09-17T08:30:00Z"


def _hit_result():
    return alert.compose_alerts(_demo_records(), AS_OF, window_days=180,
                                min_count=3, min_delta=0.05)


def _fks1_hit_result():
    return alert.compose_fks1_alerts(_demo_fks1_records(), AS_OF, window_days=180,
                                     min_count=3, min_delta=0.05)


def _empty_result(note=None):
    r = {"as_of": AS_OF, "window_days": 180, "thresholds": {},
         "coverage": {}, "alerts": []}
    if note is not None:
        r["note"] = note
    return r


# ---- WATCHING / day-0 --------------------------------------------------------

def test_none_result_renders_watching_page():
    # Default arm is the echinocandin (FKS1) product; a bare None result shows its copy.
    page = weather.render_page(None, now=NOW)
    assert "WATCHING" in page
    assert "No echinocandin-resistance variant is over threshold" in page
    assert 'class="badge"' not in page  # no tier badge element when there are no alerts


def test_erg11_gene_result_keeps_azole_calm_copy():
    # The research arm (gene=ERG11) still renders the azole wording.
    r = _empty_result()
    r["gene"] = "ERG11"
    page = weather.render_page(r, now=NOW)
    assert "No azole-resistance variant is over threshold" in page
    assert "ERG11 substitution" in page


def test_empty_result_shows_note_when_present():
    page = weather.render_page(_empty_result(note="no erg11 calls to assess"), now=NOW)
    assert "no erg11 calls to assess" in page
    assert "Last checked" in page


# ---- HIT ---------------------------------------------------------------------

def test_hit_page_renders_one_card_per_alert_with_verdict():
    result = _hit_result()
    assert result["alerts"], "fixture must produce alerts for this test to mean anything"
    page = weather.render_page(result, now=NOW)
    for a in result["alerts"]:
        assert a["mutation"] in page
        assert alert._BADGE[a["priority"]] in page
        # the structural fit verdict text is surfaced verbatim, not summarised away
        assert str(a["structural"]["fluconazole_fit_verdict"]) in page
    # a hit page references the pre-rendered pocket asset for its top alert
    assert weather.pocket_asset(result["alerts"][0]["mutation"]) in page


def test_hit_status_line_leads_with_resistance():
    page = weather.render_page(_hit_result(), now=NOW)
    assert "RESISTANCE SEEN" in page


def test_first_seen_renders_collected_with_visible_not_a_raw_dict():
    assert weather._fmt_first_seen(
        {"collected": "2023", "visible": "2024-03-01"}) == "2023 (visible 2024-03-01)"
    assert weather._fmt_first_seen({"collected": "2023", "visible": ""}) == "2023"
    assert weather._fmt_first_seen(
        {"collected": "", "visible": "2024-03-01"}) == "visible 2024-03-01"
    assert weather._fmt_first_seen({"collected": "", "visible": ""}) == "unknown"
    assert weather._fmt_first_seen(None) == "unknown"


def test_hit_page_never_shows_raw_first_seen_dict():
    # Regression: the "First seen" row once rendered str(dict) -> {'collected': ...}.
    page = weather.render_page(_hit_result(), now=NOW)
    assert "'collected'" not in page and "&#x27;collected&#x27;" not in page


# ---- FKS1 (echinocandin, detection-only) HIT --------------------------------

def test_fks1_hit_page_renders_detection_only_cards():
    result = _fks1_hit_result()
    assert result["alerts"], "FKS1 fixture must produce alerts"
    page = weather.render_page(result, now=NOW)
    assert "RESISTANCE SEEN" in page
    for a in result["alerts"]:
        assert a["mutation"] in page
    # detection-only: no structural fit row, and the caveat rides in on the headline
    assert "Structural fit" not in page
    assert "detection only" in page


def test_fks1_hit_page_has_no_pocket_image():
    # No FKS1 pocket render exists, so the card must not reference a pockets/ asset.
    page = weather.render_page(_fks1_hit_result(), now=NOW)
    assert "pockets/" not in page
    assert 'class="pocket"' not in page


# ---- QUIET framing -----------------------------------------------------------

def test_quiet_framing_counts_days_since_watch_start():
    page = weather.render_page(_empty_result(), now=NOW, watching_since="2026-09-01")
    assert "QUIET" in page
    assert "16 day(s)" in page  # 2026-09-01 -> 2026-09-17


def test_day0_line_when_no_watch_start():
    page = weather.render_page(_empty_result(), now=NOW)
    assert "WATCHING" in page
    assert "QUIET" not in page


# ---- determinism + safety ----------------------------------------------------

def test_page_is_deterministic_given_now():
    a = weather.render_page(_empty_result(), now=NOW, watching_since="2026-01-01")
    b = weather.render_page(_empty_result(), now=NOW, watching_since="2026-01-01")
    assert a == b
    assert "2026-09-17 08:30 UTC" in a


def test_now_defaults_to_wall_clock_when_omitted():
    page = weather.render_page(_empty_result())  # no `now` -> uses UTC now
    assert "Last checked" in page and "UTC" in page


def test_unparseable_watch_start_falls_back_to_watching_line():
    # A malformed watching_since must not raise; it degrades to the day-0 line.
    page = weather.render_page(_empty_result(), now=NOW, watching_since="not-a-date")
    assert "WATCHING" in page
    assert "day(s)" not in page


def test_result_text_is_html_escaped():
    page = weather.render_page(_empty_result(note="<script>alert(1)</script>"), now=NOW)
    assert "<script>alert(1)</script>" not in page
    assert "&lt;script&gt;" in page


def test_pocket_asset_is_filename_safe():
    assert weather.pocket_asset("F126L") == "pockets/F126L.png"
    assert weather.pocket_asset("TR34/L98H") == "pockets/TR34L98H.png"
    assert "/" not in weather.pocket_asset("../../etc/passwd").removeprefix("pockets/")


# ---- freshness / coverage honesty (issue #171) -------------------------------
# A public-health page must state how much of the shown window is actually re-called and
# how current those calls are -- it must degrade honestly as a fill ages out of the window,
# never present a stale/partial sample as a confident current picture.

def _cov_result(recent_total, recent_called):
    r = _empty_result()
    r["gene"] = "FKS1"
    r["coverage"] = {"recent_total": recent_total, "recent_called": recent_called,
                     "baseline_total": 0, "baseline_called": 0, "undated": 0}
    return r


def test_healthy_coverage_shows_percent_and_fill_date():
    page = weather.render_page(_cov_result(200, 150), now=NOW,
                               calls_provenance={"as_of": "2026-04-11", "n_called": 300})
    assert "150 of 200 window isolates re-called (75%)" in page
    assert "calls current through" in page and "2026-04-11" in page


def test_window_with_no_recalls_reads_as_awaiting_not_current():
    # recent_total > 0 but recent_called == 0: the fill has aged out of the window.
    page = weather.render_page(_cov_result(200, 0), now=NOW,
                               calls_provenance={"as_of": "2026-04-11", "n_called": 300})
    assert "awaiting re-call" in page
    assert "Awaiting the next echinocandin re-call batch" in page
    assert "the last batch covered through 2026-04-11" in page
    assert "window isolates re-called" not in page  # never a confident coverage %


def test_quiet_window_omits_coverage_clause():
    # recent_total == 0: a genuinely quiet window, not a stale fill -> no coverage clause,
    # the normal watching copy stands.
    page = weather.render_page(_cov_result(0, 0), now=NOW,
                               calls_provenance={"as_of": "2026-04-11", "n_called": 300})
    assert "re-called" not in page.split("<footer>")[0]  # body/meta carry no coverage %
    assert "awaiting re-call" not in page
    assert "No echinocandin-resistance variant is over threshold" in page


def test_missing_coverage_and_provenance_omit_clause_without_crashing():
    # The page-never-fails invariant: no coverage / no provenance -> clause simply absent.
    page = weather.render_page(_empty_result(), now=NOW, calls_provenance=None)
    assert "window isolates re-called" not in page
    assert "awaiting re-call" not in page


def test_footer_states_the_monthly_refresh_cadence():
    page = weather.render_page(_cov_result(200, 150), now=NOW,
                               calls_provenance={"as_of": "2026-04-11", "n_called": 300})
    assert "monthly FKS1 re-call batch" in page


# ---- certification badge + methods / provenance (track C) --------------------
# A public-health product must show that its numbers come from a certified caller, and
# expose its methods/provenance -- but only when there is a real PASS record to stand on.

_CERT_PASS = {"verdict": "PASS", "n": 98, "n_resolved": 98,
              "sensitivity": 1.0, "specificity": 1.0, "identity": 1.0,
              "benchmark": "fks1_benchmark_100.tsv"}


def test_certification_pass_renders_badge_and_methods():
    page = weather.render_page(
        _cov_result(200, 150), now=NOW,
        calls_provenance={"as_of": "2026-04-11", "n_called": 300,
                          "release_tag": "PDG000000067.709"},
        certification=_CERT_PASS)
    assert 'class="cert"' in page
    assert "certified FKS1 caller" in page
    assert 'id="methods"' in page
    # measured numbers are surfaced, and the release the fill targeted is named
    assert "100%" in page and "PASS (98/98 isolates resolved)" in page
    assert "PDG000000067.709" in page and "300" in page


def test_no_certification_omits_badge_and_methods():
    # The default (no certification) page makes no certification claim at all.
    page = weather.render_page(_cov_result(200, 150), now=NOW,
                               calls_provenance={"as_of": "2026-04-11", "n_called": 300})
    assert 'class="cert"' not in page
    assert 'id="methods"' not in page


def test_non_pass_verdict_makes_no_certification_claim():
    page = weather.render_page(_cov_result(200, 150), now=NOW,
                               certification={**_CERT_PASS, "verdict": "UNDERPOWERED"})
    assert 'class="cert"' not in page


def test_head_carries_shareable_meta():
    page = weather.render_page(_fks1_hit_result(), now=NOW,
                               page_url="https://example.org/weather/")
    for tag in ('name="description"', 'property="og:title"', 'property="og:description"',
                'name="twitter:card"'):
        assert tag in page
    assert 'property="og:url" content="https://example.org/weather/"' in page
    # the description reflects the real status, not boilerplate
    assert "RESISTANCE SEEN" in page.split("<body>")[0]


# ---- render-layer polish: design tokens, head completeness, meta hierarchy --------
# Presentation only -- no new data. The palette must live in exactly one place (:root
# tokens), the <head> must be complete (canonical, favicon, structured data), and the
# meta facts must keep their honest degradation under a cleaner hierarchy.

import json as _json  # noqa: E402


def test_palette_lives_in_root_tokens_not_inline_hex():
    page = weather.render_page(_fks1_hit_result(), now=NOW)
    # The tier palette is declared once as :root custom properties ...
    assert ":root" in page
    for token in ("--tier-act-now:", "--tier-watch:", "--tier-context:", "--accent:"):
        assert token in page
    # ... and consumed via var(), never re-injected as inline hex on the cards.
    assert 'style="--accent:#' not in page
    assert 'style="background:#' not in page


def test_alert_cards_name_their_tier_by_class():
    result = _fks1_hit_result()
    page = weather.render_page(result, now=NOW)
    for a in result["alerts"]:
        assert f'class="card tier-{a["priority"]}"' in page


def test_canonical_link_emitted_only_with_page_url():
    with_url = weather.render_page(_empty_result(), now=NOW,
                                   page_url="https://example.org/weather/")
    assert '<link rel="canonical" href="https://example.org/weather/">' in with_url
    without = weather.render_page(_empty_result(), now=NOW)
    assert 'rel="canonical"' not in without


def test_favicon_is_present_inline():
    page = weather.render_page(_empty_result(), now=NOW)
    assert 'rel="icon"' in page
    assert "data:image/svg+xml," in page


def _extract_jsonld(page):
    tag = '<script type="application/ld+json">'
    chunk = page.split(tag, 1)[1]
    return _json.loads(chunk.split("</script>", 1)[0])


def test_structured_data_is_valid_json_and_describes_the_dataset():
    page = weather.render_page(_empty_result(), now=NOW,
                               page_url="https://example.org/weather/")
    data = _extract_jsonld(page)  # must parse
    assert data["@type"] == "MedicalWebPage"
    assert data["url"] == "https://example.org/weather/"
    assert data["mainEntity"]["@type"] == "Dataset"


def test_structured_data_carries_cert_metrics_only_on_pass():
    passed = weather.render_page(_cov_result(200, 150), now=NOW, certification=_CERT_PASS)
    measured = _extract_jsonld(passed)["mainEntity"]["variableMeasured"]
    values = {m["name"]: m["value"] for m in measured}
    assert values["certified caller sensitivity"] == 100
    assert values["certified caller specificity"] == 100
    plain = weather.render_page(_cov_result(200, 150), now=NOW)
    assert "variableMeasured" not in _extract_jsonld(plain)["mainEntity"]


def test_meta_freshness_leads_on_its_own_line_when_covered():
    page = weather.render_page(_cov_result(200, 150), now=NOW,
                               calls_provenance={"as_of": "2026-04-11", "n_called": 300})
    assert 'class="meta-fresh"' in page
    assert 'class="meta-detail"' in page
    # the coverage fact and the provenance detail both survive the regroup
    assert "150 of 200 window isolates re-called (75%)" in page
    assert "Last checked" in page


def test_meta_omits_freshness_line_when_window_is_quiet():
    page = weather.render_page(_cov_result(0, 0), now=NOW,
                               calls_provenance={"as_of": "2026-04-11", "n_called": 300})
    assert 'class="meta-fresh"' not in page
    assert 'class="meta-detail"' in page


# ---- deepening: release trend + geographic breakdown -------------------------
# The page should read like weather, not a flat snapshot: a prevalence trajectory across
# committed re-call batches (reproducible) and where the current release's carriers are
# (current-release only, honestly labelled). Both degrade to nothing when their data is absent.

def _trend(*points):
    # points: (as_of, n_called, n_resistant) oldest-first; prevalence derived.
    return [{"as_of": a, "n_called": c, "n_resistant": r,
             "prevalence": (r / c) if c else None} for (a, c, r) in points]


def test_trend_section_shows_per_batch_prevalence_and_rising_delta():
    page = weather.render_page(
        _cov_result(200, 150), now=NOW,
        release_trend=_trend(("2026-04-11", 300, 15), ("2026-04-14", 3375, 197)))
    assert 'class="trend"' in page
    assert "Resistance trend" in page
    # both batches render with resolved counts and a one-decimal prevalence
    assert "15 of 300 resolved" in page and "5.0%" in page
    assert "197 of 3375 resolved" in page and "5.8%" in page
    # 5.0% -> 5.8% is a rise; the delta names the direction and magnitude
    assert "trend-up" in page and "prevalence up 0.8 pts" in page


def test_trend_section_marks_a_falling_batch():
    page = weather.render_page(
        _cov_result(200, 150), now=NOW,
        release_trend=_trend(("2026-04-11", 100, 10), ("2026-04-14", 100, 6)))
    assert "trend-down" in page and "prevalence down 4.0 pts" in page


def test_trend_section_omitted_for_a_single_batch():
    # One batch is a point, not a trend -- no section until there is a trajectory.
    page = weather.render_page(_cov_result(200, 150), now=NOW,
                               release_trend=_trend(("2026-04-14", 3375, 197)))
    assert 'class="trend"' not in page
    assert "Resistance trend" not in page


def test_geo_section_ranks_countries_and_labels_current_release():
    page = weather.render_page(
        _cov_result(200, 150), now=NOW,
        geo_breakdown={"India": 12, "USA": 5, "South Africa": 12})
    assert 'class="geo"' in page
    assert "current release" in page
    # ranked by count desc, ties broken by name -> India, South Africa (both 12), then USA (5)
    order = [page.index("India"), page.index("South Africa"), page.index("USA")]
    assert order == sorted(order)


def test_trend_and_geo_absent_by_default_without_crashing():
    page = weather.render_page(_cov_result(200, 150), now=NOW)
    assert 'class="trend"' not in page
    assert 'class="geo"' not in page


def _vus(*pairs):
    """Build an enriched VUS census (token, carriers) the way render_weather does, via the
    real geometry-free position-context tool so the fixture matches production wiring."""
    from openafr import fks1_position
    out = []
    for token, carriers in pairs:
        ctx = fks1_position.position_context(token)
        out.append({"token": token, "carriers": carriers, "position_context": ctx,
                    "position_note": f"It {fks1_position.describe(ctx)}"})
    return out


def test_vus_section_lists_each_variant_with_carriers_and_geography():
    page = weather.render_page(
        _cov_result(200, 150), now=NOW, vus_census=_vus(("F635I", 6), ("R1354H", 4)))
    assert 'class="vus"' in page
    assert "Variants of uncertain significance" in page
    assert "F635I" in page and "R1354H" in page
    assert "6 carrier(s)" in page and "4 carrier(s)" in page
    # the geography sentence rides in, and the tier label is shown
    assert "hot-spot" in page
    # detection-only boundary is restated: tracked, not counted as resistance
    assert "NOT counted as known resistance" in page


def test_vus_section_absent_without_census():
    page = weather.render_page(_cov_result(200, 150), now=NOW)
    assert 'class="vus"' not in page
    page_empty = weather.render_page(_cov_result(200, 150), now=NOW, vus_census=[])
    assert 'class="vus"' not in page_empty


def test_trend_section_surfaces_vus_separately_from_resistant():
    # #192: non-panel variants of uncertain significance are tracked, never headlined as resistant.
    trend = _trend(("2026-04-11", 300, 8), ("2026-04-14", 3375, 128))
    trend[-1]["n_vus"] = 69
    page = weather.render_page(_cov_result(200, 150), now=NOW, release_trend=trend)
    assert "128 of 3375 resolved" in page and "3.8%" in page
    assert "69 resolved isolate(s) carry a non-panel FKS1 variant of uncertain significance" in page
    assert "5.8%" not in page


def test_trend_section_omits_vus_note_when_none():
    page = weather.render_page(
        _cov_result(200, 150), now=NOW,
        release_trend=_trend(("2026-04-11", 300, 8), ("2026-04-14", 3375, 128)))
    assert "uncertain significance" not in page
