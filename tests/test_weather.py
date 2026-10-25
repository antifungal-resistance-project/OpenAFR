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
    assert "Watching" in page
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
    assert "Resistance detected" in page


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
    assert "Resistance detected" in page
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


def _fks1_contact_zone_records():
    """A snapshot where a real non-core HS3 variant (W691C) emerges over threshold.

    W691C sits in the enfumafungin drug-contact zone, so fks1_structural.structural_verdict
    returns a non-None contact-label dict for it. That dict has the FKS1 shape
    (evidence_class/proximity_tier/label/...), NOT the ERG11 fluconazole-fit shape the page's
    _alert_card reads from a `structural` field -- the collision this test guards against.
    """
    recs = []

    def add(key, created, collected, call):
        recs.append({
            "isolate_key": key, "target_acc": key + ".1",
            "target_creation_date": created, "collection_date": collected,
            "country": "Pakistan", "fks1_call": call,
            "run_acc": "SRR" + key[-3:], "fks1_resistance_source": "recalled",
        })

    for i in range(10):
        add(f"PDB{i:03d}", "2024-03-01", "2023", "")      # baseline, uncalled
    for i in range(4):
        add(f"PDR1{i:02d}", "2025-02-01", "2024", "W691C")  # genuine first-appearance
    return recs


def test_fks1_contact_zone_variant_renders_without_crashing():
    # Regression: a FKS1 alert for a contact-zone variant once put the FKS1 structural_verdict
    # dict on the `structural` channel the ERG11 renderer reads as a fluconazole fit -> KeyError
    # inside _alert_card -> failed the deploy (render_page is outside render_weather's fail-open
    # guard). The composer must not populate that channel on the FKS1 path.
    result = alert.compose_fks1_alerts(_fks1_contact_zone_records(), AS_OF,
                                       window_days=180, min_count=3, min_delta=0.05)
    assert any(a["mutation"] == "W691C" for a in result["alerts"]), \
        "fixture must flag the contact-zone variant"
    w = next(a for a in result["alerts"] if a["mutation"] == "W691C")
    assert "structural" not in w, "FKS1 path must not populate the ERG11 fit channel"
    page = weather.render_page(result, now=NOW)  # must not raise
    # detection-only discipline holds: no fit row, no pocket asset
    assert "Structural fit" not in page
    assert "pockets/" not in page and 'class="pocket"' not in page
    # the Rung-3 contact label still reaches the card via position_context
    assert "spatial co-location only" in page
    assert "at drug-contact residue" in page


def test_alert_card_shows_rung3_structural_verdict_when_present():
    # Rung 3: an alert card for a drug_contact variant must show the structural_verdict label
    # (matches the _vus_section rendering discipline). No label when no structural_verdict.
    sv = {
        "evidence_class": "proximity_label",
        "proximity_tier": "drug_contact",
        "min_dist_to_drug_A": 2.9,
        "label": "at drug-contact residue",
        "note": "Residue 635 sits within direct contact distance of the enfumafungin binding site.",
    }
    ctx = {
        "tier": "novel_substitution_at_panel_residue",
        "window": "HS1",
        "anchor": "S639",
        "anchor_distance": -4,
        "at_panel_residue": True,
        "panel_letters_at_position": ["C", "Y"],
        "proximity": {"proximity_tier": "drug_contact", "min_dist_to_drug_A": 2.9,
                      "reference_structure": "PDB 9WZU", "note": "direct contact"},
        "structural_verdict": sv,
    }
    result = _fks1_hit_result()
    result["alerts"][0]["position_context"] = ctx
    result["alerts"][0]["position_note"] = "It sits at a panel residue in HS1."
    page = weather.render_page(result, now=NOW)
    assert "Structural label" in page
    assert "at drug-contact residue" in page
    assert "spatial co-location only" in page
    # absence: an alert without structural_verdict must not show the label
    result2 = _fks1_hit_result()
    ctx2 = dict(ctx)
    del ctx2["structural_verdict"]
    result2["alerts"][0]["position_context"] = ctx2
    result2["alerts"][0]["position_note"] = "It sits at a panel residue in HS1."
    page2 = weather.render_page(result2, now=NOW)
    assert "Structural label" not in page2


# ---- QUIET framing -----------------------------------------------------------

def test_quiet_framing_counts_days_since_watch_start():
    page = weather.render_page(_empty_result(), now=NOW, watching_since="2026-09-01")
    assert "Quiet" in page
    assert "16 days" in page  # 2026-09-01 -> 2026-09-17


def test_day0_line_when_no_watch_start():
    page = weather.render_page(_empty_result(), now=NOW)
    assert "Watching" in page
    assert "Quiet" not in page


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
    assert "Watching" in page
    assert "with nothing new" not in page


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


# ---- fill-age honesty: the anchored page must not look fresher than its data ----------------
# The display window is anchored to the committed fill, so coverage reads a confident 100% with a
# static "calls current through <date>" forever. These lock the age signal that keeps the page from
# silently overstating freshness as the fill ages out from under a weekly-refreshed "Last checked".

def _aged_result(as_of, window_days):
    r = _cov_result(200, 200)
    r["window_days"] = window_days
    r["as_of"] = as_of
    return r, {"as_of": as_of, "n_called": 200}


def test_fresh_fill_shows_no_age_warning():
    # now within one window of the fill -> the "calls current through" date carries it; no warning.
    res, prov = _aged_result("2026-09-01", 180)
    page = weather.render_page(res, now="2026-10-01T00:00:00Z", calls_provenance=prov)
    assert "fill is" not in page
    assert 'class="stale"' not in page and 'class="dark"' not in page


def test_stale_fill_shows_amber_age_warning():
    # now past one window but under two -> visible amber "fill is N days old" warning.
    res, prov = _aged_result("2026-04-14", 180)
    page = weather.render_page(res, now="2026-10-22T00:00:00Z", calls_provenance=prov)
    assert "window isolates re-called" in page            # still lit/confident on coverage
    assert '<span class="stale">fill is 191 days old' in page
    assert "older than the sampling window" in page
    assert 'class="dark"' not in page


def test_aged_out_fill_shows_red_dark_warning():
    # now past two windows -> the fill's whole window has rolled off; red "refresh overdue".
    res, prov = _aged_result("2026-04-14", 180)
    page = weather.render_page(res, now="2027-05-01T00:00:00Z", calls_provenance=prov)
    assert '<span class="dark">fill is' in page
    assert "refresh overdue" in page


def test_age_warning_omitted_when_window_unknown():
    # fail-honest: no window_days -> cannot tier the age, so no warning rather than a guess.
    res = _cov_result(200, 200)
    res["as_of"] = "2026-04-14"
    res["window_days"] = None  # unknown window -> cannot tier the age
    page = weather.render_page(res, now="2027-05-01T00:00:00Z",
                               calls_provenance={"as_of": "2026-04-14", "n_called": 200})
    assert "fill is" not in page


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


def test_scope_limits_section_names_the_deferred_probability():
    # The "what we don't claim yet" block renders on every page (it is static product policy,
    # not gated on a fill or certification) and leads with the calibrated-probability deferral.
    page = weather.render_page(_fks1_hit_result(), now=NOW)
    assert 'id="scope-limits"' in page
    assert "does not claim yet" in page
    assert "calibrated resistance probability" in page


def test_head_carries_shareable_meta():
    page = weather.render_page(_fks1_hit_result(), now=NOW,
                               page_url="https://example.org/weather/")
    for tag in ('name="description"', 'property="og:title"', 'property="og:description"',
                'name="twitter:card"'):
        assert tag in page
    assert 'property="og:url" content="https://example.org/weather/"' in page
    # the description reflects the real status, not boilerplate
    assert "Resistance detected" in page.split("<body>")[0]


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


def test_meta_shows_fill_provenance_when_window_is_quiet():
    # When recent_total == 0 (NCBI release predates the window), the meta-fresh line
    # should still surface the fill provenance so the reader knows data exists.
    page = weather.render_page(_cov_result(0, 0), now=NOW,
                               calls_provenance={"as_of": "2026-04-11", "n_called": 300})
    assert 'class="meta-fresh"' in page
    assert "3,375" not in page or "300" in page  # uses provenance n_called
    assert "2026-04-11" in page
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
    assert "15 of 300 isolates" in page and "5.0%" in page
    assert "197 of 3375 isolates" in page and "5.8%" in page
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


def _fc_series(*points):
    # points: (month, n_resolved, n_resistant) oldest-first; prevalence derived.
    return [{"month": m, "n_resolved": nr, "n_resistant": np,
             "prevalence": (np / nr) if nr else None} for (m, nr, np) in points]


def test_forecast_section_shows_projection_only_when_skill_certified():
    fc = {"series": _fc_series(("2026-04", 200, 4), ("2026-05", 200, 8), ("2026-06", 200, 14)),
          "skill": True,
          "projection": {"next_month": "2026-07", "point": 0.09, "lo": 0.05, "hi": 0.15,
                         "naive_pooled": 0.043}}
    page = weather.render_page(_cov_result(200, 150), now=NOW, forecast=fc)
    assert 'class="trend forecast"' in page
    assert "How prevalence is trending" in page
    assert "4 of 200 isolates" in page and "14 of 200 isolates" in page
    assert "Projected" in page and "<code>2026-07</code>" in page and "9.0%" in page
    assert "95% PI 5.0%&ndash;15.0%" in page
    assert "trend-up" in page


def test_forecast_section_degrades_to_trajectory_when_no_skill():
    fc = {"series": _fc_series(("2026-04", 369, 19), ("2026-05", 481, 19), ("2026-06", 651, 18)),
          "skill": False,
          "projection": {"next_month": "2026-07", "point": 0.03, "lo": 0.004, "hi": 0.15,
                         "naive_pooled": 0.038}}
    page = weather.render_page(_cov_result(200, 150), now=NOW, forecast=fc)
    assert 'class="trend forecast"' in page
    assert "No forecast skill demonstrated" in page
    assert "roughly stationary near 3.8%" in page
    # no projected number is claimed on the no-skill path
    assert "Projected" not in page


def test_forecast_section_omitted_below_two_months():
    fc = {"series": _fc_series(("2026-06", 651, 18)), "skill": False, "projection": None}
    page = weather.render_page(_cov_result(200, 150), now=NOW, forecast=fc)
    assert "forecast" not in page
    assert "Prevalence trajectory" not in page


# ---- at-a-glance sparklines (page for readers) -------------------------------
# Each time-series section renders a small inline-SVG trajectory above its number table so a
# reader sees the direction at a glance. The chart is additive: the <dl> of exact values stays.

def _spark_svg(page):
    # The sparkline SVG substring, for colour-source assertions. "" when none rendered.
    i = page.find('<svg class="spark"')
    return page[i:page.find("</svg>", i) + 6] if i != -1 else ""


def test_trend_section_renders_a_sparkline_with_one_marker_per_batch():
    page = weather.render_page(
        _cov_result(200, 150), now=NOW,
        release_trend=_trend(("2026-04-11", 300, 15), ("2026-04-14", 3375, 197)))
    assert '<svg class="spark"' in page
    assert page.count('class="spark-dot"') == 2  # one marker per batch
    assert 'role="img"' in page and "aria-label=" in page
    # the number table is still there -- the chart is additive, not a replacement
    assert "15 of 300 isolates" in page


def test_forecast_sparkline_carries_projection_only_when_skill_certified():
    series = _fc_series(("2026-04", 200, 4), ("2026-05", 200, 8), ("2026-06", 200, 14))
    skilled = {"series": series, "skill": True,
               "projection": {"next_month": "2026-07", "point": 0.09, "lo": 0.05, "hi": 0.15,
                              "naive_pooled": 0.043}}
    page = weather.render_page(_cov_result(200, 150), now=NOW, forecast=skilled)
    assert '<svg class="spark"' in page
    assert 'class="spark-proj"' in page  # the projected-point marker
    assert 'class="spark-pi"' in page    # its 95% PI whisker

    noskill = dict(skilled, skill=False)
    page2 = weather.render_page(_cov_result(200, 150), now=NOW, forecast=noskill)
    assert '<svg class="spark"' in page2           # trajectory still charts
    assert 'class="spark-proj"' not in page2       # but no forecast marker
    assert 'class="spark-pi"' not in page2


def test_sparkline_omitted_below_two_points():
    # One batch is a point, not a trajectory -- no chart (the section itself is also omitted).
    page = weather.render_page(_cov_result(200, 150), now=NOW,
                               release_trend=_trend(("2026-04-14", 3375, 197)))
    assert '<svg class="spark"' not in page


def test_sparkline_colour_comes_from_tokens_not_hardcoded_hex():
    # Every colour on the chart must flow through the :root CSS tokens so light/dark/print work;
    # a stray hex in the SVG markup would escape the theme.
    page = weather.render_page(
        _cov_result(200, 150), now=NOW,
        release_trend=_trend(("2026-04-11", 300, 15), ("2026-04-14", 3375, 197)))
    svg = _spark_svg(page)
    assert svg and "#" not in svg


def test_geo_section_ranks_countries_and_labels_current_release():
    page = weather.render_page(
        _cov_result(200, 150), now=NOW,
        geo_breakdown={"India": 12, "USA": 5, "South Africa": 12})
    assert 'class="geo"' in page
    assert "latest batch" in page
    # ranked by count desc, ties broken by name -> India, South Africa (both 12), then USA (5)
    order = [page.index("India"), page.index("South Africa"), page.index("USA")]
    assert order == sorted(order)


def test_trend_and_geo_absent_by_default_without_crashing():
    page = weather.render_page(_cov_result(200, 150), now=NOW)
    assert 'class="trend"' not in page
    assert 'class="geo"' not in page
    assert 'class="geo geo-history"' not in page


def test_geo_history_renders_one_ranked_block_per_release():
    page = weather.render_page(
        _cov_result(200, 150), now=NOW,
        geo_history=[
            {"release_tag": "PDG000000067.709", "as_of": "2026-04-11",
             "geography": {"USA": 8, "China": 1}},
            {"release_tag": "PDG000000067.710", "as_of": "2026-04-14",
             "geography": {"USA": 111, "China": 14, "Colombia": 1}},
        ])
    assert 'class="geo geo-history"' in page
    assert "over time" in page
    # both releases appear, each labelled with its as-of date
    assert "PDG000000067.709 &middot; as of 2026-04-11" in page
    assert "PDG000000067.710 &middot; as of 2026-04-14" in page
    # within a release, countries rank by count desc -> USA (111) before China (14)
    block = page[page.index("PDG000000067.710"):]
    assert block.index("USA") < block.index("China") < block.index("Colombia")


def test_geo_history_absent_without_data_and_skips_empty_releases():
    assert 'class="geo geo-history"' not in weather.render_page(
        _cov_result(200, 150), now=NOW, geo_history=None)
    # a release whose table carried no country (empty breakdown) contributes nothing
    assert 'class="geo geo-history"' not in weather.render_page(
        _cov_result(200, 150), now=NOW,
        geo_history=[{"release_tag": "PDG000000067.709", "as_of": "2026-04-11",
                      "geography": {}}])


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
    assert "uncertain significance" in page
    assert "F635I" in page and "R1354H" in page
    assert "6 carriers" in page and "4 carriers" in page
    # the variant list stacks (header line + full-width note) via .vus-list, never the
    # two-column .grid: its long token+tags dt makes the max-content grid squeeze the note
    # into a cramped far-right column. Regression guard for that layout bug (/qa 2026-10-25).
    assert 'class="vus-list"' in page
    vus_block = page[page.index('class="vus"'):]
    assert 'class="grid"' not in vus_block[:vus_block.index("</section>")]
    # the geography sentence rides in, and the tier label is shown
    assert "hot-spot" in page
    # detection-only boundary is restated: tracked, not counted as resistance
    assert "not counted as resistant" in page


def test_vus_section_absent_without_census():
    page = weather.render_page(_cov_result(200, 150), now=NOW)
    assert 'class="vus"' not in page
    page_empty = weather.render_page(_cov_result(200, 150), now=NOW, vus_census=[])
    assert 'class="vus"' not in page_empty


def test_panel_breakdown_lists_each_marker_with_carriers():
    census = [{"token": "S639F", "carriers": 20, "low_ppv": False},
              {"token": "D642Y", "carriers": 50, "low_ppv": True}]
    page = weather.render_page(_cov_result(200, 150), now=NOW, panel_census=census)
    assert 'class="panel-breakdown"' in page
    assert "Which resistance variant is circulating" in page
    assert "S639F" in page and "D642Y" in page
    assert "20 carriers" in page and "50 carriers" in page
    # the low-PPV marker is tagged; the high-PPV one is not
    assert '<code>D642Y</code> <span class="low-ppv">weaker marker (low PPV)</span>' in page
    assert '<code>S639F</code></dt>' in page
    # the low-PPV caveat only appears because a low-PPV marker is present
    assert "split resistant/susceptible" in page


def test_panel_breakdown_omits_low_ppv_note_when_all_core():
    census = [{"token": "S639F", "carriers": 20, "low_ppv": False}]
    page = weather.render_page(_cov_result(200, 150), now=NOW, panel_census=census)
    assert 'class="panel-breakdown"' in page
    assert "low PPV" not in page
    assert "split resistant/susceptible" not in page


def test_panel_breakdown_absent_without_census():
    page = weather.render_page(_cov_result(200, 150), now=NOW)
    assert 'class="panel-breakdown"' not in page
    page_empty = weather.render_page(_cov_result(200, 150), now=NOW, panel_census=[])
    assert 'class="panel-breakdown"' not in page_empty


def test_geo_prevalence_ranks_places_and_flags_above_national():
    geo_prev = {
        "headline_prevalence": 0.10,
        "places": [
            {"place": "USA: California", "n_resolved": 40, "n_panel_positive": 12,
             "event_frequency": 0.30, "ci95": (0.18, 0.45)},   # CI low 18% > 10% -> flagged
            {"place": "United Kingdom", "n_resolved": 60, "n_panel_positive": 6,
             "event_frequency": 0.10, "ci95": (0.05, 0.20)},   # CI low 5% < 10% -> not flagged
        ],
        "pooled_small": {"n_places": 2, "n_resolved": 3, "n_panel_positive": 1},
        "n_unplaced": 4,
    }
    page = weather.render_page(_cov_result(200, 150), now=NOW, geo_prevalence=geo_prev)
    assert 'class="geo-prevalence"' in page
    assert "Where resistance is highest (by share, not just counts)" in page
    assert "USA: California" in page and "30.0%" in page
    assert "95% CI 18.0&ndash;45.0%" in page and "12 of 40 isolates" in page
    # The high-rate place whose CI lower bound clears the national rate is flagged; the other isn't.
    assert 'USA: California <span class="elevated">above national</span>' in page
    assert 'United Kingdom <span class="elevated">' not in page
    # The footnote names the national rate and the pooled-out small places (nothing hidden).
    assert "national rate is 10.0%" in page
    assert "2 smaller places" in page
    # Ranked: the higher-prevalence place comes first.
    assert page.index("USA: California") < page.index("United Kingdom")


def test_geo_prevalence_absent_without_data():
    assert 'class="geo-prevalence"' not in weather.render_page(_cov_result(200, 150), now=NOW)
    empty = {"places": [], "pooled_small": {}, "headline_prevalence": None}
    assert 'class="geo-prevalence"' not in weather.render_page(
        _cov_result(200, 150), now=NOW, geo_prevalence=empty)


def test_trend_section_surfaces_vus_separately_from_resistant():
    # #192: non-panel variants of uncertain significance are tracked, never headlined as resistant.
    trend = _trend(("2026-04-11", 300, 8), ("2026-04-14", 3375, 128))
    trend[-1]["n_vus"] = 69
    page = weather.render_page(_cov_result(200, 150), now=NOW, release_trend=trend)
    assert "128 of 3375 isolates" in page and "3.8%" in page
    assert "69 resolved isolates carry a non-panel FKS1 variant of uncertain significance" in page
    assert "5.8%" not in page


def test_trend_section_omits_vus_note_when_none():
    page = weather.render_page(
        _cov_result(200, 150), now=NOW,
        release_trend=_trend(("2026-04-11", 300, 8), ("2026-04-14", 3375, 128)))
    assert "uncertain significance" not in page


# ---- Presentation / accessibility / plain-language polish --------------------

def test_page_offers_light_scheme_and_print_styles():
    # The page must be legible on a bright clinical screen and on paper, not dark-only.
    page = weather.render_page(_empty_result(), now=NOW)
    assert "color-scheme: light dark" in page
    assert "@media (prefers-color-scheme: light)" in page
    assert "@media print" in page


def test_page_uses_a_main_landmark():
    # Content sits inside a <main> landmark, not a bare <div>, for screen-reader navigation.
    page = weather.render_page(_empty_result(), now=NOW)
    assert "<main class=\"wrap\">" in page
    assert "</main>" in page
    assert "<div class=\"wrap\">" not in page


def test_key_value_blocks_are_semantic_description_lists():
    # Alert-card and section key/value rows render as <dl>/<dt>/<dd>, not styled <div>s.
    page = weather.render_page(_fks1_hit_result(), now=NOW)
    assert "<dl class=\"grid\">" in page
    assert "<dt>" in page and "<dd>" in page
    assert 'class="k"' not in page and 'class="v"' not in page


def test_plain_language_lede_frames_the_page_for_a_lay_reader():
    # A zero-background reader gets a plain "what this is / not a clinical result" lede up front.
    page = weather.render_page(_empty_result(), now=NOW)
    assert 'class="lede"' in page
    assert "early-warning triage aid" in page
    assert "not a clinical test result" in page


def test_jargon_terms_carry_a_plain_language_gloss():
    # "Prevalence" and "variants of uncertain significance" each get a one-line plain gloss.
    trend_page = weather.render_page(
        _cov_result(200, 150), now=NOW,
        release_trend=_trend(("2026-04-11", 300, 15), ("2026-04-14", 3375, 197)))
    assert 'class="gloss"' in trend_page
    assert "share of successfully called isolates" in trend_page
    vus_page = weather.render_page(
        _cov_result(200, 150), now=NOW, vus_census=_vus(("F635I", 6)))
    assert 'class="gloss"' in vus_page
    assert "<em>not yet</em> known" in vus_page
    assert "not counted as resistant" in vus_page
