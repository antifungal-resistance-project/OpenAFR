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
