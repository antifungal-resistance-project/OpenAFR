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

AS_OF = "2024-12-31"
NOW = "2026-09-17T08:30:00Z"


def _hit_result():
    return alert.compose_alerts(_demo_records(), AS_OF, window_days=180,
                                min_count=3, min_delta=0.05)


def _empty_result(note=None):
    r = {"as_of": AS_OF, "window_days": 180, "thresholds": {},
         "coverage": {}, "alerts": []}
    if note is not None:
        r["note"] = note
    return r


# ---- WATCHING / day-0 --------------------------------------------------------

def test_none_result_renders_watching_page():
    page = weather.render_page(None, now=NOW)
    assert "WATCHING" in page
    assert "No azole-resistance variant is over threshold" in page
    assert 'class="badge"' not in page  # no tier badge element when there are no alerts


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
