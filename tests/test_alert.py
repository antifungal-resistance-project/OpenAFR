"""Tests for alert composition (openafr/alert.py, issue #25).

These lock the delivery tip of the early-warning loop -- the one screen a
non-bioinformatician acts on. The guarantees they enforce:

  * the join is complete: every emergence signal becomes exactly one alert carrying
    its structural fit verdict, the regions it spans, and a first-seen date,
  * first-seen is computed over ALL carriers of a mutation (not just the recent
    window), and keeps collected vs NCBI-visible dates separate -- the backlog tell,
  * priority is derived, never invented: a genuine weaker-fit signal is ACT-NOW, a
    genuine no-fit-story signal is WATCH, an archival-backlog signal is CONTEXT,
  * the honest "no calls yet" path yields a note and zero alerts, not a false all-clear,
  * both renderers (Markdown, JSON) round-trip the same self-describing dict.
No network: everything runs against the synthetic #22 fixture + checked-in structure.
"""
import json
import pathlib
import sys

from openafr import alert

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
from detect_emergence import _demo_records  # noqa: E402  (the shared #22 fixture)

AS_OF = "2024-12-31"


def _compose():
    return alert.compose_alerts(_demo_records(), AS_OF, window_days=180,
                                min_count=3, min_delta=0.05)


def _by_mut(result):
    return {a["mutation"]: a for a in result["alerts"]}


# ---- the join is complete ----------------------------------------------------

def test_every_signal_becomes_one_alert_with_structural_verdict():
    from openafr import emergence as em
    signals = em.detect_emergence(_demo_records(), AS_OF, window_days=180,
                                  min_count=3, min_delta=0.05)["signals"]
    result = _compose()
    assert len(result["alerts"]) == len(signals)
    assert {a["mutation"] for a in result["alerts"]} == {s["mutation"] for s in signals}
    for a in result["alerts"]:
        # every alert carries a real fit verdict with an explicit confidence
        assert a["structural"]["fluconazole_fit_verdict"]
        assert a["structural"]["confidence"] in ("measured", "estimate", "none")
        assert a["priority"] in alert.PRIORITY


def test_alert_carries_the_epidemiology_the_detector_dropped():
    a = _by_mut(_compose())["F126L"]
    # F126L is on 4 freshly-collected Pakistan isolates in the fixture
    assert a["emergence"]["recent_count"] == 4
    assert a["regions"] == ["Pakistan"]
    assert a["n_regions"] == 1


# ---- first-seen: all carriers, collected vs visible kept apart ---------------

def test_first_seen_spans_all_carriers_not_just_recent_window():
    # Y132F co-occurs on the 2019-collected backlog isolates as well as the 2024
    # baseline; first-collected must reflect the EARLIEST sample we hold across every
    # carrier (2019), not the recent window or the baseline-only view.
    a = _by_mut(_compose())["Y132F"]
    assert a["first_seen"]["collected"] == "2019"
    assert a["first_seen"]["visible"] == "2024-03-01"
    # carrier total counts across the whole store, not just the flagged window
    assert a["n_carriers_total"] >= a["emergence"]["recent_count"]


def test_first_seen_collected_and_visible_are_distinct_fields():
    # TR34's archival tell: collected 2019 but not visible until 2025.
    a = _by_mut(_compose())["TR34"]
    assert a["first_seen"]["collected"] == "2019"
    assert a["first_seen"]["visible"].startswith("2025")


def test_min_date_handles_bare_years_and_blanks():
    assert alert._min_date(["2024-03-01", "2019", ""]) == "2019"
    assert alert._min_date(["", "  "]) == ""
    assert alert._min_date([]) == ""


# ---- priority is derived, never invented -------------------------------------

def test_priority_tiers_match_upstream_verdicts():
    by = _by_mut(_compose())
    # Y132F: genuine + measured weaker-fit dock -> ACT-NOW
    assert by["Y132F"]["priority"] == "act-now"
    assert by["Y132F"]["structural"]["direction"] == "weaker-fit"
    # F126L / K143R: genuine but no confident fit call -> WATCH
    assert by["F126L"]["priority"] == "watch"
    assert by["K143R"]["priority"] == "watch"
    # TR34: possible archival backlog -> CONTEXT regardless of fit
    assert by["TR34"]["priority"] == "context"
    assert by["TR34"]["emergence"]["possible_backlog"] is True


def test_alerts_ranked_act_now_first():
    ranks = [alert.PRIORITY[a["priority"]] for a in _compose()["alerts"]]
    assert ranks == sorted(ranks)
    assert _compose()["alerts"][0]["priority"] == "act-now"


def test_backlog_does_not_suppress_the_alert_only_deprioritises_it():
    # We never hide a backlog signal; it is reported in full, just as CONTEXT.
    tr34 = _by_mut(_compose())["TR34"]
    assert tr34["emergence"]["backlog_reason"]
    assert tr34["priority"] == "context"


# ---- the honest "no calls yet" path ------------------------------------------

def test_no_calls_yields_note_and_no_alerts():
    recs = [{"isolate_key": f"P{i}", "target_acc": f"P{i}.1",
             "target_creation_date": "2025-02-01", "collection_date": "2024",
             "country": "India", "erg11_call": "", "run_acc": f"SRR{i}",
             "resistance_source": "pending:sra-recaller"} for i in range(5)]
    result = alert.compose_alerts(recs, AS_OF)
    assert result["alerts"] == []
    assert "note" in result and "cannot detect" in result["note"]
    md = alert.render_markdown(result)
    assert "No calls to assess" in md
    assert "No emergence alerts" in md


# ---- rendering round-trips the self-describing dict --------------------------

def test_render_json_is_valid_and_lossless():
    result = _compose()
    parsed = json.loads(alert.render_json(result))
    assert parsed["alerts"][0]["mutation"] == result["alerts"][0]["mutation"]
    assert parsed["coverage"] == result["coverage"]


def test_render_markdown_is_one_screen_signal_not_dump():
    md = alert.render_markdown(_compose())
    assert md.startswith("# C. auris azole early-warning alert")
    # the headline so-what for each mutation is present
    for mut in ("Y132F", "F126L", "K143R", "TR34"):
        assert f"`{mut}`" in md
    # the act-now badge and the measured Y132F verdict surface
    assert "ACT-NOW" in md
    assert "WEAKER azole fit" in md
    # caveats are carried, not stripped
    assert "caveat" in md
