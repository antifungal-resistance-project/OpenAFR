"""Tests for the FKS1/echinocandin DETECTION-ONLY alert path (openafr/alert.py, v2).

The FKS1 early-warning track deliberately claims no structural so-what: echinocandins
coordinate no metal, so the CYP51/heme-iron geometry that powers the ERG11 verdict does
not transfer to a glucan synthase. These lock that the composer therefore:

  * detects emergence over `fks1_call` (not `erg11_call`),
  * marks the result gene=FKS1 / detection_only=True and carries NO structural block,
  * derives only the watch/context tiers -- never act-now, which requires a structural
    weaker-fit verdict FKS1 does not have,
  * states the detection-only caveat on every alert and in the rendered brief, so an
    absent fit verdict is never mistaken for an all-clear,
  * keeps the honest "no calls yet" note path,
  * leaves the ERG11 compose_alerts path untouched.
No network: everything runs against a synthetic called snapshot.
"""
import json
import pathlib
import sys

from openafr import alert

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
from compose_alert import _demo_fks1_records  # noqa: E402

AS_OF = "2024-12-31"


def _compose():
    return alert.compose_fks1_alerts(_demo_fks1_records(), AS_OF, window_days=180,
                                     min_count=3, min_delta=0.05)


def _by_mut(result):
    return {a["mutation"]: a for a in result["alerts"]}


def test_result_is_marked_detection_only_for_fks1():
    r = _compose()
    assert r["gene"] == "FKS1"
    assert r["detection_only"] is True


def test_reads_the_fks1_call_column_not_erg11():
    # The fixture puts calls only in fks1_call; reading erg11_call would find nothing.
    recs = _demo_fks1_records()
    assert all("erg11_call" not in r for r in recs)
    assert _compose()["alerts"], "should flag FKS1 signals from the fks1_call column"


def test_no_alert_carries_a_structural_verdict():
    for a in _compose()["alerts"]:
        assert "structural" not in a
        assert a["detection_only"] is True


def test_no_act_now_tier_only_watch_and_context():
    tiers = {a["priority"] for a in _compose()["alerts"]}
    assert "act-now" not in tiers
    assert tiers <= {"watch", "context"}


def test_genuine_signal_is_watch_backlog_is_context():
    by = _by_mut(_compose())
    assert by["S639F"]["priority"] == "watch"          # genuine, freshly collected
    assert by["S639P"]["priority"] == "context"        # archival backlog
    assert by["S639P"]["emergence"]["possible_backlog"] is True


def test_alerts_ranked_watch_before_context():
    ranks = [alert.PRIORITY[a["priority"]] for a in _compose()["alerts"]]
    assert ranks == sorted(ranks)


def test_every_headline_states_the_detection_only_caveat():
    for a in _compose()["alerts"]:
        assert "detection only" in a["headline"]
        assert "no structural verdict" in a["headline"]


def test_epidemiology_is_carried_like_the_erg11_path():
    s639f = _by_mut(_compose())["S639F"]
    assert s639f["regions"] == ["Pakistan"]
    assert s639f["n_carriers_total"] == 4
    assert s639f["first_seen"]["collected"] == "2024"


def test_no_calls_yields_note_and_no_alerts():
    recs = [dict(r, fks1_call="") for r in _demo_fks1_records()]
    r = alert.compose_fks1_alerts(recs, AS_OF)
    assert "note" in r
    assert r["alerts"] == []
    md = alert.render_fks1_markdown(r)
    assert "No calls to assess" in md


def test_render_fks1_markdown_is_echinocandin_titled_and_detection_only():
    md = alert.render_fks1_markdown(_compose())
    assert md.startswith("# C. auris echinocandin (FKS1) early-warning alert")
    assert "Detection only." in md
    assert "Structural fit:** none claimed" in md
    # never emits an act-now badge for the detection-only track
    assert "ACT-NOW" not in md


def test_render_json_round_trips_the_fks1_dict():
    r = _compose()
    back = json.loads(alert.render_json(r))
    assert back["gene"] == "FKS1"
    assert back == r
