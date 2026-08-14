"""Tests for the delivery layer (openafr/delivery.py, issue #26).

These lock the promise the schedule makes: the signal finds a human, on a cadence,
*minimally*. The guarantees they enforce:

  * new news is defined narrowly: a never-before-delivered mutation, or a genuine
    escalation above the most urgent tier ever delivered for that mutation,
  * an unchanged re-run is a QUIET NO-OP -- no alerts surface, nothing is delivered,
  * a de-escalation (act-now -> watch) is NOT news and stays quiet,
  * state is most-urgent-ever and byte-stable, so its git history is the delivery log,
  * the honest "no calls yet" path is a no-op, not a false all-clear,
  * the digest names what's new and carries the full brief; the issue title is scannable.
No network: everything runs against the synthetic #22 fixture.
"""
import json
import pathlib
import sys

from openafr import alert, delivery

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
from detect_emergence import _demo_records  # noqa: E402  (the shared #22 fixture)

AS_OF = "2024-12-31"


def _compose(records=None):
    return alert.compose_alerts(records or _demo_records(), AS_OF, window_days=180,
                                min_count=3, min_delta=0.05)


def _clear_tr34_backlog(records):
    """The fixture with TR34 freshly collected -> backlog guard clears -> context rises."""
    out = []
    for r in records:
        r = dict(r)
        if "TR34" in r["erg11_call"]:
            r["collection_date"] = "2024"
        out.append(r)
    return out


# ---- new-news detection ------------------------------------------------------

def test_empty_state_makes_every_alert_new():
    result = _compose()
    fresh = delivery.new_alerts(result, delivery.load_state(""))
    assert {a["mutation"] for a in fresh} == {a["mutation"] for a in result["alerts"]}
    assert not delivery.is_noop(result, delivery.load_state(""))


def test_unchanged_rerun_is_a_quiet_noop():
    result = _compose()
    state = delivery.next_state(result, delivery.load_state(""))
    assert delivery.new_alerts(result, state) == []
    assert delivery.is_noop(result, state)


def test_escalation_refires_context_to_watch():
    base = _compose()
    state = delivery.next_state(base, delivery.load_state(""))
    # sanity: TR34 was delivered as context
    assert state["delivered"]["TR34"] == "context"

    escalated = _compose(_clear_tr34_backlog(_demo_records()))
    fresh = delivery.new_alerts(escalated, state)
    muts = {a["mutation"] for a in fresh}
    assert muts == {"TR34"}
    assert next(a for a in fresh if a["mutation"] == "TR34")["priority"] == "watch"


def test_deescalation_is_not_news():
    # A mutation already delivered at act-now that now reads watch must stay quiet.
    result = _compose()
    y = next(a for a in result["alerts"] if a["mutation"] == "Y132F")
    assert y["priority"] == "act-now"
    state = {"delivered": {"Y132F": "act-now"}, "updated": AS_OF}
    demoted = {"alerts": [dict(y, priority="watch")], "as_of": AS_OF}
    assert delivery.new_alerts(demoted, state) == []


# ---- state semantics ---------------------------------------------------------

def test_state_keeps_most_urgent_tier_ever():
    # deliver act-now, then a later run reads the same mutation at watch: memory holds.
    state = {"delivered": {"Y132F": "act-now"}, "updated": AS_OF}
    watch_run = {"alerts": [{"mutation": "Y132F", "priority": "watch"}], "as_of": AS_OF}
    nxt = delivery.next_state(watch_run, state)
    assert nxt["delivered"]["Y132F"] == "act-now"


def test_state_remembers_vanished_mutations():
    state = {"delivered": {"GONE": "watch"}, "updated": AS_OF}
    result = _compose()  # does not contain GONE
    nxt = delivery.next_state(result, state)
    assert nxt["delivered"]["GONE"] == "watch"


def test_state_roundtrips_byte_stable():
    result = _compose()
    nxt = delivery.next_state(result, delivery.load_state(""))
    text = delivery.dumps_state(nxt)
    assert text.endswith("\n")
    assert delivery.dumps_state(delivery.load_state(text)) == text
    assert json.loads(text)["delivered"]  # valid JSON with content


def test_load_state_tolerates_blank_and_bad_tokens():
    assert delivery.load_state("")["delivered"] == {}
    assert delivery.load_state(None)["delivered"] == {}
    got = delivery.load_state('{"delivered": {"X": "bogus", "Y": "watch"}}')
    assert got["delivered"] == {"Y": "watch"}  # unknown tier dropped, not trusted


# ---- honest no-call path -----------------------------------------------------

def test_no_calls_yet_is_a_noop_not_a_false_all_clear():
    # Uncalled records -> composer emits a note and zero alerts -> nothing delivered.
    records = [dict(r, erg11_call="") for r in _demo_records()]
    result = _compose(records)
    assert result["alerts"] == []
    assert "note" in result
    assert delivery.is_noop(result, delivery.load_state(""))


# ---- rendering ---------------------------------------------------------------

def test_digest_names_whats_new_and_carries_full_brief():
    result = _compose()
    fresh = delivery.new_alerts(result, delivery.load_state(""))
    digest = delivery.render_digest(result, fresh, delivery.load_state(""))
    assert "What's new this run" in digest
    assert "new alert" in digest              # brand-new mutations this run
    assert "# C. auris azole early-warning alert" in digest  # the full composed brief
    for a in result["alerts"]:
        assert a["mutation"] in digest


def test_digest_shows_escalation_movement():
    base = _compose()
    state = delivery.next_state(base, delivery.load_state(""))
    escalated = _compose(_clear_tr34_backlog(_demo_records()))
    fresh = delivery.new_alerts(escalated, state)
    digest = delivery.render_digest(escalated, fresh, state)
    assert "CONTEXT → WATCH" in digest


def test_issue_title_is_scannable():
    result = _compose()
    fresh = delivery.new_alerts(result, delivery.load_state(""))
    title = delivery.issue_title(result, fresh)
    assert title.startswith("[early-warning]")
    assert "ACT-NOW" in title  # most urgent tier surfaced
    assert AS_OF in title
