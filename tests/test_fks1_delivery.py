"""Tests for the FKS1/echinocandin delivery path (openafr/delivery.py, v2).

The delivery decision logic (new-vs-known, escalation, quiet no-op, committed state) is
gene-agnostic and covered by test_delivery.py. These lock the gene-aware half: an FKS1
result must render an echinocandin-titled, DETECTION-ONLY digest/issue -- never the
azole/ERG11 wording or a structural fit line -- while the ERG11 path stays unchanged.
No network: runs against the synthetic FKS1 fixture.
"""
import pathlib
import sys

from openafr import alert, delivery

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
from compose_alert import _demo_fks1_records  # noqa: E402

AS_OF = "2024-12-31"


def _compose(records=None):
    return alert.compose_fks1_alerts(records or _demo_fks1_records(), AS_OF,
                                     window_days=180, min_count=3, min_delta=0.05)


def test_decision_logic_still_works_for_fks1():
    result = _compose()
    fresh = delivery.new_alerts(result, delivery.load_state(""))
    assert {a["mutation"] for a in fresh} == {"S639F", "S639P"}
    # unchanged re-run is a quiet no-op
    state = delivery.next_state(result, delivery.load_state(""))
    assert delivery.is_noop(result, state)


def test_digest_is_echinocandin_titled_and_detection_only():
    result = _compose()
    fresh = delivery.new_alerts(result, delivery.load_state(""))
    digest = delivery.render_digest(result, fresh, delivery.load_state(""))
    assert "# C. auris echinocandin (FKS1) early-warning digest" in digest
    assert "# C. auris echinocandin (FKS1) early-warning alert" in digest  # full brief
    assert "Detection only." in digest
    assert "azole early-warning" not in digest
    assert "Structural fit** [" not in digest  # no ERG11 structural line


def test_issue_title_names_fks1_not_erg11():
    result = _compose()
    fresh = delivery.new_alerts(result, delivery.load_state(""))
    title = delivery.issue_title(result, fresh)
    assert title.startswith("[early-warning]")
    assert "C. auris FKS1" in title
    assert "ERG11" not in title
    assert AS_OF in title


def test_no_act_now_ever_in_the_fks1_title():
    # FKS1 is detection-only, so the loudest tier it can deliver is WATCH.
    result = _compose()
    fresh = delivery.new_alerts(result, delivery.load_state(""))
    title = delivery.issue_title(result, fresh)
    assert "ACT-NOW" not in title


def test_no_calls_yet_is_a_quiet_noop_for_fks1():
    recs = [dict(r, fks1_call="") for r in _demo_fks1_records()]
    result = _compose(recs)
    assert delivery.is_noop(result, delivery.load_state(""))


def test_erg11_digest_unchanged_by_the_gene_routing():
    # Regression guard: an ERG11 result (no gene marker) still renders the azole wording.
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
    from detect_emergence import _demo_records
    result = alert.compose_alerts(_demo_records(), AS_OF, window_days=180,
                                  min_count=3, min_delta=0.05)
    fresh = delivery.new_alerts(result, delivery.load_state(""))
    digest = delivery.render_digest(result, fresh, delivery.load_state(""))
    assert "# C. auris azole early-warning digest" in digest
    assert "C. auris ERG11" in delivery.issue_title(result, fresh)
