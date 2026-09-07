"""Tests for the caller concordance metrics (openafr/concordance.py).

These lock the guarantees the external-truth validation write-up leans on:
  * unresolved isolates (uncalled / refused / failed) are EXCLUDED from every denominator
    and surfaced as a separate count -- never scored as a negative,
  * the panel confusion matrix (TP/TN/FP/FN) is built from expected-vs-called token
    presence, and sensitivity/specificity/accuracy are the exact ratios with Wilson CIs,
  * identity is measured only among detected positives and demands every expected token
    (S639F != S639P is a miss, not a pass),
  * per-token sensitivity is per expected-carrier, and
  * resistance-detection is anchored on the PHENOTYPE and legitimately sits below panel
    sensitivity when a resistant isolate carries a non-panel mechanism.
No network, no reads: everything runs against synthetic per-isolate evaluations.
"""
from openafr import concordance as cc


def _rec(expected, called, resolved=True, resistant=None):
    return {"expected": expected, "called": called,
            "resolved": resolved, "resistant": resistant}


def test_norm_treats_dash_and_blank_as_wildtype():
    assert cc._norm("-") == set()
    assert cc._norm("") == set()
    assert cc._norm(None) == set()
    assert cc._norm("S639F") == {"S639F"}
    assert cc._norm("S639F, s639p") == {"S639F", "S639P"}
    assert cc._norm(["S639F"]) == {"S639F"}


def test_unresolved_excluded_from_denominators_and_counted():
    recs = [
        _rec("S639F", "S639F"),          # TP
        _rec("-", "-"),                  # TN
        _rec("S639P", "", resolved=False),   # unresolved -> excluded, NOT a FN
    ]
    s = cc.evaluate(recs)
    assert s["n_total"] == 3
    assert s["n_resolved"] == 2
    assert s["n_unresolved"] == 1
    # the unresolved S639P isolate must not appear as a false negative
    assert s["confusion"] == {"tp": 1, "tn": 1, "fp": 0, "fn": 0}
    assert s["sensitivity"]["n"] == 1        # only the resolved positive
    assert s["specificity"]["n"] == 1


def test_confusion_and_sensitivity_specificity():
    recs = [
        _rec("S639F", "S639F"),   # TP
        _rec("S639P", ""),        # FN (missed)
        _rec("-", "S639Y"),       # FP (spurious)
        _rec("-", "-"),           # TN
        _rec("-", "-"),           # TN
    ]
    s = cc.evaluate(recs)
    assert s["confusion"] == {"tp": 1, "tn": 2, "fp": 1, "fn": 1}
    assert s["sensitivity"]["k"] == 1 and s["sensitivity"]["n"] == 2
    assert s["sensitivity"]["point"] == 0.5
    assert s["specificity"]["k"] == 2 and s["specificity"]["n"] == 3
    assert s["accuracy"]["k"] == 3 and s["accuracy"]["n"] == 5
    # Wilson CI present and inside [0,1]
    lo, hi = s["sensitivity"]["ci95"]
    assert 0.0 <= lo <= 0.5 <= hi <= 1.0


def test_identity_requires_exact_token_among_detected():
    recs = [
        _rec("S639F", "S639F"),   # detected + exact -> identity match
        _rec("S639P", "S639Y"),   # detected but wrong token -> detected TP, identity MISS
        _rec("S639Y", ""),        # missed entirely -> not counted in identity
    ]
    s = cc.evaluate(recs)
    # both detected positives count toward identity denominator; only one matches
    assert s["identity"]["k"] == 1 and s["identity"]["n"] == 2
    # the wrong-token call is still a detection TP (a panel hit was emitted)
    assert s["confusion"]["tp"] == 2 and s["confusion"]["fn"] == 1


def test_per_token_sensitivity_is_per_carrier():
    recs = [
        _rec("S639F", "S639F"),
        _rec("S639F", ""),        # a second S639F carrier, missed
        _rec("S639P", "S639P"),
    ]
    s = cc.evaluate(recs)
    assert s["per_token"]["S639F"]["k"] == 1 and s["per_token"]["S639F"]["n"] == 2
    assert s["per_token"]["S639P"]["k"] == 1 and s["per_token"]["S639P"]["n"] == 1


def test_resistance_detection_below_panel_sensitivity_for_nonpanel_mechanism():
    # a resistant isolate whose truth token is a NON-panel mechanism (empty panel set):
    # the HS1 caller correctly emits no panel hit, so it is an honest resistance miss but
    # NOT a panel false-negative.
    recs = [
        _rec("S639F", "S639F", resistant=True),   # panel TP + resistance flagged
        _rec("-", "", resistant=True),            # non-panel resistant -> resistance miss
        _rec("-", "-", resistant=False),          # susceptible wild-type
    ]
    s = cc.evaluate(recs)
    # panel sensitivity sees only the one panel-positive isolate
    assert s["sensitivity"]["k"] == 1 and s["sensitivity"]["n"] == 1
    # resistance detection is over BOTH resistant isolates; only one is flagged
    assert s["resistance_detection"]["k"] == 1 and s["resistance_detection"]["n"] == 2
    # the non-panel resistant isolate is a TN at the panel (no expected panel token), not FN;
    # together with the susceptible wild-type that is two panel TNs
    assert s["confusion"]["fn"] == 0 and s["confusion"]["tn"] == 2


def test_resistance_detection_absent_when_no_phenotype_given():
    s = cc.evaluate([_rec("S639F", "S639F"), _rec("-", "-")])
    assert s["resistance_detection"]["n"] == 0
    assert s["resistance_detection"]["point"] is None


def test_format_report_runs_and_mentions_key_lines():
    s = cc.evaluate([_rec("S639F", "S639F", resistant=True),
                     _rec("-", "-", resistant=False),
                     _rec("S639P", "", resolved=False)])
    out = cc.format_report(s)
    assert "unresolved" in out
    assert "sensitivity" in out
    assert "resistance detection" in out
