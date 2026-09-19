"""Tests for the verdict-level accuracy metrics (openafr/accuracy.py, issue #137).

These lock the guarantees the diagnostic-accuracy write-up leans on:
  * the verdict enum -> R/S call mapping is exactly the pre-registered one
    (RESISTANCE_MARKER_DETECTED->R, NO_KNOWN_MARKER->not-R, UNCHARACTERIZED->abstain,
    UNRESOLVED->exclude);
  * abstained and unresolved verdicts, and isolates with no phenotype, are EXCLUDED from the
    confusion matrix and surfaced as separate counts -- never scored as a negative, so a low
    VME cannot be bought by reclassifying hard isolates;
  * VME is denominated over the RESISTANT isolates and ME over the SUSCEPTIBLE ones (the CLSI
    convention), and sensitivity = 1 - VME, specificity = 1 - ME;
  * metrics are grouped per drug class (azole VME != echinocandin VME);
  * the pre-committed pass bar (apply_bar) gates on the Wilson UPPER bound, and reports
    UNDERPOWERED when an arm is too thin to read.
No network, no reads: everything runs against synthetic verdicts.
"""
import pytest

from openafr import accuracy as ac
from openafr import verdict as vd


def _rec(v, phenotype, drug_class="azole"):
    return {"verdict": v, "phenotype": phenotype, "drug_class": drug_class}


def test_phenotype_normalisation():
    assert ac._phenotype("R") == "R"
    assert ac._phenotype("susceptible") == "S"
    assert ac._phenotype(True) == "R"
    assert ac._phenotype(False) == "S"
    assert ac._phenotype(None) is None
    assert ac._phenotype("-") is None
    with pytest.raises(ValueError):
        ac._phenotype("maybe")


def test_disposition_covers_every_verdict_value():
    # every enum value in the contract must have a scoring disposition (no silent drop)
    assert set(ac._DISPOSITION) == set(vd.VERDICTS)


def test_marker_detected_on_resistant_is_true_positive():
    s = ac.evaluate([_rec(vd.RESISTANCE_MARKER_DETECTED, "R")])["by_drug_class"]["azole"]
    assert s["confusion"] == {"tp": 1, "fn": 0, "fp": 0, "tn": 0}
    assert s["very_major_error"]["k"] == 0 and s["very_major_error"]["n"] == 1


def test_no_marker_on_resistant_is_very_major_error():
    # the load-bearing dangerous error: NO_KNOWN_MARKER on a phenotypically resistant isolate
    s = ac.evaluate([_rec(vd.NO_KNOWN_MARKER, "R")])["by_drug_class"]["azole"]
    assert s["confusion"] == {"tp": 0, "fn": 1, "fp": 0, "tn": 0}
    assert s["very_major_error"]["k"] == 1 and s["very_major_error"]["n"] == 1
    assert s["very_major_error"]["point"] == 1.0
    assert s["sensitivity"]["point"] == 0.0        # sensitivity = 1 - VME


def test_marker_on_susceptible_is_major_error():
    s = ac.evaluate([_rec(vd.RESISTANCE_MARKER_DETECTED, "S")])["by_drug_class"]["azole"]
    assert s["confusion"] == {"tp": 0, "fn": 0, "fp": 1, "tn": 0}
    assert s["major_error"]["k"] == 1 and s["major_error"]["n"] == 1
    assert s["specificity"]["point"] == 0.0        # specificity = 1 - ME


def test_abstention_and_unresolved_and_no_phenotype_excluded():
    recs = [
        _rec(vd.RESISTANCE_MARKER_DETECTED, "R"),      # scored TP
        _rec(vd.UNCHARACTERIZED_VARIANT, "R"),         # abstain: excluded, NOT a VME
        _rec(vd.UNRESOLVED, "R"),                      # unresolved: excluded, NOT a VME
        _rec(vd.NO_KNOWN_MARKER, None),                # no phenotype: excluded
    ]
    s = ac.evaluate(recs)["by_drug_class"]["azole"]
    assert s["n_total"] == 4
    assert s["n_scored"] == 1
    assert s["n_abstained"] == 1
    assert s["n_unresolved"] == 1
    assert s["n_no_phenotype"] == 1
    # the abstained/unresolved resistant isolates must NOT inflate the FN count
    assert s["confusion"]["fn"] == 0
    assert s["very_major_error"]["n"] == 1          # only the one scored resistant isolate
    # abstention denominator = scored + abstained (had phenotype + a categorical verdict)
    assert s["abstention"]["k"] == 1 and s["abstention"]["n"] == 2


def test_vme_me_denominators_follow_clsi_arms():
    recs = [
        _rec(vd.RESISTANCE_MARKER_DETECTED, "R"),   # TP
        _rec(vd.NO_KNOWN_MARKER, "R"),              # FN (VME)
        _rec(vd.NO_KNOWN_MARKER, "R"),              # FN (VME)
        _rec(vd.RESISTANCE_MARKER_DETECTED, "S"),   # FP (ME)
        _rec(vd.NO_KNOWN_MARKER, "S"),              # TN
    ]
    s = ac.evaluate(recs)["by_drug_class"]["azole"]
    # VME over the 3 resistant isolates; ME over the 2 susceptible isolates
    assert s["very_major_error"]["k"] == 2 and s["very_major_error"]["n"] == 3
    assert s["major_error"]["k"] == 1 and s["major_error"]["n"] == 2
    assert s["sensitivity"]["point"] == pytest.approx(1 / 3)
    assert s["specificity"]["point"] == pytest.approx(1 / 2)
    assert s["categorical_agreement"]["k"] == 2 and s["categorical_agreement"]["n"] == 5


def test_grouped_per_drug_class_and_not_pooled():
    recs = [
        _rec(vd.NO_KNOWN_MARKER, "R", drug_class="azole"),          # azole VME
        _rec(vd.RESISTANCE_MARKER_DETECTED, "R", drug_class="echinocandin"),  # echino TP
    ]
    out = ac.evaluate(recs)
    assert out["drug_classes"] == ["azole", "echinocandin"]
    assert out["by_drug_class"]["azole"]["very_major_error"]["point"] == 1.0
    assert out["by_drug_class"]["echinocandin"]["very_major_error"]["point"] == 0.0


def test_evaluate_accepts_full_verdict_objects():
    v = vd.echinocandin_verdict(
        call={"tokens": ["S639F"], "panel_hits": ["S639F"], "uncalled_panel": []})
    assert v["verdict"] == vd.RESISTANCE_MARKER_DETECTED
    out = ac.evaluate([ac.record(v, "R")])
    assert out["by_drug_class"]["echinocandin"]["confusion"]["tp"] == 1


def test_record_requires_drug_class():
    with pytest.raises(ValueError):
        ac.evaluate([{"verdict": vd.NO_KNOWN_MARKER, "phenotype": "R"}])


def test_apply_bar_underpowered_when_arm_thin():
    s = ac.evaluate([_rec(vd.RESISTANCE_MARKER_DETECTED, "R")])["by_drug_class"]["azole"]
    bar = {"vme_point": 0.03, "vme_upper": 0.10, "me_point": 0.03,
           "min_resistant": 15, "min_susceptible": 15, "max_abstention": 0.30}
    outcome, reasons = ac.apply_bar(s, bar)
    assert outcome == "UNDERPOWERED"
    assert any("too thin" in r for r in reasons)


def test_apply_bar_fail_on_high_vme_and_pass_on_clean():
    bar = {"vme_point": 0.05, "vme_upper": 0.20, "me_point": 0.10,
           "min_resistant": 3, "min_susceptible": 3, "max_abstention": 0.50}
    # many resistant isolates all missed -> VME = 100% -> FAIL
    bad = ac.evaluate([_rec(vd.NO_KNOWN_MARKER, "R")] * 5
                      + [_rec(vd.NO_KNOWN_MARKER, "S")] * 5)["by_drug_class"]["azole"]
    assert ac.apply_bar(bad, bar)[0] == "FAIL"
    # all resistant caught, all susceptible clean -> PASS (with a lenient CI bar)
    lenient = dict(bar, vme_upper=0.60)
    good = ac.evaluate([_rec(vd.RESISTANCE_MARKER_DETECTED, "R")] * 5
                       + [_rec(vd.NO_KNOWN_MARKER, "S")] * 5)["by_drug_class"]["azole"]
    assert ac.apply_bar(good, lenient)[0] == "PASS"


def test_format_report_runs_and_mentions_error_names():
    out = ac.format_report(ac.evaluate([
        _rec(vd.RESISTANCE_MARKER_DETECTED, "R"),
        _rec(vd.NO_KNOWN_MARKER, "S"),
        _rec(vd.UNCHARACTERIZED_VARIANT, "R"),
    ]))
    assert "very major error" in out
    assert "major error" in out
    assert "abstention" in out
