"""Tests for the RUO verdict report (openafr/report.py, issue #138).

These lock the three things #138 exists to guarantee, on top of the #133 verdict contract:
  * the RUO / no-clinical-claim disclaimer is baked into every report (JSON + Markdown),
    never fine print, and a verdict whose scope disagrees is refused;
  * the UNCHARACTERIZED verdict is surfaced prominently -- its own top-level JSON list and a
    callout near the TOP of the Markdown, not buried down a table (#135 moat);
  * NO_KNOWN_MARKER is never rendered as 'susceptible', and no probability/score/MIC field is
    ever introduced (the day-one contract is concordance-only, #132).
Everything runs against verdict objects built by openafr.verdict -- no reads, no network.
"""
import json

from openafr import report as rpt
from openafr import verdict as v


def _call(tokens=(), panel_hits=(), uncalled_panel=()):
    return {"tokens": list(tokens), "panel_hits": list(panel_hits),
            "uncalled_panel": list(uncalled_panel)}


def test_disclaimer_baked_into_dict_and_both_renderings():
    r = rpt.build_report("iso-1", [v.azole_verdict(_call())])
    assert r["disclaimer"] == rpt.DISCLAIMER
    assert r["scope"] == v.RUO_SCOPE
    # JSON carries it verbatim.
    assert "RESEARCH USE ONLY" in rpt.render_json(r)
    # Markdown puts it in the header block, before any per-class verdict.
    md = rpt.render_markdown(r)
    assert "RESEARCH USE ONLY" in md
    assert md.index("RESEARCH USE ONLY") < md.index("## ")


def test_scope_drift_is_refused():
    good = v.azole_verdict(_call())
    bad = dict(v.echinocandin_verdict(_call()))
    bad["scope"] = "totally clinical, trust me"
    try:
        rpt.build_report("iso-1", [good, bad])
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError on disagreeing scope")


def test_empty_verdicts_rejected():
    try:
        rpt.build_report("iso-1", [])
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError on no verdicts")


def test_uncharacterized_surfaced_prominently():
    az = v.azole_verdict(_call(tokens=["T123I"]))          # UNCHARACTERIZED_VARIANT
    ec = v.echinocandin_verdict(_call())                    # NO_KNOWN_MARKER
    r = rpt.build_report("iso-2", [az, ec])
    # Top-level summary list, not something you have to scan the verdicts to find.
    assert r["summary"]["uncharacterized"] == ["azole"]
    md = rpt.render_markdown(r)
    # The callout appears before the per-class detail blocks.
    assert "Uncharacterized variant(s) present" in md
    assert md.index("Uncharacterized variant(s) present") < md.index("## azole")


def test_no_known_marker_never_reads_as_susceptible():
    r = rpt.build_report("iso-3", [v.azole_verdict(_call())])
    rationale = r["verdicts"][0]["rationale"].lower()
    assert "not establish susceptibility" in rationale or "not a susceptibility" in \
        r["verdicts"][0]["headline"].lower()
    assert "suscept" in rationale                # it must explicitly address susceptibility
    # ...and only to deny it -- the word 'susceptible' never stands as a positive call.
    md = rpt.render_markdown(r)
    assert "NOT a susceptibility call" in md


def test_resistance_marker_rationale_names_the_token():
    r = rpt.build_report("iso-4",
                         [v.azole_verdict(_call(tokens=["Y132F"], panel_hits=["Y132F"]))])
    item = r["verdicts"][0]
    assert item["verdict"] == v.RESISTANCE_MARKER_DETECTED
    assert "Y132F" in item["rationale"]
    assert "Y132F" in rpt.render_markdown(r)


def test_unresolved_is_reported_as_excluded_not_negative():
    r = rpt.build_report("iso-5",
                         [v.azole_verdict(unresolved_reason="consensus not multiple of 3")])
    item = r["verdicts"][0]
    assert item["verdict"] == v.UNRESOLVED
    assert r["summary"]["unresolved"] == ["azole"]
    md = rpt.render_markdown(r)
    assert "excluded" in md.lower()
    assert "multiple of 3" in md


def test_no_probability_or_score_field_anywhere_in_report():
    az = v.azole_verdict(_call(tokens=["T123I"]))
    r = rpt.build_report("iso-6", [az])
    # No calibrated NUMBER or score/MIC field sneaks in via the report layer (contract
    # rule 1). The word 'probability' is allowed ONLY in the disclaimer/confidence_model,
    # where it is explicitly denied ("does not emit a calibrated probability").
    blob = rpt.render_json(r).lower()
    for banned in ('"score"', "mic_value", '"mic"', "0.9", "0.8"):
        assert banned not in blob
    # Every mention of 'probability' must be a denial, never a reported value.
    for line in blob.splitlines():
        if "probability" in line:
            assert "no calibrated probability" in line or "not emit a calibrated" in line


def test_build_report_does_not_mutate_inputs():
    az = v.azole_verdict(_call(tokens=["Y132F"], panel_hits=["Y132F"]))
    before = json.dumps(az, sort_keys=True)
    rpt.build_report("iso-7", [az])
    assert json.dumps(az, sort_keys=True) == before   # no headline/rationale leaked back in


def test_structural_best_guess_surfaced_when_present():
    # A verdict carrying an explicit ERG11 structural block renders it, flagged calibrated-low.
    structural = {
        "kind": "uncharacterized_best_guess",
        "flag": "uncharacterized",
        "confidence": "low",
        "basis": "modeled pocket (#13)",
        "summary": "best-guess for T123I; low-confidence mechanism inference",
        "variants": [{
            "token": "T123I", "evidence_class": "no-call", "confidence": "none",
            "fluconazole_fit_verdict": "THR123ILE: NO CONFIDENT CALL (low prior).",
            "lost_contacts": [],
        }],
    }
    az = v.azole_verdict(_call(tokens=["T123I"]), structural=structural)
    md = rpt.render_markdown(rpt.build_report("iso-8", [az]))
    assert "Structural best-guess" in md
    assert "confidence=low" in md
    assert "T123I" in md
