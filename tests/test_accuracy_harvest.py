"""Tests for the #137 accuracy-fixture HARVEST in scripts/validate_fks1_concordance.py.

The concordance grader can, opt-in (`--emit-accuracy-fixture`), emit the paired
verdict/phenotype fixture the #137 accuracy grader (scripts/validate_diagnostic_accuracy.py)
needs, from the SAME reads->call_windows pass. The execution needs the cloud host, but the
row-assembly logic -- caller result -> `called_verdict` enum + fixture row -- is pure and is
pinned here against synthetic call_windows dicts, so a wrong column or a mis-mapped verdict
can't reach the expensive run undetected.
"""
import importlib.util
import pathlib

from openafr import verdict as V

_SCRIPT = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "validate_fks1_concordance.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("validate_fks1_concordance", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _result(panel_hits=(), tokens=(), uncalled_panel=()):
    """A synthetic openafr.fks1_caller.call_windows result (the fields the harvest reads)."""
    return {
        "panel_hits": list(panel_hits),
        "tokens": list(tokens),
        "uncalled_panel": list(uncalled_panel),
        "windows": {"HS1": {"status": "called"}, "HS2": {"status": "called"},
                    "HS3": {"status": "called"}},
    }


def _row_dict(row, script):
    return dict(zip(script.ACCURACY_HEADER, row))


def test_panel_hit_resistant_marker_detected():
    s = _load_script()
    row = s._accuracy_row("B11211", "SRR28270963", _result(panel_hits=["S639F"],
                          tokens=["S639F"]), "R", "PMC12323592")
    d = _row_dict(row, s)
    assert d["called_verdict"] == V.RESISTANCE_MARKER_DETECTED
    assert d["susceptibility"] == "R"
    assert (d["isolate"], d["species"], d["gene"], d["drug_class"]) == (
        "B11211", "Candida auris", "FKS1", "echinocandin")
    assert d["run_acc"] == "SRR28270963"
    assert d["citation"] == "PMC12323592"
    assert "panel_hits=S639F" in d["note"]


def test_wild_type_no_known_marker():
    s = _load_script()
    d = _row_dict(s._accuracy_row("B11098", "SRR3883438", _result(), "S", "cite"), s)
    assert d["called_verdict"] == V.NO_KNOWN_MARKER
    assert d["susceptibility"] == "S"


def test_non_panel_token_is_uncharacterized():
    s = _load_script()
    # W691C is emitted verbatim in the HS3 window but is NOT panel-tagged (dropped for
    # circularity -- its only report is the benchmark itself): reported, not tagged.
    d = _row_dict(s._accuracy_row("B20717", "SRRx", _result(tokens=["W691C"]), "R", "c"), s)
    assert d["called_verdict"] == V.UNCHARACTERIZED_VARIANT


def test_uncalled_panel_residue_is_unresolved():
    s = _load_script()
    d = _row_dict(s._accuracy_row("Bxx", "SRRy", _result(uncalled_panel=[639]), "R", "c"), s)
    assert d["called_verdict"] == V.UNRESOLVED


def test_orchestration_failure_is_unresolved():
    s = _load_script()
    d = _row_dict(s._accuracy_row("Bzz", "SRRz", None, "R", "c"), s)
    assert d["called_verdict"] == V.UNRESOLVED
    assert "orchestration_failed" in d["note"]


def test_written_fixture_parses_with_accuracy_schema(tmp_path):
    s = _load_script()
    rows = [
        s._accuracy_row("B11211", "SRR1", _result(panel_hits=["S639F"], tokens=["S639F"]), "R", "c"),
        s._accuracy_row("B11098", "SRR2", _result(), "S", "c"),
    ]
    out = s._write_accuracy_fixture(tmp_path / "acc.tsv", rows)
    # Parse exactly like the accuracy grader: skip '#'-comments, first row is the header.
    parsed, header = [], None
    for line in pathlib.Path(out).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if header is None:
            header = f
            continue
        parsed.append(dict(zip(header, f)))
    assert header == s.ACCURACY_HEADER
    assert len(parsed) == 2
    assert parsed[0]["called_verdict"] == V.RESISTANCE_MARKER_DETECTED
    assert parsed[1]["called_verdict"] == V.NO_KNOWN_MARKER
    # every graded column the accuracy grader reads is populated
    for r in parsed:
        assert r["drug_class"] == "echinocandin"
        assert r["susceptibility"] in ("R", "S")
        assert r["called_verdict"] in V.VERDICTS
