"""Tests for the reproducibility / drift harness (scripts/repro.py, issue #89).

The harness is the thing that fails a PR when the frozen surface drifts, so its own
pass/fail logic has to be trustworthy. Two things are pinned here:

  * on the real repo, all three checks pass (the committed tree is internally consistent),
  * each check FAILS when its target is tampered — a drift-detector that cannot detect
    drift is worse than none, because it launders a false sense of safety.

The tamper cases use temp fixtures / monkeypatch so they never touch the committed files.
"""
import hashlib
import pathlib

import pytest

from scripts import repro


def _sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


# --- the real repo is internally consistent (integration) -------------------------

def test_protocol_check_passes_on_real_repo():
    ok, _ = repro.check_protocol()
    assert ok


def test_preregistrations_check_passes_on_real_repo():
    ok, lines = repro.check_preregistrations()
    assert ok, lines


def test_scorecard_reproduces_on_real_repo():
    ok, lines = repro.check_scorecard_reproduces()
    assert ok, lines


# --- each check detects its own tamper --------------------------------------------

def test_protocol_check_detects_drift(monkeypatch):
    monkeypatch.setattr(repro.protocol, "PROTOCOL_SHA", "0" * 64)
    ok, lines = repro.check_protocol()
    assert not ok
    assert any("DRIFTED" in l for l in lines)


def test_preregistration_check_detects_a_wrong_hash(tmp_path, monkeypatch):
    # A sha256 file whose pinned hash does not match its (real) target = tamper.
    bad = tmp_path / "PREREG_fake.sha256"
    bad.write_text(f"{'0' * 64}  README.md\n")
    monkeypatch.setattr(repro.glob, "glob", lambda _pat: [str(bad)])
    ok, lines = repro.check_preregistrations()
    assert not ok
    assert any("DRIFTED" in l for l in lines)


def test_preregistration_check_detects_a_missing_target(tmp_path, monkeypatch):
    bad = tmp_path / "PREREG_fake.sha256"
    bad.write_text(f"{'0' * 64}  work/does_not_exist.md\n")
    monkeypatch.setattr(repro.glob, "glob", lambda _pat: [str(bad)])
    ok, lines = repro.check_preregistrations()
    assert not ok
    assert any("MISSING" in l for l in lines)


def test_preregistration_check_passes_with_correct_hash(tmp_path, monkeypatch):
    good = tmp_path / "PREREG_ok.sha256"
    readme = repro.ROOT / "README.md"
    good.write_text(f"{_sha(readme)}  README.md\n")
    monkeypatch.setattr(repro.glob, "glob", lambda _pat: [str(good)])
    ok, _ = repro.check_preregistrations()
    assert ok


def test_sha256_file_parser_handles_two_space_separator():
    entries = list(repro._parse_sha256_file(
        _write_lines("abc123  path/one.md\n\ndef456  data/two.smi\n")))
    assert entries == [("abc123", "path/one.md"), ("def456", "data/two.smi")]


def test_main_returns_nonzero_when_a_check_fails(monkeypatch, capsys):
    monkeypatch.setattr(repro, "CHECKS", [("always fails", lambda: (False, ["  nope"]))])
    rc = repro.main(["--quiet"])
    assert rc == 1
    assert "FAIL" in capsys.readouterr().out


def test_main_returns_zero_when_all_pass(monkeypatch, capsys):
    monkeypatch.setattr(repro, "CHECKS", [("always ok", lambda: (True, ["  fine"]))])
    rc = repro.main([])
    assert rc == 0
    assert "PASS" in capsys.readouterr().out


# --- helper ------------------------------------------------------------------------

_TMP = []


def _write_lines(text):
    import tempfile
    f = tempfile.NamedTemporaryFile("w", suffix=".sha256", delete=False)
    f.write(text)
    f.close()
    _TMP.append(f.name)
    return f.name


def teardown_module(_module):
    import os
    for p in _TMP:
        os.unlink(p)
