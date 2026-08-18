"""Tests for the frozen run protocol (openafr/protocol.py, task T2).

The protocol carries the numbers most tempting to silently tune into a pass — the
docking box, the search effort, and the gate's pass bar — so these lock both the
values and the tamper-detection that guards them.
"""
import json
import pathlib

import pytest

from openafr import protocol


def test_frozen_values_match_prereg():
    # These must equal work/PREREGISTRATION_run2.md exactly; the gate and the
    # docking script both read them from here.
    p = protocol.load()
    assert p["docking"] == {
        "center_x": 70.61, "center_y": 66.28, "center_z": 4.18,
        "size": 26, "exhaustiveness": 32, "num_modes": 20, "seed": 42,
    }
    assert p["gate"]["min_auc"] == 0.70
    assert p["gate"]["min_ef1"] == 5.0
    assert p["gate"]["ef_fractions"] == [0.01, 0.05, 0.10]
    assert p["gate"]["bedroc_alpha"] == 20.0


def test_frozen_file_verifies():
    assert protocol.verify() is True
    assert protocol._digest() == protocol.PROTOCOL_SHA


def test_tampered_protocol_fails_verification(tmp_path):
    # Editing any value must make verify() report the file as modified — the
    # whole anti-tuning guarantee. A pinned hash of a mutable file is worthless
    # if a one-character edit still verifies.
    edited = tmp_path / "protocol.yaml"
    text = protocol.PROTOCOL.read_text().replace("min_auc: 0.70", "min_auc: 0.50")
    edited.write_text(text)
    assert protocol.verify(edited) is False
    assert protocol._digest(edited) != protocol.PROTOCOL_SHA


def test_parser_handles_scalars_lists_and_comments():
    parsed = protocol._parse(
        "# header comment\n"
        "sec:\n"
        "  i: 26            # trailing comment stripped\n"
        "  f: 0.70\n"
        "  lst: [0.01, 0.05, 0.10]\n"
        "\n"                       # blank line ignored
    )
    assert parsed == {"sec": {"i": 26, "f": 0.70, "lst": [0.01, 0.05, 0.10]}}
    assert isinstance(parsed["sec"]["i"], int)
    assert isinstance(parsed["sec"]["f"], float)


def test_shell_emits_docking_vars(capsys):
    protocol._emit_shell("docking")
    out = capsys.readouterr().out.strip()
    # exactly the variables scripts/run_screen2.sh eval's, with frozen values
    assert out == "CX=70.61; CY=66.28; CZ=4.18; SIZE=26; EXHAUST=32; MODES=20; SEED=42"


# --- parser rejects malformed input --------------------------------------
# The parser is the sole reader of the frozen protocol. If it silently mis-reads
# a malformed file instead of refusing it, the tamper-detection above is moot: a
# broken protocol could still yield numbers the gate treats as authoritative.
# Every raise below is a "refuse, don't guess" guarantee, so each gets a test.

def test_parser_preserves_bare_string_values():
    # A value that is neither int, float, nor list stays a string verbatim
    # (e.g. a mode name), never coerced or dropped.
    parsed = protocol._parse("sec:\n  mode: axial\n")
    assert parsed == {"sec": {"mode": "axial"}}
    assert isinstance(parsed["sec"]["mode"], str)


def test_parser_rejects_top_level_line_without_colon():
    with pytest.raises(ValueError, match="expected 'section:' header"):
        protocol._parse("docking\n")


def test_parser_rejects_top_level_header_carrying_a_value():
    # A section header must be bare ("docking:"), not "docking: 26" -- otherwise a
    # value silently attaches to nothing.
    with pytest.raises(ValueError, match="expected 'section:' header"):
        protocol._parse("docking: 26\n")


def test_parser_rejects_value_before_any_section():
    # An indented key: value with no section above it has nowhere to live.
    with pytest.raises(ValueError, match="before any section"):
        protocol._parse("  center_x: 70.61\n")


def test_parser_rejects_indented_line_without_colon():
    with pytest.raises(ValueError, match="expected 'key: value'"):
        protocol._parse("docking:\n  center_x\n")


# --- the CLI dispatch (main) ---------------------------------------------

def test_main_show_prints_full_protocol_json(capsys):
    assert protocol.main([]) == 0                     # default subcommand is --show
    loaded = json.loads(capsys.readouterr().out)
    assert loaded == protocol.load()


def test_main_shell_emits_docking_vars(capsys):
    assert protocol.main(["--shell"]) == 0
    assert "EXHAUST=32" in capsys.readouterr().out


def test_main_freeze_prints_current_digest(capsys):
    assert protocol.main(["--freeze"]) == 0
    assert capsys.readouterr().out.strip() == protocol._digest()


def test_main_verify_returns_zero_for_unmodified_protocol(capsys):
    assert protocol.main(["--verify"]) == 0
    assert "verified unmodified" in capsys.readouterr().out
