"""Tests for the frozen run protocol (openafr/protocol.py, task T2).

The protocol carries the numbers most tempting to silently tune into a pass — the
docking box, the search effort, and the gate's pass bar — so these lock both the
values and the tamper-detection that guards them.
"""
import pathlib

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
