"""Tests for the lab-user genotype->verdict CLI (scripts/interpret_isolate.py, issue #163).

The CLI is thin glue over openafr.interpret + openafr.report (both unit-tested elsewhere), so
these lock the glue, not the engine:
  * each of the four verdict values is reachable from the command line and lands in the rendered
    brief (RESISTANCE_MARKER_DETECTED, UNCHARACTERIZED_VARIANT, NO_KNOWN_MARKER, UNRESOLVED);
  * one invocation can type an isolate for BOTH drug-classes;
  * --json emits a parseable report carrying the RUO disclaimer and summary buckets;
  * top-level flags (--json/--isolate) are honored even when written AFTER a --gene block
    (the parse_known_args split) -- the exact bug the CLI is built to avoid;
  * malformed/unsupported input exits non-zero with a message, never a traceback;
  * the consensus-FASTA path reads a file and routes through the real caller.
Offline, stdlib + pinned references only; no network, no external tools.
"""
import json

import pytest

from scripts import interpret_isolate as cli


def _run(argv):
    """Drive main() in-process; returns the exit code (0 on success)."""
    return cli.main(argv)


def test_resistance_marker_detected(capsys):
    assert _run(["--isolate", "B1", "--gene", "FKS1", "--variants", "S639P"]) == 0
    out = capsys.readouterr().out
    assert "RESISTANCE MARKER" in out and "S639P" in out
    assert "RESEARCH USE ONLY" in out          # RUO banner is never dropped


def test_uncharacterized_variant(capsys):
    # A real change off the resistance panel -> the explicit "I don't know" verdict.
    assert _run(["--isolate", "B2", "--gene", "FKS1", "--variants", "A100T"]) == 0
    assert "UNCHARACTERIZED" in capsys.readouterr().out


def test_no_known_marker_from_empty_variants(capsys):
    # An empty variant list is a valid "nothing found" input, not missing input.
    assert _run(["--isolate", "B3", "--gene", "FKS1", "--variants"]) == 0
    assert "NO KNOWN MARKER" in capsys.readouterr().out


def test_uncalled_panel_residue_forces_unresolved(capsys):
    # A panel residue declared uncalled can't read as wild-type (contract rule 4).
    assert _run(["--isolate", "B4", "--gene", "FKS1", "--variants", "--uncalled", "639"]) == 0
    assert "UNRESOLVED" in capsys.readouterr().out


def test_both_drug_classes_in_one_call(capsys):
    assert _run(["--isolate", "B5",
                 "--gene", "FKS1", "--variants", "S639P",
                 "--gene", "ERG11", "--variants", "Y132F"]) == 0
    out = capsys.readouterr().out
    assert "(FKS1)" in out and "(ERG11)" in out


def test_json_flag_after_gene_block(capsys):
    # The parse_known_args split: --json written last must still be a top-level flag, not
    # swallowed by the preceding --gene block.
    assert _run(["--isolate", "B6", "--gene", "FKS1", "--variants", "S639P", "--json"]) == 0
    doc = json.loads(capsys.readouterr().out)
    assert "RESEARCH USE ONLY" in doc["disclaimer"]
    assert "summary" in doc
    assert doc["verdicts"][0]["verdict"] == "RESISTANCE_MARKER_DETECTED"


def test_consensus_fasta_path(tmp_path, capsys):
    """A single-record consensus FASTA routes through the real FKS1 caller. The pinned
    wild-type CDS must call clean (NO_KNOWN_MARKER), proving the file->caller bridge works."""
    import pathlib
    root = pathlib.Path(cli.__file__).resolve().parent.parent
    wt = (root / "data/earlywarning/demo/fks1_wildtype.fasta").read_text()
    f = tmp_path / "cds.fasta"
    f.write_text(wt)
    assert _run(["--isolate", "B7", "--gene", "FKS1", "--cds", str(f)]) == 0
    assert "NO KNOWN MARKER" in capsys.readouterr().out


def test_unsupported_gene_exits(capsys):
    with pytest.raises(SystemExit) as e:
        _run(["--isolate", "B8", "--gene", "CYP51A", "--variants", "Y132F"])
    assert e.value.code != 0


def test_malformed_variant_token_exits():
    with pytest.raises(SystemExit) as e:
        _run(["--isolate", "B9", "--gene", "FKS1", "--variants", "not-a-token"])
    assert e.value.code != 0


def test_no_gene_block_exits():
    with pytest.raises(SystemExit) as e:
        _run(["--isolate", "B10"])
    assert e.value.code != 0
