"""Tests for the ERG11 re-caller ORCHESTRATION layer (scripts/recall_erg11.py).

The deterministic core (openafr/recaller.py) is well covered by test_recaller.py. The
messy half -- reads -> aligned -> consensus CDS, and the batch `fill` that writes the
`erg11_call` column back into a snapshot -- shells out to bioinformatics tools and so is
never exercised by the unit suite. But its *logic* (which commands get issued, how the
consensus FASTA is parsed, which rows `fill` selects and how it marks them) is pure Python
and testable offline by faking the tool calls.

Why this matters: `fill` runs AFTER multi-GB SRA downloads on a remote host. A bug in row
selection or column writing corrupts the snapshot at the most expensive possible moment --
exactly where you least want to discover it by running it for real. These tests pin that
logic without a single external tool or byte of network:

  * _read_consensus_fasta joins a wrapped record, and refuses empty / multi-record input,
  * reads_to_consensus issues the exact fetch/align/consensus pipeline (short-read preset,
    min-depth threaded through to `samtools consensus`) and fails loudly if no reads arrive,
  * cmd_fill re-calls only pending rows that carry a run_acc, picks the first linked run,
    writes call+source in place, leaves already-resolved rows untouched, and honors
    --limit / --dry-run.
"""
import argparse
import importlib.util
import pathlib
import types

import pytest

from openafr import recaller as R
from openafr import earlywarning as ew

_SCRIPT = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "recall_erg11.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("recall_erg11", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def script():
    return _load_script()


# --- _read_consensus_fasta: the bridge from tool output to the core ----------

def test_read_consensus_joins_wrapped_record(script, tmp_path):
    p = tmp_path / "c.fasta"
    p.write_text(">isolate erg11\nACGTACGT\nTTAACC\n  GGGG  \n")
    assert script._read_consensus_fasta(str(p)) == "ACGTACGTTTAACCGGGG"


def test_read_consensus_empty_file_is_refused(script, tmp_path):
    p = tmp_path / "empty.fasta"
    p.write_text("")
    with pytest.raises(SystemExit) as ei:
        script._read_consensus_fasta(str(p))
    assert "no sequence" in str(ei.value)


def test_read_consensus_header_only_is_refused(script, tmp_path):
    # A header with no bases must not read as an empty (silently wild-type) sequence.
    p = tmp_path / "header_only.fasta"
    p.write_text(">isolate\n")
    with pytest.raises(SystemExit) as ei:
        script._read_consensus_fasta(str(p))
    assert "no sequence" in str(ei.value)


def test_read_consensus_multi_record_is_refused(script, tmp_path):
    # `samtools consensus` on a single-contig reference should emit one record; more than
    # one means something is wrong upstream, and picking one arbitrarily would be a guess.
    p = tmp_path / "two.fasta"
    p.write_text(">a\nACGT\n>b\nTTTT\n")
    with pytest.raises(SystemExit) as ei:
        script._read_consensus_fasta(str(p))
    assert "found 2" in str(ei.value)


# --- reads_to_consensus: the exact tool pipeline -----------------------------

def test_reads_to_consensus_issues_expected_pipeline(script, tmp_path, monkeypatch):
    run_acc = "SRR9999999"
    ref_fasta = "/path/to/erg11_cds.fasta"
    calls = []

    def fake_run(cmd, **kw):
        calls.append(list(cmd))
        if cmd[0] == "fasterq-dump":
            outdir = pathlib.Path(cmd[cmd.index("-O") + 1])
            # paired reads, so the glob must pick up _1 and _2 in R1/R2 order
            (outdir / f"{run_acc}_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
            (outdir / f"{run_acc}_2.fastq").write_text("@r\nTGCA\n+\nIIII\n")
        return types.SimpleNamespace(returncode=0)

    monkeypatch.setattr(script, "_run", fake_run)
    cons = script.reads_to_consensus(run_acc, str(tmp_path), ref_fasta, min_depth=7)

    # fetch: split-spot, skip-technical, into the workdir, for this run
    fetch = calls[0]
    assert fetch[0] == "fasterq-dump"
    assert "--split-spot" in fetch and "--skip-technical" in fetch
    assert fetch[fetch.index("-O") + 1] == str(tmp_path)
    assert fetch[-1] == run_acc

    # align: short-read preset, reference, then both FASTQs in R1,R2 order
    mm = next(c for c in calls if c[0] == "minimap2")
    assert mm[:4] == ["minimap2", "-ax", "sr", ref_fasta]
    assert mm[4:] == [str(tmp_path / f"{run_acc}_1.fastq"),
                      str(tmp_path / f"{run_acc}_2.fastq")]

    # sort, index, then consensus with min-depth threaded through + 0.5 call fraction
    assert any(c[0] == "samtools" and c[1] == "sort" for c in calls)
    assert any(c[0] == "samtools" and c[1] == "index" for c in calls)
    consensus = next(c for c in calls if c[0] == "samtools" and c[1] == "consensus")
    assert consensus[consensus.index("--min-depth") + 1] == "7"
    assert consensus[consensus.index("--call-fract") + 1] == "0.5"

    assert cons == str(tmp_path / "erg11_consensus.fasta")


def test_reads_to_consensus_raises_when_no_reads_downloaded(script, tmp_path, monkeypatch):
    # fasterq-dump "succeeds" but yields no FASTQ (e.g. a withdrawn/empty run): must raise,
    # so _recall_one marks the isolate failed rather than aligning nothing into a call.
    def fake_run(cmd, **kw):
        return types.SimpleNamespace(returncode=0)   # never writes any FASTQ

    monkeypatch.setattr(script, "_run", fake_run)
    with pytest.raises(RuntimeError, match="no FASTQ"):
        script.reads_to_consensus("SRR0000000", str(tmp_path), "/ref.fasta")


# --- cmd_fill: batch row selection + in-place marking ------------------------

def _snapshot(path, rows):
    """Write a minimal but schema-valid snapshot from partial rows."""
    records = []
    for r in rows:
        rec = {c: "" for c in ew.SNAPSHOT_COLUMNS}
        rec.update(r)
        records.append(rec)
    ew.write_snapshot(str(path), records)


def _fill_args(snapshot, limit=0, dry_run=False):
    return argparse.Namespace(snapshot=str(snapshot), limit=limit, dry_run=dry_run)


@pytest.fixture
def patched_fill(script, monkeypatch):
    """cmd_fill with the tool gate and reference load stubbed, and a scripted
    _recall_one that records the run_acc it was handed."""
    monkeypatch.setattr(script, "_require_tools", lambda: None)
    monkeypatch.setattr(R, "load_reference", lambda *a, **k: ("CDS", "PROT"))
    seen = []

    def fake_recall_one(run_acc, reference, reference_fasta, keep_tmp=False):
        seen.append(run_acc)
        return {
            "SRR_CALLED": ("Y132F", "sra-recaller:called"),
            "SRR_WT": ("", "sra-recaller:wild-type"),
            "SRR_MULTI_A": ("K143R", "sra-recaller:called"),
        }[run_acc]

    monkeypatch.setattr(script, "_recall_one", fake_recall_one)
    return script, seen


def test_fill_recalls_pending_and_leaves_resolved_rows_alone(patched_fill, tmp_path):
    script, seen = patched_fill
    snap = tmp_path / "snap.tsv"
    _snapshot(snap, [
        {"isolate_key": "K1", "run_acc": "SRR_CALLED",
         "resistance_source": "pending:sra-recaller"},
        # multi-run cell: the first linked run must be the one re-called
        {"isolate_key": "K2", "run_acc": "SRR_MULTI_A,SRR_MULTI_B",
         "resistance_source": "pending:sra-recaller"},
        # pending but no run_acc -> nothing to re-call, must be skipped and left pending
        {"isolate_key": "K3", "run_acc": "",
         "resistance_source": "pending:sra-recaller"},
        # already resolved -> must not be touched
        {"isolate_key": "K4", "run_acc": "SRR_OLD", "erg11_call": "Y132F",
         "resistance_source": "sra-recaller:called"},
    ])

    rc = script.cmd_fill(_fill_args(snap))
    assert rc == 0
    assert seen == ["SRR_CALLED", "SRR_MULTI_A"]      # only pending+run_acc, first run

    by_key = {r["isolate_key"]: r for r in ew.read_snapshot(str(snap))}
    assert (by_key["K1"]["erg11_call"], by_key["K1"]["resistance_source"]) == \
        ("Y132F", "sra-recaller:called")
    assert (by_key["K2"]["erg11_call"], by_key["K2"]["resistance_source"]) == \
        ("K143R", "sra-recaller:called")
    # untouched rows keep their exact prior values
    assert by_key["K3"]["resistance_source"] == "pending:sra-recaller"
    assert by_key["K3"]["erg11_call"] == ""
    assert (by_key["K4"]["erg11_call"], by_key["K4"]["resistance_source"]) == \
        ("Y132F", "sra-recaller:called")


def test_fill_limit_stops_after_n(patched_fill, tmp_path):
    script, seen = patched_fill
    snap = tmp_path / "snap.tsv"
    _snapshot(snap, [
        {"isolate_key": "K1", "run_acc": "SRR_CALLED",
         "resistance_source": "pending:sra-recaller"},
        {"isolate_key": "K2", "run_acc": "SRR_WT",
         "resistance_source": "pending:sra-recaller"},
    ])

    script.cmd_fill(_fill_args(snap, limit=1))
    assert len(seen) == 1                              # only the first pending isolate


def test_fill_dry_run_recalls_but_does_not_write(patched_fill, tmp_path):
    script, seen = patched_fill
    snap = tmp_path / "snap.tsv"
    _snapshot(snap, [
        {"isolate_key": "K1", "run_acc": "SRR_CALLED",
         "resistance_source": "pending:sra-recaller"},
    ])
    before = snap.read_text()

    script.cmd_fill(_fill_args(snap, dry_run=True))
    assert seen == ["SRR_CALLED"]                      # the re-call still happened
    assert snap.read_text() == before                 # ...but the file is unchanged
