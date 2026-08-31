"""Tests for the FKS1 re-caller ORCHESTRATION layer (scripts/recall_fks1.py).

The deterministic core (openafr/fks1_caller.py) is well covered by test_fks1_caller.py. The
messy half -- reads -> aligned -> PER-WINDOW consensus, and the deterministic `call` that
turns a consensus into (fks1_call, source) -- either shells out to bioinformatics tools or is
pure Python. These tests pin the pure-Python logic without a single external tool or byte of
network, and fake the tool calls to assert the exact windowed pipeline gets issued.

Why this matters: `recall` runs AFTER multi-GB SRA downloads on a remote host. The windowed
departure from ERG11 -- one `samtools consensus -r <region>` per hot-spot -- is the one novel
piece, and a wrong region string or a dropped window would mis-call an isolate at the most
expensive possible moment. These tests lock:

  * _read_consensus_fasta joins a wrapped record and refuses empty / multi-record input,
  * _reference_contig / _window_region compute the BAM contig name and 1-based nt regions
    from the core's residue numbering (so they can never drift from fks1_caller),
  * reads_to_window_consensus issues fetch/align then ONE region-scoped consensus per window
    (short-read preset, -a, min-depth + call-fract threaded through) and fails loudly with no
    reads,
  * _recall_one maps a windowed result to (call, source) and marks tool failures failed(...),
  * cmd_call calls from a full CDS consensus and from per-window FASTAs, and refuses a
    full consensus that cannot be anchored.
"""
import argparse
import importlib.util
import pathlib
import types

import pytest

from openafr import fks1_caller as F
from openafr import runlog

_SCRIPT = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "recall_fks1.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("recall_fks1", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def script():
    return _load_script()


# --- helpers to build reference-anchored consensus with a known substitution --

_F_CODON = "TTT"   # translates to F -- the S639F panel mutant used across these tests


def _ref_cds():
    return F.load_reference()[0]


def _mutate_residue(cds, residue_pos, codon):
    """Return `cds` with the codon at 1-based `residue_pos` replaced (for building a known
    substitution consensus straight off the pinned reference)."""
    i = (residue_pos - 1) * 3
    assert len(codon) == 3
    return cds[:i] + codon + cds[i + 3:]


# --- _read_consensus_fasta: the bridge from tool output to the core ----------

def test_read_consensus_joins_wrapped_record(script, tmp_path):
    p = tmp_path / "c.fasta"
    p.write_text(">isolate hs1\nACGTACGT\nTTAACC\n  GGGG  \n")
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
    # A region-scoped `samtools consensus` should emit one record; more than one means
    # something is wrong upstream, and picking one arbitrarily would be a guess.
    p = tmp_path / "two.fasta"
    p.write_text(">a\nACGT\n>b\nTTTT\n")
    with pytest.raises(SystemExit) as ei:
        script._read_consensus_fasta(str(p))
    assert "found 2" in str(ei.value)


# --- contig + region math: derived from the core, never hand-typed -----------

def test_reference_contig_is_header_first_token(script, tmp_path):
    p = tmp_path / "ref.fasta"
    p.write_text(">XM_085597048.1 Candidozyma auris B8441 GSC1/FKS1 CDS\nACGTACGT\n")
    assert script._reference_contig(str(p)) == "XM_085597048.1"


def test_reference_contig_uses_the_real_reference(script):
    # The header pinned in the repo -> the contig minimap2 emits into the BAM.
    assert script._reference_contig(F.DEFAULT_REFERENCE) == "XM_085597048.1"


def test_window_region_matches_core_numbering(script):
    # residue r -> nt (r-1)*3+1 .. r*3 (1-based inclusive), and the span is nt_len long.
    hs1 = F.FKS1_WINDOWS["HS1"]
    hs2 = F.FKS1_WINDOWS["HS2"]
    assert script._window_region("REF", hs1) == "REF:1903-1929"   # residues 635-643
    assert script._window_region("REF", hs2) == "REF:4048-4074"   # residues 1350-1358
    for w in (hs1, hs2):
        start, end = script._window_region("REF", w).split(":")[1].split("-")
        assert int(end) - int(start) + 1 == w.nt_len


# --- reads_to_window_consensus: the windowed tool pipeline -------------------

def test_reads_to_window_consensus_issues_per_window_pipeline(script, tmp_path, monkeypatch):
    run_acc = "SRR9999999"
    ref_fasta = "/path/to/fks1_cds.fasta"
    calls = []

    monkeypatch.setattr(script, "_reference_contig", lambda p: "XM_085597048.1")

    def fake_run(cmd, **kw):
        calls.append(list(cmd))
        if cmd[0] == "fasterq-dump":
            outdir = pathlib.Path(cmd[cmd.index("-O") + 1])
            (outdir / f"{run_acc}_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
            (outdir / f"{run_acc}_2.fastq").write_text("@r\nTGCA\n+\nIIII\n")
        elif cmd[0] == "samtools" and cmd[1] == "consensus":
            # write the region-named record the consensus step is expected to emit
            region = cmd[cmd.index("-r") + 1]
            kw["stdout"].write(f">{region}\nACGACGACGACGACGACGACGACGACG\n")
        return types.SimpleNamespace(returncode=0)

    monkeypatch.setattr(script, "_run", fake_run)
    window_seqs = script.reads_to_window_consensus(run_acc, str(tmp_path), ref_fasta,
                                                   min_depth=7)

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
    assert any(c[0] == "samtools" and c[1] == "sort" for c in calls)
    assert any(c[0] == "samtools" and c[1] == "index" for c in calls)

    # ONE region-scoped consensus per window, min-depth + call-fract threaded through, -a set
    cons = [c for c in calls if c[0] == "samtools" and c[1] == "consensus"]
    assert len(cons) == len(F.FKS1_WINDOWS)
    regions = {c[c.index("-r") + 1] for c in cons}
    assert regions == {"XM_085597048.1:1903-1929", "XM_085597048.1:4048-4074"}
    for c in cons:
        assert "-a" in c
        assert c[c.index("--min-depth") + 1] == "7"
        assert c[c.index("--call-fract") + 1] == "0.5"

    # returns a dict keyed by window name, each the parsed consensus string
    assert set(window_seqs) == set(F.FKS1_WINDOWS)
    assert all(seq == "ACGACGACGACGACGACGACGACGACG" for seq in window_seqs.values())


def test_reads_to_window_consensus_raises_when_no_reads(script, tmp_path, monkeypatch):
    # fasterq-dump "succeeds" but yields no FASTQ: must raise so _recall_one marks the
    # isolate failed rather than aligning nothing into a call.
    monkeypatch.setattr(script, "_reference_contig", lambda p: "REF")
    monkeypatch.setattr(script, "_run",
                        lambda cmd, **kw: types.SimpleNamespace(returncode=0))
    with pytest.raises(RuntimeError, match="no FASTQ"):
        script.reads_to_window_consensus("SRR0000000", str(tmp_path), "/ref.fasta")


# --- _recall_one: windowed result -> (call, source), failures marked ---------

def test_recall_one_maps_window_seqs_to_call(script, monkeypatch):
    reference = F.load_reference()
    ref_cds = reference[0]
    # HS1 window nt with residue 639 mutated S->F, HS2 window wild-type off the reference.
    mutated = _mutate_residue(ref_cds, 639, _F_CODON)
    window_seqs = {
        "HS1": F.FKS1_WINDOWS["HS1"].nt_slice(mutated),
        "HS2": F.FKS1_WINDOWS["HS2"].nt_slice(ref_cds),
    }
    monkeypatch.setattr(script, "reads_to_window_consensus",
                        lambda *a, **k: window_seqs)
    call, source = script._recall_one("SRR1", reference, F.DEFAULT_REFERENCE)
    assert call == "S639F"
    assert source == "sra-fks1-recaller:HS1=called,HS2=wild-type"


def test_recall_one_marks_tool_failure_failed(script, monkeypatch):
    reference = F.load_reference()

    def boom(*a, **k):
        raise RuntimeError("fasterq-dump produced no FASTQ for SRRX")

    monkeypatch.setattr(script, "reads_to_window_consensus", boom)
    call, source = script._recall_one("SRRX", reference, F.DEFAULT_REFERENCE)
    assert call == ""
    assert source.startswith("sra-fks1-recaller:failed(")
    assert "no FASTQ" in source


# --- cmd_call: deterministic call from full CDS and from per-window FASTAs ----

def _call_args(consensus_fasta=None, window=None):
    return argparse.Namespace(consensus_fasta=consensus_fasta, window=window)


def test_cmd_call_from_full_cds(script, tmp_path, capsys):
    ref_cds = _ref_cds()
    mutated = _mutate_residue(ref_cds, 639, _F_CODON)
    p = tmp_path / "cons.fasta"
    p.write_text(">isolate\n" + mutated + "\n")
    rc = script.cmd_call(_call_args(consensus_fasta=str(p)))
    assert rc == 0
    out = capsys.readouterr().out
    assert "call:   S639F" in out
    assert "sra-fks1-recaller:HS1=called,HS2=wild-type" in out


def test_cmd_call_from_per_window_fasta(script, tmp_path, capsys):
    ref_cds = _ref_cds()
    hs1_nt = F.FKS1_WINDOWS["HS1"].nt_slice(_mutate_residue(ref_cds, 639, _F_CODON))
    p = tmp_path / "hs1.fasta"
    p.write_text(">hs1\n" + hs1_nt + "\n")
    # only HS1 supplied -> HS2 reported missing, never assumed wild-type
    rc = script.cmd_call(_call_args(window=[f"HS1={p}"]))
    assert rc == 0
    out = capsys.readouterr().out
    assert "call:   S639F" in out
    assert "HS2=missing" in out


def test_cmd_call_unanchorable_full_consensus_is_refused(script, tmp_path, capsys):
    # A consensus that is not the reference length cannot anchor the windows -> refuse (2),
    # never guess.
    p = tmp_path / "short.fasta"
    p.write_text(">iso\nACGTACGTACGT\n")
    rc = script.cmd_call(_call_args(consensus_fasta=str(p)))
    assert rc == 2
    assert "refused" in capsys.readouterr().out


def test_cmd_call_unknown_window_name_is_refused(script, tmp_path):
    p = tmp_path / "w.fasta"
    p.write_text(">w\nACGACGACG\n")
    with pytest.raises(SystemExit) as ei:
        script.cmd_call(_call_args(window=[f"HSX={p}"]))
    assert "unknown window" in str(ei.value)


# --- cmd_recall: runlog-wrapped single-isolate orchestration -----------------

def test_cmd_recall_logs_and_reports(script, monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(script, "_require_tools", lambda: None)
    monkeypatch.setattr(F, "load_reference", lambda *a, **k: ("CDS", "PROT"))
    monkeypatch.setattr(runlog, "DEFAULT_LOG_DIR", tmp_path / "runlog")
    monkeypatch.setattr(script, "_recall_one",
                        lambda *a, **k: ("S639F", "sra-fks1-recaller:HS1=called,HS2=wild-type"))
    args = argparse.Namespace(run_acc="SRR42", keep_tmp=False)
    rc = script.cmd_recall(args)
    assert rc == 0
    out = capsys.readouterr().out
    assert "S639F" in out
    # the run was logged
    logged = (tmp_path / "runlog" / "recall-fks1.jsonl").read_text()
    assert "SRR42" in logged and "S639F" in logged


def test_cmd_recall_returns_nonzero_on_failure(script, monkeypatch, tmp_path):
    monkeypatch.setattr(script, "_require_tools", lambda: None)
    monkeypatch.setattr(F, "load_reference", lambda *a, **k: ("CDS", "PROT"))
    monkeypatch.setattr(runlog, "DEFAULT_LOG_DIR", tmp_path / "runlog")
    monkeypatch.setattr(script, "_recall_one",
                        lambda *a, **k: ("", "sra-fks1-recaller:failed(no FASTQ)"))
    args = argparse.Namespace(run_acc="SRRbad", keep_tmp=False)
    assert script.cmd_recall(args) == 1


# --- batch surface (T3): plan / fill / prevalence / migrate over a snapshot ----

from openafr import earlywarning as _ew   # noqa: E402


def _snapshot(tmp_path, rows):
    """Write a current-schema snapshot from (target_acc, run) pairs; return its path."""
    raw = [{"target_acc": t, "Run": run, "geo_loc_name": "Pakistan",
            "target_creation_date": "2024-01-01", "collection_date": "2023",
            "biosample_acc": "SAMN1", "asm_acc": "NULL", "lat_lon": "NULL"}
           for t, run in rows]
    records = _ew.normalize(raw)[0]
    p = tmp_path / "snap.tsv"
    _ew.write_snapshot(p, records)
    return p


def test_cmd_plan_matches_fill_selection_without_tools(script, tmp_path, capsys):
    # Two isolates with reads (fill-eligible) + one without (skipped, stays pending).
    p = _snapshot(tmp_path, [("PDT1.1", "SRR1"), ("PDT2.1", "SRR2"), ("PDT3.1", "")])
    args = argparse.Namespace(snapshot=str(p), limit=0, random=False, seed=0,
                              list_runs=True)
    assert script.cmd_plan(args) == 0
    out = capsys.readouterr().out
    # PDT3 has no reads -> source is pending:no-reads, not the fill-eligible
    # pending:sra-fks1-recaller state, so it is not in this count (as in the ERG11 plan).
    assert "pending:sra-fks1-recaller 2" in out
    assert "fill would re-call:      2 isolate(s)" in out
    assert "SRR1" in out and "SRR2" in out


def test_cmd_prevalence_is_not_measured_on_fresh_snapshot(script, tmp_path, capsys):
    p = _snapshot(tmp_path, [("PDT1.1", "SRR1"), ("PDT2.1", "SRR2")])
    args = argparse.Namespace(snapshot=str(p))
    assert script.cmd_prevalence(args) == 0
    out = capsys.readouterr().out
    assert "NOT MEASURED YET" in out


def test_cmd_prevalence_reports_frequency_after_a_fill(script, tmp_path, capsys):
    p = _snapshot(tmp_path, [("PDT1.1", "SRR1"), ("PDT2.1", "SRR2")])
    records = _ew.read_snapshot(p)
    records[0]["fks1_call"] = "S639F"
    records[0]["fks1_resistance_source"] = "sra-fks1-recaller:HS1=called,HS2=wild-type"
    records[1]["fks1_call"] = ""
    records[1]["fks1_resistance_source"] = "sra-fks1-recaller:HS1=wild-type,HS2=wild-type"
    _ew.write_snapshot(p, records)
    assert script.cmd_prevalence(argparse.Namespace(snapshot=str(p))) == 0
    out = capsys.readouterr().out
    assert "resolved panel:    2" in out
    assert "event frequency:   50.0%" in out
    assert "S639F x1" in out


def test_cmd_fill_marks_calls_and_logs(script, monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(script, "_require_tools", lambda: None)
    monkeypatch.setattr(F, "load_reference", lambda *a, **k: ("CDS", "PROT"))
    monkeypatch.setattr(runlog, "DEFAULT_LOG_DIR", tmp_path / "runlog")
    monkeypatch.setattr(script, "_recall_one",
                        lambda *a, **k: ("S639F",
                                         "sra-fks1-recaller:HS1=called,HS2=wild-type"))
    p = _snapshot(tmp_path, [("PDT1.1", "SRR1"), ("PDT2.1", "")])
    args = argparse.Namespace(snapshot=str(p), limit=0, random=False, seed=0,
                              dry_run=False)
    assert script.cmd_fill(args) == 0
    records = {r["isolate_key"]: r for r in _ew.read_snapshot(p)}
    # the with-reads isolate got the call; the no-reads isolate stayed pending, not faked
    assert records["PDT1"]["fks1_call"] == "S639F"
    assert records["PDT2"]["fks1_call"] == "" \
        and records["PDT2"]["fks1_resistance_source"] == "pending:no-reads"
    assert "S639F" in (tmp_path / "runlog" / "fill-fks1.jsonl").read_text()


def test_cmd_migrate_upgrades_a_pre_v2_snapshot(script, tmp_path, capsys):
    p = tmp_path / "old.tsv"
    p.write_text(
        "isolate_key\ttarget_acc\trun_acc\terg11_call\tresistance_source\n"
        "PDT1\tPDT1.1\tSRR9\t\tpending:sra-recaller\n"
    )
    assert script.cmd_migrate(argparse.Namespace(snapshot=str(p))) == 0
    assert "migrated" in capsys.readouterr().out
    assert "fks1_resistance_source" in p.read_text()
    # idempotent second run
    assert script.cmd_migrate(argparse.Namespace(snapshot=str(p))) == 0
    assert "already current" in capsys.readouterr().out
