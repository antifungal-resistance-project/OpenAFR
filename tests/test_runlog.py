"""Tests for the append-only run log (openafr/runlog.py).

The log is the durable memory of the re-caller's un-reproducible tool+network runs, so
these lock the guarantees a reader leans on:
  * a run writes exactly one append-only line carrying its provenance + outcome,
  * items and the run-level summary round-trip,
  * a crashed run body is still logged (status='error') and the exception re-raised --
    the whole point is that a failed run is recorded, not lost,
  * malformed lines never break read-back,
  * provenance probes are best-effort (absent tool -> None, never a raise),
  * logging never touches the network or shells out to the real bio tools.
"""
import hashlib
import json

import pytest

from openafr import runlog


def _lines(path):
    return [json.loads(l) for l in path.read_text().splitlines() if l.strip()]


def test_run_writes_one_append_only_line_with_provenance(tmp_path):
    ref = tmp_path / "erg11_ref.fasta"
    ref.write_text(">ref\nATGAAA\n")
    with runlog.record_run("sanity", reference_path=ref, tools=[],
                           log_dir=tmp_path, extra={"fixture": "known.tsv"}) as run:
        run.add_item(strain="B8441", outcome="PASS")
        run.set(status="ok", summary={"pass": 1, "fail": 0})

    log = tmp_path / "sanity.jsonl"
    rows = _lines(log)
    assert len(rows) == 1
    e = rows[0]
    assert e["kind"] == "sanity"
    assert e["status"] == "ok"
    assert e["schema"] == runlog.SCHEMA
    assert e["fixture"] == "known.tsv"
    assert e["summary"] == {"pass": 1, "fail": 0}
    assert e["items"] == [{"strain": "B8441", "outcome": "PASS"}]
    # reference digest is the sha256 of the file bytes, so a call is tied to what it called against
    assert e["reference"]["name"] == "erg11_ref.fasta"
    assert e["reference"]["sha256"] == hashlib.sha256(ref.read_bytes()).hexdigest()
    assert "started_at" in e and "finished_at" in e and len(e["run_id"]) == 12


def test_two_runs_append_not_overwrite(tmp_path):
    for i in range(2):
        with runlog.record_run("fill", log_dir=tmp_path) as run:
            run.set(summary={"attempted": i})
    rows = _lines(tmp_path / "fill.jsonl")
    assert [r["summary"]["attempted"] for r in rows] == [0, 1]
    # distinct runs get distinct ids
    assert rows[0]["run_id"] != rows[1]["run_id"]


def test_crashed_run_is_logged_and_reraises(tmp_path):
    with pytest.raises(ValueError, match="boom"):
        with runlog.record_run("fill", log_dir=tmp_path) as run:
            run.add_item(isolate_key="PDT1", source="pending")
            raise ValueError("boom mid-batch")

    rows = _lines(tmp_path / "fill.jsonl")
    assert len(rows) == 1
    assert rows[0]["status"] == "error"
    assert "boom" in rows[0]["error"]
    # whatever accumulated before the crash is preserved -- that's the audit value
    assert rows[0]["items"] == [{"isolate_key": "PDT1", "source": "pending"}]


def test_read_runs_oldest_first_and_skips_malformed(tmp_path):
    with runlog.record_run("sanity", log_dir=tmp_path) as run:
        run.set(summary={"pass": 1})
    # a corrupt half-written line must not sink the whole read
    with open(tmp_path / "sanity.jsonl", "a") as fh:
        fh.write("{not valid json\n")
    with runlog.record_run("sanity", log_dir=tmp_path) as run:
        run.set(summary={"pass": 2})

    runs = runlog.read_runs("sanity", log_dir=tmp_path)
    assert [r["summary"]["pass"] for r in runs] == [1, 2]


def test_read_runs_missing_file_is_empty(tmp_path):
    assert runlog.read_runs("recall", log_dir=tmp_path) == []


def test_append_failure_is_swallowed(tmp_path, capsys):
    # point the log at a path whose parent is a file, so mkdir/open cannot succeed;
    # the run must still complete (best-effort logging never breaks the science run).
    blocker = tmp_path / "blocker"
    blocker.write_text("i am a file, not a dir")
    with runlog.record_run("sanity", log_dir=blocker / "subdir") as run:
        run.set(status="ok")
    assert "could not append" in capsys.readouterr().err


def test_file_digest_and_missing(tmp_path):
    f = tmp_path / "x"
    f.write_bytes(b"abc")
    assert runlog.file_digest(f) == hashlib.sha256(b"abc").hexdigest()
    assert runlog.file_digest(tmp_path / "nope") is None


def test_tool_version_absent_tool_is_none():
    # best-effort provenance: a tool that isn't installed yields None, never a raise
    assert runlog.tool_version("definitely-not-a-real-binary-xyz") is None
    versions = runlog.tool_versions(["definitely-not-a-real-binary-xyz"])
    assert versions == {"definitely-not-a-real-binary-xyz": None}
