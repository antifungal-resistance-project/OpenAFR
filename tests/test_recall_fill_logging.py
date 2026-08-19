"""Integration test for the `fill` path's run-log wiring (scripts/recall_erg11.py).

The runlog unit tests cover the logger; this covers the *wiring* the logger can't see:
that a `fill` run records the before/after snapshot digest, the per-isolate transcript,
and the called/partial/failed summary -- the fields that let a re-`fill` be audited
instead of silently overwriting its predecessor. The real reads->consensus orchestration
(tools + network) is stubbed, so no bio tools or network are touched.
"""
import pathlib
import sys
import types

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / "scripts"))

from openafr import earlywarning as ew  # noqa: E402
from openafr import runlog              # noqa: E402
import recall_erg11 as rc               # noqa: E402


def _rec(key, run_acc, source="pending:sra-recaller"):
    return {c: "" for c in ew.SNAPSHOT_COLUMNS} | {
        "isolate_key": key, "target_acc": f"{key}.1", "run_acc": run_acc,
        "erg11_call": "", "resistance_source": source,
    }


def _stub_tools(monkeypatch, tmp_path, outcomes):
    """Wire cmd_fill to run without bio tools/network. `outcomes` maps run_acc ->
    (call, source); an unlisted run_acc fails honestly."""
    monkeypatch.setattr(rc, "_require_tools", lambda: None)
    monkeypatch.setattr(rc.R, "load_reference", lambda: "REF")
    monkeypatch.setattr(rc, "_recall_one",
                        lambda run_acc, ref, reff, keep_tmp=False:
                        outcomes.get(run_acc, ("", "sra-recaller:failed(no reads)")))
    monkeypatch.setattr(runlog, "DEFAULT_LOG_DIR", tmp_path / "runlog")


def test_fill_logs_run_with_before_after_digest(tmp_path, monkeypatch):
    snap = tmp_path / "PDG_test.tsv"
    ew.write_snapshot(snap, [
        _rec("PDT1", "SRR1"),
        _rec("PDT2", "SRR2"),
        _rec("PDT3", "SRR3", source="pending:no-reads"),  # not pending:sra-recaller -> skipped
    ])
    _stub_tools(monkeypatch, tmp_path,
                {"SRR1": ("Y132F", "sra-recaller:panel(Y132F)")})  # PDT2 (SRR2) fails

    args = types.SimpleNamespace(snapshot=str(snap), limit=0, dry_run=False,
                                 random=False, seed=0)
    assert rc.cmd_fill(args) == 0

    runs = runlog.read_runs("fill", log_dir=tmp_path / "runlog")
    assert len(runs) == 1
    e = runs[0]
    assert e["kind"] == "fill" and e["status"] == "ok"
    assert e["n_pending"] == 2                       # PDT3 (no-reads) is not attempted
    assert e["wrote"] is True
    assert e["snapshot_sha256_before"] != e["snapshot_sha256_after"]  # PDT1 got a call
    assert e["snapshot"] == "PDG_test.tsv"
    assert e["summary"] == {"called": 1, "partial": 0, "failed": 1, "attempted": 2}
    assert {it["isolate_key"] for it in e["items"]} == {"PDT1", "PDT2"}
    # the snapshot itself now carries PDT1's call -- the log and the file agree
    on_disk = {r["isolate_key"]: r["erg11_call"] for r in ew.read_snapshot(snap)}
    assert on_disk["PDT1"] == "Y132F" and on_disk["PDT2"] == ""


def test_dry_run_records_no_write_but_still_logs(tmp_path, monkeypatch):
    snap = tmp_path / "PDG_dry.tsv"
    before = ew.write_snapshot(snap, [_rec("PDT9", "SRR9")])
    _stub_tools(monkeypatch, tmp_path,
                {"SRR9": ("Y132F", "sra-recaller:panel(Y132F)")})

    args = types.SimpleNamespace(snapshot=str(snap), limit=0, dry_run=True,
                                 random=False, seed=0)
    assert rc.cmd_fill(args) == 0

    # the run is logged, but the snapshot on disk is untouched
    assert ew.snapshot_sha256(ew.read_snapshot(snap)) == before
    e = runlog.read_runs("fill", log_dir=tmp_path / "runlog")[0]
    assert e["wrote"] is False
    assert e["dry_run"] is True
    # a dry run computes what it *would* write, so before/after digests differ...
    assert e["snapshot_sha256_after"] != e["snapshot_sha256_before"]
    # ...but nothing hit disk
