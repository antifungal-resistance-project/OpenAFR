"""Tests for the non-billable FKS1 fill-due watchdog (openafr.fill_status + fill_due_check).

The resistance-weather page overlays a committed, billable FKS1 call table onto a live NCBI pull;
nothing re-runs that fill on a cadence, so the page drifts stale. `fill_status.fill_due` measures
the in-window re-call gap from a FRESH pull against the committed table -- the exact delta a
refill would add -- so a scheduled watchdog can ring the bell only when the paid step pays off.

These lock the delta arithmetic (uncalled_new == new in-window pending-with-run after overlay),
the two due-triggers (threshold crossed, or the last fill's whole window rolled off the page),
the quiet/current no-op, and the never-filled day-0 case. Pure functions -- no network/tools.
"""
import importlib.util
import json
import pathlib
import subprocess
import sys

from openafr import earlywarning as ew
from openafr import fill_status

_SCRIPT = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "fill_due_check.py"
ROOT = pathlib.Path(__file__).resolve().parent.parent


def _rec(key, source, call="", run="SRR1", created="2025-02-01"):
    return {"isolate_key": key, "target_acc": key + ".1", "run_acc": run,
            "target_creation_date": created, "fks1_call": call,
            "fks1_resistance_source": source}


def _write_committed(calls_dir, resolved_records, *, release="PDG000000067.700",
                     as_of="2025-01-01", window_days=180):
    """Commit a call table + INDEX row for `release` so fill_due overlays it, mirroring a real
    fill's durable output (earlywarning.write_fks1_call_table + update_fks1_calls_index)."""
    calls_dir.mkdir(parents=True, exist_ok=True)
    table = calls_dir / f"{release}.tsv"
    sha, n = ew.write_fks1_call_table(str(table), resolved_records)
    ew.update_fks1_calls_index(str(calls_dir / "INDEX.tsv"), {
        "release_tag": release, "as_of": as_of, "window_days": window_days,
        "n_called": n, "sha256": sha, "written_at_utc": ew.utcnow_iso(),
    })


# --- the incremental delta: uncalled_new == new in-window pending-with-run after overlay ------

def test_uncalled_new_is_the_post_overlay_in_window_pending(tmp_path):
    # A fresh pull: 5 in-window pending-with-run isolates (A..E). The committed table resolved
    # two of them (A, B), so a refill would add the remaining 3.
    pull = [_rec(k, "pending:sra-fks1-recaller") for k in ("A", "B", "C", "D", "E")]
    committed = [_rec("A", "sra-fks1-recaller:called", call="S639F"),
                 _rec("B", "sra-fks1-recaller:called")]  # re-called wild-type still resolved
    _write_committed(tmp_path, committed)

    status = fill_status.fill_due(pull, str(tmp_path), as_of="2025-01-01",
                                  window_days=180, threshold=25)
    assert status["total_in_window"] == 5
    assert status["resolved_in_window"] == 2      # A, B overlaid as resolved
    assert status["uncalled_new"] == 3            # C, D, E still pending-with-run
    assert status["newest_release"] == "PDG000000067.700"


def test_pending_without_run_is_not_counted_as_refill_workload(tmp_path):
    # An in-window pending isolate with no SRA run can't be re-called, so it is NOT refill work.
    pull = [_rec("A", "pending:sra-fks1-recaller", run=""),
            _rec("B", "pending:sra-fks1-recaller", run="SRR9")]
    _write_committed(tmp_path, [])
    status = fill_status.fill_due(pull, str(tmp_path), as_of="2025-01-01", threshold=25)
    assert status["uncalled_new"] == 1            # only B (A has no run_acc)


# --- due triggers -----------------------------------------------------------------------------

def test_due_when_uncalled_new_crosses_threshold(tmp_path):
    pull = [_rec(f"K{i}", "pending:sra-fks1-recaller", run=f"SRR{i}") for i in range(30)]
    _write_committed(tmp_path, [])
    due = fill_status.fill_due(pull, str(tmp_path), as_of="2025-01-01", threshold=25)
    assert due["due"] is True and due["uncalled_new"] == 30 and due["window_aged_out"] is False
    quiet = fill_status.fill_due(
        [_rec(f"K{i}", "pending:sra-fks1-recaller", run=f"SRR{i}") for i in range(10)],
        str(tmp_path), as_of="2025-01-01", threshold=25)
    assert quiet["due"] is False and quiet["uncalled_new"] == 10


def test_due_when_last_fill_window_rolled_off_the_page(tmp_path):
    # Committed fill is a year stale: its window (2024-01-01 +180d) ends before the current
    # window starts (2025-01-01), so NONE of its calls land on a displayed isolate -> aged out,
    # due even though almost nothing new is uncalled.
    pull = [_rec("A", "pending:sra-fks1-recaller", run="SRR1")]
    _write_committed(tmp_path, [_rec("Z", "sra-fks1-recaller:called")],
                     as_of="2024-01-01", window_days=180)
    status = fill_status.fill_due(pull, str(tmp_path), as_of="2025-01-01", threshold=25)
    assert status["window_aged_out"] is True and status["due"] is True


def test_quiet_when_window_current_and_delta_small(tmp_path):
    # Committed window still covers the displayed window and only one new isolate is uncalled.
    pull = [_rec("A", "pending:sra-fks1-recaller", run="SRR1"),
            _rec("B", "pending:sra-fks1-recaller", run="SRR2")]
    _write_committed(tmp_path, [_rec("A", "sra-fks1-recaller:called")],
                     as_of="2025-01-01", window_days=180)
    status = fill_status.fill_due(pull, str(tmp_path), as_of="2025-01-01", threshold=25)
    assert status["due"] is False and status["window_aged_out"] is False
    assert status["uncalled_new"] == 1


def test_due_when_never_filled(tmp_path):
    # Empty calls dir -> no committed table/INDEX: the page is at day-0, always due.
    (tmp_path).mkdir(parents=True, exist_ok=True)
    pull = [_rec("A", "pending:sra-fks1-recaller", run="SRR1")]
    status = fill_status.fill_due(pull, str(tmp_path), as_of="2025-01-01", threshold=25)
    assert status["due"] is True and status["newest_release"] is None
    assert "no committed" in status["reason"].lower()


# --- the CLI handshake: issue body on due, github-output flags --------------------------------

def test_cli_writes_issue_and_outputs_when_due(tmp_path):
    snap = tmp_path / "PDG000000067.800.tsv"
    records = [_rec(f"K{i}", "pending:sra-fks1-recaller", run=f"SRR{i}") for i in range(30)]
    _write_snapshot(snap, records)

    calls_dir = tmp_path / "calls"
    calls_dir.mkdir()
    issue_out = tmp_path / "issue.md"
    gh_out = tmp_path / "gh.txt"
    r = subprocess.run(
        [sys.executable, str(_SCRIPT), "--snapshot", str(snap), "--calls-dir", str(calls_dir),
         "--as-of", "2025-01-01", "--window", "180", "--threshold", "25",
         "--issue-out", str(issue_out), "--github-output", str(gh_out)],
        capture_output=True, text=True, cwd=str(ROOT))
    assert r.returncode == 0, r.stderr
    out = gh_out.read_text()
    assert "due=true" in out and "uncalled=30" in out and "label=fks1-fill-due" in out
    body = issue_out.read_text()
    assert "FKS1 fill is due" in body and "RUNBOOK_fks1_fill.md" in body


def test_plan_json_matches_text_counts(tmp_path):
    snap = tmp_path / "PDG000000067.801.tsv"
    records = [_rec("A", "pending:sra-fks1-recaller", run="SRR1"),
               _rec("B", "sra-fks1-recaller:called", call="S639F")]
    _write_snapshot(snap, records)
    script = ROOT / "scripts" / "recall_fks1.py"
    r = subprocess.run(
        [sys.executable, str(script), "plan", str(snap), "--json", "--no-window"],
        capture_output=True, text=True, cwd=str(ROOT))
    assert r.returncode == 0, r.stderr
    data = json.loads(r.stdout)
    assert data["total_isolates"] == 2
    assert data["pending"] == 1 and data["fill_would_recall"] == 1


def _write_snapshot(path, records):
    """Write a snapshot TSV the earlywarning reader round-trips (all SNAPSHOT_COLUMNS present)."""
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(ew.SNAPSHOT_COLUMNS), delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in records:
            w.writerow({c: r.get(c, "") for c in ew.SNAPSHOT_COLUMNS})
