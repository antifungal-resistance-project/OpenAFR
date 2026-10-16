"""Tests for the FKS1 (echinocandin) real-data backtest additions to openafr/backtest.py.

These lock the guarantees the FKS1 backtest leans on
(work/PREREGISTRATION_fks1_backtest.md):
  * walk_forward reads the gene selected by call_field -- on the SAME records it runs the
    real FKS1 detector (fks1_call) but stays the honest null on the empty azole column,
  * load_fks1_call_table turns a committed call table into detector records that
    fks1_panel_prevalence scores with the right numerator/denominator,
  * the two CLI subcommands run end-to-end on the committed release-710 table and exit 0.
No network: synthetic records plus the committed, offline call table.
"""
import pathlib
import subprocess
import sys

from openafr import backtest as bt

ROOT = pathlib.Path(__file__).resolve().parent.parent
REAL_TABLE = ROOT / "data" / "earlywarning" / "fks1_calls" / "PDG000000067.710.tsv"

# A resolved FKS1 source: every panel window (HS1/HS2/HS3) called or wild-type.
_RESOLVED = "sra-fks1-recaller:HS1=called,HS2=wild-type,HS3=wild-type"
_RESOLVED_WT = "sra-fks1-recaller:HS1=wild-type,HS2=wild-type,HS3=wild-type"


def _fks1_rec(key, created, call="", source=_RESOLVED_WT):
    return {
        "isolate_key": key,
        "target_creation_date": created,
        "fks1_call": call,
        "fks1_resistance_source": source,
    }


# ---- walk_forward is gene-selectable -----------------------------------------

def test_walk_forward_reads_fks1_calls_not_the_empty_azole_column():
    # A cohort where S639F keeps appearing across months, called on the FKS1 column only.
    recs = []
    for i in range(8):
        recs.append(_fks1_rec(f"WT{i:02d}", f"2026-0{1 + i % 6}-10"))
    for i in range(4):
        recs.append(_fks1_rec(f"R{i:02d}", f"2026-0{3 + i}-10",
                               call="S639F", source=_RESOLVED))
    steps = bt.month_steps(*bt.creation_date_range(recs))

    fks1 = bt.walk_forward(recs, steps, target="S639F", call_field="fks1_call")
    assert fks1["any_callable"] is True            # real FKS1 calls -> detector runs
    assert fks1["first_flag"] is not None          # S639F is flagged

    # Same records, default (azole) column is empty -> the honest null, no flag.
    azole = bt.walk_forward(recs, steps, target="S639F")   # call_field defaults to erg11_call
    assert azole["any_callable"] is False
    assert azole["first_flag"] is None


# ---- load_fks1_call_table + fks1_panel_prevalence ----------------------------

def test_load_table_round_trips_and_prevalence_scores(tmp_path):
    table = tmp_path / "PDG000000067.999.tsv"
    rows = [
        "isolate_key\tfks1_call\tfks1_resistance_source\trun_acc\ttarget_creation_date",
        f"A\tS639F\t{_RESOLVED}\tSRR1\t2026-04-10",   # panel-positive
        f"B\tD642Y\t{_RESOLVED}\tSRR2\t2026-05-10",   # panel-positive
        f"C\t\t{_RESOLVED_WT}\tSRR3\t2026-05-11",     # resolved, wild-type
        f"D\tR1354H\t{_RESOLVED}\tSRR4\t2026-06-10",  # resolved, NON-panel only
    ]
    table.write_text("\n".join(rows) + "\n")

    recs = bt.load_fks1_call_table(str(table))
    assert len(recs) == 4
    assert {r["isolate_key"] for r in recs} == {"A", "B", "C", "D"}

    p = bt.fks1_panel_prevalence(recs)
    assert p["n_resolved"] == 4
    assert p["n_panel_positive"] == 2                 # S639F, D642Y
    assert p["n_called_nonpanel"] == 1                # R1354H is not a known-resistance panel token
    assert abs(p["event_frequency"] - 2 / 4) < 1e-9
    assert set(p["per_mutation"]) == {"S639F", "D642Y"}


# ---- CLI smoke on the committed table ----------------------------------------

def _run(*args):
    return subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "backtest_earlywarning.py"), *args],
        capture_output=True, text=True, cwd=str(ROOT))


def test_cli_fks1_prevalence_runs_on_committed_table():
    r = _run("fks1-prevalence", str(REAL_TABLE))
    assert r.returncode == 0, r.stderr
    assert "echinocandin event frequency" in r.stdout.lower()


def test_cli_fks1_replay_runs_on_committed_table():
    r = _run("fks1-replay", str(REAL_TABLE))
    assert r.returncode == 0, r.stderr
    assert "walk-forward" in r.stdout.lower()
    assert "retired for fks1" in r.stdout.lower()
