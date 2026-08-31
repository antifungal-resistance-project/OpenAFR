"""Tests for the echinocandin (FKS1) event frequency (openafr/backtest.py, v2/T3).

The FKS1 analog of the azole panel_prevalence, with the same honesty contract:
  * only isolates whose FKS1 panel could be RESOLVED enter the denominator,
  * per-window statuses (HS1=..., HS2=...) are classified over the panel-bearing
    window(s) only -- a partial/refused/missing HS1 is excluded, never counted negative,
  * panel positivity is decided from the tokens via fks1_caller.is_panel_token, so an
    isolate called only for a clade/non-panel SNP is resolved-but-negative,
  * an un-filled snapshot reports 0 resolved / event_frequency None -- NOT MEASURED YET,
    never a false 0%.
No network: synthetic snapshot records only.
"""
from openafr import backtest as bt


def _rec(key, fks1_call, source, run="SRR9"):
    return {"isolate_key": key, "target_acc": key + ".1", "run_acc": run,
            "fks1_call": fks1_call, "fks1_resistance_source": source}


def test_resolution_buckets_over_panel_window():
    b = bt._fks1_resolution_bucket
    assert b("sra-fks1-recaller:HS1=called,HS2=wild-type") == "resolved"
    assert b("sra-fks1-recaller:HS1=wild-type,HS2=wild-type") == "resolved"
    # HS2 uncalled does NOT block resolution: the panel lives in HS1 only
    assert b("sra-fks1-recaller:HS1=called,HS2=missing") == "resolved"
    # a partial/missing/refused HS1 (the panel window) is not resolved
    assert b("sra-fks1-recaller:HS1=partial(uncalled:639),HS2=wild-type") == "partial"
    assert b("sra-fks1-recaller:HS1=missing,HS2=wild-type") == "partial"
    assert b("sra-fks1-recaller:HS1=refused(indel),HS2=wild-type") == "failed"
    # whole-isolate failure (no per-window part) and pending states
    assert b("sra-fks1-recaller:failed(no reads)") == "failed"
    assert b("pending:sra-fks1-recaller") == "pending"
    assert b("pending:no-reads") == "pending"
    assert b("") == "pending"


def test_prevalence_counts_only_resolved_and_reads_panel_from_tokens():
    records = [
        _rec("A", "S639F", "sra-fks1-recaller:HS1=called,HS2=wild-type"),   # panel-positive
        _rec("B", "S639P", "sra-fks1-recaller:HS1=called,HS2=wild-type"),   # panel-positive
        _rec("C", "", "sra-fks1-recaller:HS1=wild-type,HS2=wild-type"),     # resolved, negative
        _rec("D", "D642E", "sra-fks1-recaller:HS1=called,HS2=wild-type"),   # resolved, non-panel
        _rec("E", "", "sra-fks1-recaller:HS1=partial(uncalled:639),HS2=wild-type"),  # partial
        _rec("F", "", "sra-fks1-recaller:failed(align)"),                   # failed
        _rec("G", "", "pending:sra-fks1-recaller"),                         # pending
    ]
    p = bt.fks1_panel_prevalence(records)
    assert p["n_total"] == 7
    assert p["n_resolved"] == 4              # A,B,C,D
    assert p["n_panel_positive"] == 2        # A,B
    assert p["n_called_nonpanel"] == 1       # D (D642E is not on the panel)
    assert p["event_frequency"] == 0.5       # 2/4
    assert p["excluded"] == {"partial": 1, "failed": 1, "pending": 1}
    assert p["per_mutation"] == {"S639F": 1, "S639P": 1}
    assert p["ci95"] is not None


def test_prevalence_is_not_measured_until_filled():
    # A fresh snapshot: every isolate pending -> nothing resolved -> undefined, not 0%.
    records = [_rec(str(i), "", "pending:sra-fks1-recaller") for i in range(50)]
    p = bt.fks1_panel_prevalence(records)
    assert p["n_resolved"] == 0
    assert p["event_frequency"] is None and p["ci95"] is None
    assert p["excluded"]["pending"] == 50
