"""Tests for the echinocandin (FKS1) event frequency (openafr/backtest.py, v2/T3).

The FKS1 analog of the azole panel_prevalence, with the same honesty contract:
  * only isolates whose FKS1 panel could be RESOLVED enter the denominator,
  * per-window statuses (HS1=..., HS2=..., HS3=...) are classified over the panel-bearing
    windows -- all three now carry markers, so a partial/refused/missing panel window is
    excluded, never counted negative,
  * panel positivity is decided from the tokens via fks1_caller.is_panel_token, so an
    isolate called only for a clade/non-panel SNP is resolved-but-negative,
  * an un-filled snapshot reports 0 resolved / event_frequency None -- NOT MEASURED YET,
    never a false 0%.
No network: synthetic snapshot records only.
"""
from openafr import backtest as bt


def _rec(key, fks1_call, source, run="SRR9", geo=None):
    rec = {"isolate_key": key, "target_acc": key + ".1", "run_acc": run,
           "fks1_call": fks1_call, "fks1_resistance_source": source}
    if geo is not None:
        rec["geo_loc_name"] = geo
    return rec


def test_resolution_buckets_over_panel_window():
    b = bt._fks1_resolution_bucket
    assert b("sra-fks1-recaller:HS1=called,HS2=wild-type,HS3=wild-type") == "resolved"
    assert b("sra-fks1-recaller:HS1=wild-type,HS2=wild-type,HS3=wild-type") == "resolved"
    # all three windows now carry panel markers, so an uncalled/missing panel window blocks
    # resolution -- HS2 or HS3 missing is no longer "resolved" the way it was pre-extension
    assert b("sra-fks1-recaller:HS1=called,HS2=missing,HS3=wild-type") == "partial"
    assert b("sra-fks1-recaller:HS1=called,HS2=wild-type,HS3=missing") == "partial"
    # a partial/missing/refused panel window is not resolved
    assert b("sra-fks1-recaller:HS1=partial(uncalled:639),HS2=wild-type,HS3=wild-type") == "partial"
    assert b("sra-fks1-recaller:HS1=missing,HS2=wild-type,HS3=wild-type") == "partial"
    assert b("sra-fks1-recaller:HS1=refused(indel),HS2=wild-type,HS3=wild-type") == "failed"
    assert b("sra-fks1-recaller:HS1=wild-type,HS2=wild-type,HS3=refused(indel)") == "failed"
    # whole-isolate failure (no per-window part) and pending states
    assert b("sra-fks1-recaller:failed(no reads)") == "failed"
    assert b("pending:sra-fks1-recaller") == "pending"
    assert b("pending:no-reads") == "pending"
    assert b("") == "pending"


def test_prevalence_counts_only_resolved_and_reads_panel_from_tokens():
    records = [
        _rec("A", "S639F", "sra-fks1-recaller:HS1=called,HS2=wild-type,HS3=wild-type"),  # panel-positive
        _rec("B", "S639P", "sra-fks1-recaller:HS1=called,HS2=wild-type,HS3=wild-type"),  # panel-positive
        _rec("C", "", "sra-fks1-recaller:HS1=wild-type,HS2=wild-type,HS3=wild-type"),    # resolved, negative
        _rec("D", "D642E", "sra-fks1-recaller:HS1=called,HS2=wild-type,HS3=wild-type"),  # resolved, non-panel (D642E != D642Y)
        _rec("E", "", "sra-fks1-recaller:HS1=partial(uncalled:639),HS2=wild-type,HS3=wild-type"),  # partial
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


_RES = "sra-fks1-recaller:HS1=called,HS2=wild-type,HS3=wild-type"
_WT = "sra-fks1-recaller:HS1=wild-type,HS2=wild-type,HS3=wild-type"


def test_prevalence_by_place_uses_resolved_denominator_not_carrier_count():
    # USA: 2 resolved, 1 panel-positive (50%). UK: 5 resolved, 1 positive (20%). A region's
    # denominator is its RESOLVED isolates, never its carrier count: the partial/pending USA
    # rows below are excluded from USA's denominator, so volume can't inflate the rate.
    records = [
        _rec("a", "S639F", _RES, geo="USA: California"),
        _rec("b", "", _WT, geo="USA: California"),
        _rec("c", "", "sra-fks1-recaller:HS1=partial(uncalled:639),HS2=wild-type,HS3=wild-type",
             geo="USA: California"),                       # excluded, not a negative
    ] + [_rec(f"u{i}", "S639P" if i == 0 else "", _RES, geo="United Kingdom")
         for i in range(5)]
    out = bt.fks1_panel_prevalence_by_place(records, finest=True, min_resolved=2)
    places = {p["place"]: p for p in out["places"]}
    assert places["USA: California"]["n_resolved"] == 2
    assert places["USA: California"]["n_panel_positive"] == 1
    assert places["USA: California"]["event_frequency"] == 0.5
    assert places["United Kingdom"]["n_resolved"] == 5
    assert places["United Kingdom"]["event_frequency"] == 0.2
    # Each place carries a per_mutation tally of which marker its carriers hold.
    assert places["USA: California"]["per_mutation"] == {"S639F": 1}
    assert places["United Kingdom"]["per_mutation"] == {"S639P": 1}
    # Ranked by prevalence descending: the higher-rate place leads despite its smaller N.
    assert out["places"][0]["place"] == "USA: California"
    assert all(p["ci95"] is not None for p in out["places"])


def test_prevalence_by_place_per_mutation_counts_each_marker_and_sorts_by_position():
    # A multi-marker isolate is counted under EACH token (so per_mutation is a token-carrier
    # tally, not a partition of n_panel_positive), and tokens sort by position number.
    records = [
        _rec("a", "M690I/S639F", _RES, geo="USA: Mountain"),   # two panel markers, one isolate
        _rec("b", "S639F", _RES, geo="USA: Mountain"),
        _rec("c", "", _WT, geo="USA: Mountain"),
    ]
    out = bt.fks1_panel_prevalence_by_place(records, finest=True, min_resolved=2)
    place = out["places"][0]
    assert place["n_panel_positive"] == 2                       # two isolates, not three tokens
    assert place["per_mutation"] == {"S639F": 2, "M690I": 1}
    assert list(place["per_mutation"]) == ["S639F", "M690I"]    # sorted by position (639 < 690)


def test_prevalence_by_place_pools_small_places_and_drops_unplaced():
    records = [
        _rec("a", "S639F", _RES, geo="USA"),
        _rec("b", "", _WT, geo="USA"),
        _rec("small", "S639F", _RES, geo="Tinyland"),      # 1 resolved < min -> pooled
        _rec("nowhere", "S639F", _RES),                    # no geo -> unplaced, dropped
    ]
    out = bt.fks1_panel_prevalence_by_place(records, finest=True, min_resolved=2)
    assert [p["place"] for p in out["places"]] == ["USA"]
    assert out["pooled_small"] == {"n_places": 1, "n_resolved": 1, "n_panel_positive": 1}
    assert out["n_unplaced"] == 1
    # SSoT invariant: per-place + pooled positives == the national headline positives.
    headline = bt.fks1_panel_prevalence(records)
    summed = (sum(p["n_panel_positive"] for p in out["places"])
              + out["pooled_small"]["n_panel_positive"]
              + sum(1 for r in records                       # unplaced positives
                    if r.get("geo_loc_name") is None
                    and "S639" in (r.get("fks1_call") or "")))
    assert summed == headline["n_panel_positive"]
