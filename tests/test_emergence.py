"""Tests for emergence detection (openafr/emergence.py, issue #22).

These lock the guarantees the alert layer (#24/#25) will lean on:
  * a call is parsed into substitutions, never invented or silently dropped,
  * baseline vs recent is split by NCBI visibility date, undated isolates set aside,
  * proportions are taken over CALLED isolates only (an empty call is not
    counted as susceptible),
  * first-appearance and rising both fire, and both respect the min-count floor,
  * the backlog guard catches an archival spike (old collection, fresh visibility),
  * a real-shaped snapshot with no calls yields an honest "cannot detect yet".
No network: everything runs against synthetic snapshot records.
"""
from openafr import emergence as em


def _rec(key, created, call, collected="2024", country="Pakistan"):
    return {
        "isolate_key": key, "target_acc": key + ".1",
        "target_creation_date": created, "collection_date": collected,
        "country": country, "erg11_call": call,
        "run_acc": "SRR" + key[-1], "resistance_source": "recalled",
    }


# ---- parse_call --------------------------------------------------------------

def test_parse_call_splits_and_uppercases():
    assert em.parse_call("Y132F,K143R") == {"Y132F", "K143R"}
    assert em.parse_call("y132f; k143r") == {"Y132F", "K143R"}
    assert em.parse_call("  ") == set()
    assert em.parse_call(None) == set()


def test_parse_call_keeps_noncanonical_tokens():
    # We never silently drop a nonempty call the re-caller emitted.
    toks = em.parse_call("TR34/L98H")
    assert toks == {"TR34", "L98H"}
    assert em.is_substitution("Y132F") and not em.is_substitution("TR34")


# ---- windowing ---------------------------------------------------------------

def test_split_by_creation_buckets_and_sets_aside_undated():
    recs = [
        _rec("PDT1", "2024-01-01", "Y132F"),          # baseline
        _rec("PDT2", "2025-02-01", "Y132F"),          # recent
        _rec("PDT3", "2025-12-01", "Y132F"),          # beyond window -> dropped
        _rec("PDT4", "NULL", "Y132F"),                # undated
    ]
    base, recent, undated = em.split_by_creation(recs, "2024-12-31", window_days=180)
    assert [r["isolate_key"] for r in base] == ["PDT1"]
    assert [r["isolate_key"] for r in recent] == ["PDT2"]
    assert [r["isolate_key"] for r in undated] == ["PDT4"]


# ---- called denominator honesty ---------------------------------------------

def test_uncalled_isolate_is_not_a_susceptible_isolate():
    recs = [_rec("PDT1", "2025-01-01", ""), _rec("PDT2", "2025-01-01", "Y132F")]
    assert [r["isolate_key"] for r in em.called(recs)] == ["PDT2"]
    carriers = em.mutation_carriers(recs)
    assert set(carriers) == {"Y132F"}
    # denominator is the called count (1), never the total (2)
    res = em.detect_emergence(recs, "2024-12-31", min_count=1)
    assert res["coverage"]["recent_called"] == 1
    assert res["signals"][0]["recent_prop"] == 1.0


# ---- first-appearance & rising ----------------------------------------------

def test_first_appearance_fires_and_respects_min_count():
    base = [_rec(f"PDTB{i}", "2024-01-01", "Y132F") for i in range(5)]
    recent = [_rec(f"PDTR{i}", "2025-02-01", "VF125A") for i in range(3)]
    res = em.detect_emergence(base + recent, "2024-12-31", min_count=3)
    sigs = {s["mutation"]: s for s in res["signals"]}
    assert "VF125A" in sigs and sigs["VF125A"]["first_appearance"]
    assert sigs["VF125A"]["baseline_count"] == 0

    # one fewer carrier -> below min_count -> no signal
    res2 = em.detect_emergence(base + recent[:2], "2024-12-31", min_count=3)
    assert res2["signals"] == []


def test_rising_fires_on_proportion_increase():
    # K143R: 1/5 baseline (20%) -> 4/4 recent (100%), delta +80%
    base = ([_rec("PDTB0", "2024-01-01", "K143R")]
            + [_rec(f"PDTB{i}", "2024-01-01", "Y132F") for i in range(1, 5)])
    recent = [_rec(f"PDTR{i}", "2025-02-01", "K143R") for i in range(4)]
    res = em.detect_emergence(base + recent, "2024-12-31", min_count=3, min_delta=0.05)
    k = next(s for s in res["signals"] if s["mutation"] == "K143R")
    assert k["rising"] and not k["first_appearance"]
    assert k["delta"] > 0.5


# ---- backlog guard -----------------------------------------------------------

def test_backlog_guard_flags_archival_spike():
    base = [_rec(f"PDTB{i}", "2024-01-01", "Y132F") for i in range(5)]
    # TR34 spike, but every carrier was COLLECTED in 2019 (archival), just now visible
    recent = [_rec(f"PDTR{i}", "2025-03-01", "TR34", collected="2019")
              for i in range(5)]
    res = em.detect_emergence(base + recent, "2024-12-31", min_count=3)
    tr = next(s for s in res["signals"] if s["mutation"] == "TR34")
    assert tr["possible_backlog"]
    assert "backlog" in tr["backlog_reason"].lower()


def test_genuine_recent_collection_is_not_flagged_backlog():
    base = [_rec(f"PDTB{i}", "2024-01-01", "Y132F") for i in range(5)]
    recent = [_rec(f"PDTR{i}", "2025-02-01", "VF125A", collected="2024")
              for i in range(4)]
    res = em.detect_emergence(base + recent, "2024-12-31", min_count=3)
    vf = next(s for s in res["signals"] if s["mutation"] == "VF125A")
    assert not vf["possible_backlog"]


# ---- honest empty state on real (uncalled) data -----------------------------

def test_no_calls_yields_honest_cannot_detect():
    recs = [_rec(f"PDT{i}", "2025-02-01", "") for i in range(20)]
    res = em.detect_emergence(recs, "2024-12-31")
    assert res["signals"] == []
    assert "note" in res and "cannot detect" in res["note"].lower()
    assert res["coverage"]["recent_called"] == 0


def test_signal_ordering_genuine_before_backlog():
    base = [_rec(f"PDTB{i}", "2024-01-01", "Y132F") for i in range(5)]
    genuine = [_rec(f"PDTG{i}", "2025-02-01", "VF125A", collected="2024")
               for i in range(3)]
    backlog = [_rec(f"PDTK{i}", "2025-03-01", "TR34", collected="2019")
               for i in range(6)]
    res = em.detect_emergence(base + genuine + backlog, "2024-12-31", min_count=3)
    # genuine signal must sort ahead of the (larger) backlog spike
    assert res["signals"][0]["mutation"] == "VF125A"
    assert res["signals"][-1]["possible_backlog"]
