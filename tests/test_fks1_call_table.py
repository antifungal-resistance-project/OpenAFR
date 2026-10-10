"""Tests for the FKS1 call-table durability layer + windowed fill scope (issue #171).

The echinocandin re-caller runs billably on GCP and its calls are NOT regenerable without
re-paying; the snapshot blob they live in is gitignored. So the calls get their own small,
committed home (openafr.earlywarning.{write,read}_fks1_call_table + overlay_fks1_calls), and
the fill is scoped to the trailing-180-day window the weather page shows
(scripts/recall_fks1.py _in_window / _fill_order). These lock both halves and the end-to-end
transition: a live pull that overlays the committed calls LEAVES the day-0 watching state.

No network, no external tools -- pure functions + the shared #22-style FKS1 fixture.
"""
import importlib.util
import pathlib
import sys

from openafr import alert
from openafr import earlywarning as ew

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
from compose_alert import _demo_fks1_records  # noqa: E402  (the FKS1 fixture)

_SCRIPT = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "recall_fks1.py"


def _load_recall_fks1():
    spec = importlib.util.spec_from_file_location("recall_fks1", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _rec(key, source, call="", run="SRR1", created="2025-02-01"):
    return {"isolate_key": key, "target_acc": key + ".1", "run_acc": run,
            "target_creation_date": created, "fks1_call": call,
            "fks1_resistance_source": source}


# --- call table: write / read / round-trip -----------------------------------

def test_call_table_holds_only_resolved_rows(tmp_path):
    records = [
        _rec("A", "sra-fks1-recaller:called", call="S639F"),
        _rec("B", "pending:sra-fks1-recaller"),          # not spoken for -> excluded
        _rec("C", "pending:no-reads"),                    # excluded
        _rec("D", "sra-fks1-recaller:failed(fetch)"),     # resolved (re-caller spoke) -> kept
    ]
    path = tmp_path / "PDG000000067.700.tsv"
    sha, n = ew.write_fks1_call_table(str(path), records)
    assert n == 2
    table = ew.read_fks1_call_table(str(path))
    assert set(table) == {"A", "D"}
    assert table["A"]["fks1_call"] == "S639F"
    # byte-stable: same resolved rows -> same sha
    sha2, _ = ew.write_fks1_call_table(str(path), records)
    assert sha2 == sha


def test_overlay_replays_calls_and_leaves_others_pending():
    live = [_rec("A", "pending:sra-fks1-recaller"),
            _rec("B", "pending:sra-fks1-recaller"),
            _rec("C", "pending:no-reads")]
    calls = {"A": {"isolate_key": "A", "fks1_call": "S639F",
                   "fks1_resistance_source": "sra-fks1-recaller:called"}}
    n = ew.overlay_fks1_calls(live, calls)
    assert n == 1
    by = {r["isolate_key"]: r for r in live}
    assert by["A"]["fks1_call"] == "S639F"
    assert by["A"]["fks1_resistance_source"] == "sra-fks1-recaller:called"
    # unmatched rows keep their honest pending state -- never silently wild-type
    assert by["B"]["fks1_resistance_source"] == "pending:sra-fks1-recaller"
    assert by["C"]["fks1_resistance_source"] == "pending:no-reads"


def test_overlay_accepts_a_path(tmp_path):
    path = tmp_path / "PDG000000067.700.tsv"
    ew.write_fks1_call_table(str(path), [_rec("A", "sra-fks1-recaller:called", call="S639F")])
    live = [_rec("A", "pending:sra-fks1-recaller")]
    assert ew.overlay_fks1_calls(live, str(path)) == 1
    assert live[0]["fks1_call"] == "S639F"


def test_latest_call_table_picks_highest_release(tmp_path):
    for rel in (670, 700, 686):
        ew.write_fks1_call_table(str(tmp_path / f"PDG000000067.{rel}.tsv"), [])
    assert ew.latest_fks1_call_table(str(tmp_path)).endswith("PDG000000067.700.tsv")


def test_latest_call_table_none_when_empty(tmp_path):
    assert ew.latest_fks1_call_table(str(tmp_path)) is None


def test_calls_index_upserts_per_release(tmp_path):
    idx = tmp_path / "INDEX.tsv"
    ew.update_fks1_calls_index(str(idx), {"release_tag": "PDG000000067.700", "as_of": "2025-01-01",
                                          "window_days": 180, "n_called": 3, "sha256": "aa",
                                          "written_at_utc": "2026-10-07T00:00:00Z"})
    ew.update_fks1_calls_index(str(idx), {"release_tag": "PDG000000067.700", "as_of": "2025-02-01",
                                          "window_days": 180, "n_called": 5, "sha256": "bb",
                                          "written_at_utc": "2026-10-07T01:00:00Z"})
    rows = list(__import__("csv").DictReader(idx.open(), delimiter="\t"))
    assert len(rows) == 1 and rows[0]["n_called"] == "5"   # re-fill replaced the row


def test_latest_index_entry_picks_highest_release(tmp_path):
    idx = tmp_path / "INDEX.tsv"
    for rel, as_of in (("300", "2025-01-01"), ("709", "2026-04-11"), ("12", "2024-01-01")):
        ew.update_fks1_calls_index(str(idx), {"release_tag": f"PDG000000067.{rel}",
                                              "as_of": as_of, "window_days": 180,
                                              "n_called": 1, "sha256": "x",
                                              "written_at_utc": "2026-10-07T00:00:00Z"})
    entry = ew.latest_fks1_calls_index_entry(str(tmp_path))
    assert entry["release_tag"] == "PDG000000067.709"   # highest release number, not lexical
    assert entry["as_of"] == "2026-04-11"


def test_latest_index_entry_none_when_missing_or_empty(tmp_path):
    assert ew.latest_fks1_calls_index_entry(str(tmp_path)) is None   # no INDEX.tsv
    (tmp_path / "INDEX.tsv").write_text(
        "\t".join(ew.FKS1_CALLS_INDEX_COLUMNS) + "\n")                # header only
    assert ew.latest_fks1_calls_index_entry(str(tmp_path)) is None


# --- windowed fill scope (recall_fks1.py) ------------------------------------

def test_in_window_restricts_to_recent_partition():
    mod = _load_recall_fks1()
    records = [_rec("OLD", "pending:sra-fks1-recaller", created="2020-01-01"),
               _rec("NEW", "pending:sra-fks1-recaller", created="2025-02-01")]
    # as_of 2025-01-01, window 180d -> recent = created in (2025-01-01, 2025-06-30]
    keys = mod._in_window(records, "2025-01-01", 180)
    assert keys == {"NEW"}


def test_fill_order_scoped_to_window_keys():
    mod = _load_recall_fks1()
    records = [_rec("OLD", "pending:sra-fks1-recaller", created="2020-01-01"),
               _rec("NEW", "pending:sra-fks1-recaller", created="2025-02-01"),
               _rec("NOREADS", "pending:no-reads", created="2025-02-01", run="")]
    keys = mod._in_window(records, "2025-01-01", 180)
    attempted = mod._fill_order(records, window_keys=keys)
    assert [r["isolate_key"] for r in attempted] == ["NEW"]
    # no window filter -> both pending-with-run rows
    assert {r["isolate_key"] for r in mod._fill_order(records)} == {"OLD", "NEW"}


# --- end-to-end: overlay lifts a live pull out of day-0 ----------------------

def test_overlay_lights_up_the_page(tmp_path):
    """A raw live pull carries no FKS1 calls (day-0); overlaying the committed table makes
    compose_fks1_alerts fire -- the loop issue #171 closes."""
    called = _demo_fks1_records()
    as_of, window = "2024-12-31", 180

    # the committed call table, built from a prior (billable) fill
    table = tmp_path / "PDG000000067.700.tsv"
    ew.write_fks1_call_table(str(table), called)

    # a fresh live pull: same isolates, every FKS1 call reset to the honest pending state
    live = []
    for r in called:
        p = dict(r)
        p["fks1_call"] = ""
        p["fks1_resistance_source"] = "pending:sra-fks1-recaller"
        live.append(p)

    before = alert.compose_fks1_alerts(live, as_of, window_days=window,
                                       min_count=3, min_delta=0.05)
    assert before["alerts"] == []            # day-0: nothing to show

    ew.overlay_fks1_calls(live, str(table))
    after = alert.compose_fks1_alerts(live, as_of, window_days=window,
                                      min_count=3, min_delta=0.05)
    assert after["alerts"], "overlaying the committed FKS1 calls should light up the page"
