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


# --- country column + trend readers (deepening) ------------------------------

def test_call_table_carries_country_forward():
    # `country` is the newest schema column: a fill that has it writes it, so a future
    # geo breakdown is reproducible from the committed table.
    assert "country" in ew.FKS1_CALLS_COLUMNS
    rec = _rec("A", "sra-fks1-recaller:called", call="S639F")
    rec["country"] = "India"
    text = ew.serialize_fks1_calls([rec])
    assert "country" in text.splitlines()[0]   # header
    assert "India" in text


def test_reader_tolerates_pre_country_tables(tmp_path):
    # Releases 709/710 were written before `country` existed: a 5-column table must still
    # read, overlay and count cleanly (missing country simply absent).
    path = tmp_path / "PDG000000067.700.tsv"
    path.write_text("isolate_key\tfks1_call\tfks1_resistance_source\trun_acc\t"
                    "target_creation_date\n"
                    "A\tS639F\tsra-fks1-recaller:called\tSRR1\t2025-02-01\n"
                    "B\t\tsra-fks1-recaller:called\tSRR2\t2025-02-02\n")
    table = ew.read_fks1_call_table(str(path))
    assert set(table) == {"A", "B"} and "country" not in table["A"]
    live = [_rec("A", "pending:sra-fks1-recaller")]
    assert ew.overlay_fks1_calls(live, str(path)) == 1 and live[0]["fks1_call"] == "S639F"


def test_read_calls_index_returns_full_series_ascending(tmp_path):
    idx = tmp_path / "INDEX.tsv"
    for rel, as_of in (("710", "2026-04-14"), ("709", "2026-04-11"), ("12", "2024-01-01")):
        ew.update_fks1_calls_index(str(idx), {"release_tag": f"PDG000000067.{rel}",
                                              "as_of": as_of, "window_days": 180,
                                              "n_called": 1, "sha256": "x",
                                              "written_at_utc": "2026-10-07T00:00:00Z"})
    rows = ew.read_fks1_calls_index(str(tmp_path))
    assert [r["release_tag"] for r in rows] == [       # ascending by release number, not lexical
        "PDG000000067.12", "PDG000000067.709", "PDG000000067.710"]


def test_read_calls_index_empty_when_missing(tmp_path):
    assert ew.read_fks1_calls_index(str(tmp_path)) == []


def test_panel_positive_separates_known_resistance_from_uncertain_variants(tmp_path):
    # count_fks1_panel_positive counts only KNOWN-resistance panel substitutions, carrying
    # non-panel hot-spot variants as VUS so the page never headlines a variant of uncertain
    # significance as "resistant".
    path = tmp_path / "PDG000000067.800.tsv"
    path.write_text("isolate_key\tfks1_call\tfks1_resistance_source\trun_acc\t"
                    "target_creation_date\n"
                    "A\tS639F\tsra-fks1-recaller:called\tSRR1\t2026-02-01\n"   # known resistance
                    "B\tD642Y\tsra-fks1-recaller:called\tSRR2\t2026-02-02\n"   # known resistance
                    "C\tR1354H\tsra-fks1-recaller:called\tSRR3\t2026-02-03\n"  # non-panel (VUS)
                    "D\t\tsra-fks1-recaller:wild-type\tSRR4\t2026-02-04\n")    # wild-type
    assert ew.count_fks1_panel_positive(str(path)) == (2, 1)  # (S639F+D642Y known, R1354H VUS)


def test_count_panel_positive_missing_table_is_zero(tmp_path):
    assert ew.count_fks1_panel_positive(str(tmp_path / "nope.tsv")) == (0, 0)


def test_vus_census_breaks_down_non_panel_variants_by_token(tmp_path):
    # The per-variant detail behind count_fks1_panel_positive's n_vus: distinct non-panel
    # tokens with carrier counts, descending by count then token. Panel tokens (incl. the
    # low-PPV D642Y/M690I letters) are excluded.
    path = tmp_path / "PDG000000067.800.tsv"
    path.write_text("isolate_key\tfks1_call\tfks1_resistance_source\trun_acc\t"
                    "target_creation_date\n"
                    "A\tS639F\tsra-fks1-recaller:called\tSRR1\t2026-02-01\n"   # panel (excluded)
                    "B\tD642Y\tsra-fks1-recaller:called\tSRR2\t2026-02-02\n"   # low-PPV panel (excluded)
                    "C\tR1354H\tsra-fks1-recaller:called\tSRR3\t2026-02-03\n"  # VUS
                    "D\tR1354H\tsra-fks1-recaller:called\tSRR4\t2026-02-04\n"  # VUS (same token)
                    "E\tD642H\tsra-fks1-recaller:called\tSRR5\t2026-02-05\n"   # VUS
                    "F\t\tsra-fks1-recaller:wild-type\tSRR6\t2026-02-06\n")    # wild-type
    assert ew.fks1_vus_census(str(path)) == [
        {"token": "R1354H", "carriers": 2},  # highest count first
        {"token": "D642H", "carriers": 1},   # ties broken by token
    ]


def test_vus_census_missing_table_is_empty(tmp_path):
    assert ew.fks1_vus_census(str(tmp_path / "nope.tsv")) == []


def test_panel_census_breaks_down_known_resistance_markers_by_token(tmp_path):
    # The per-marker detail behind count_fks1_panel_positive's n_resistant: distinct PANEL
    # tokens with carrier counts, descending by count then token, each tagged with its PPV
    # tier. Non-panel (VUS) tokens are excluded. The low-PPV D642Y/M690I letters stay on the
    # panel (counted in the headline) but carry low_ppv=True.
    path = tmp_path / "PDG000000067.800.tsv"
    path.write_text("isolate_key\tfks1_call\tfks1_resistance_source\trun_acc\t"
                    "target_creation_date\n"
                    "A\tS639F\tsra-fks1-recaller:called\tSRR1\t2026-02-01\n"   # core panel
                    "B\tS639F\tsra-fks1-recaller:called\tSRR2\t2026-02-02\n"   # core panel (same)
                    "C\tD642Y\tsra-fks1-recaller:called\tSRR3\t2026-02-03\n"   # low-PPV panel
                    "D\tR1354H\tsra-fks1-recaller:called\tSRR4\t2026-02-04\n"  # VUS (excluded)
                    "E\t\tsra-fks1-recaller:wild-type\tSRR5\t2026-02-05\n")    # wild-type
    assert ew.fks1_panel_census(str(path)) == [
        {"token": "S639F", "carriers": 2, "low_ppv": False},  # highest count first, core tier
        {"token": "D642Y", "carriers": 1, "low_ppv": True},   # low-PPV panel member tagged
    ]


def test_panel_census_missing_table_is_empty(tmp_path):
    assert ew.fks1_panel_census(str(tmp_path / "nope.tsv")) == []


def test_panel_census_carriers_reconcile_with_headline(tmp_path):
    # Invariant: when no isolate carries more than one panel marker, the per-token carrier
    # counts sum to the headline n_resistant (count_fks1_panel_positive[0]). Same source of
    # truth (iter_resolved_calls), so page and headline can never drift.
    path = tmp_path / "PDG000000067.800.tsv"
    path.write_text("isolate_key\tfks1_call\tfks1_resistance_source\trun_acc\t"
                    "target_creation_date\n"
                    "A\tS639F\tsra-fks1-recaller:called\tSRR1\t2026-02-01\n"
                    "B\tD642Y\tsra-fks1-recaller:called\tSRR2\t2026-02-02\n"
                    "C\tW691L\tsra-fks1-recaller:called\tSRR3\t2026-02-03\n"
                    "D\tR1354H\tsra-fks1-recaller:called\tSRR4\t2026-02-04\n")  # VUS, not counted
    n_resistant, _n_vus = ew.count_fks1_panel_positive(str(path))
    census = ew.fks1_panel_census(str(path))
    assert sum(e["carriers"] for e in census) == n_resistant == 3


def test_iter_resolved_calls_splits_panel_from_vus():
    # The shared primitive all three FKS1 counters + the geo map consume: it splits each
    # record's parsed call tokens into (panel-positive, VUS) on fks1_caller.is_panel_token.
    records = [
        {"isolate_key": "A", "fks1_call": "S639F"},          # panel only
        {"isolate_key": "B", "fks1_call": "R1354H"},         # VUS only
        {"isolate_key": "C", "fks1_call": "S639F,R1354H"},   # mixed row
        {"isolate_key": "D", "fks1_call": ""},               # unresolved
    ]
    out = {r["isolate_key"]: (panel, vus) for r, panel, vus in ew.iter_resolved_calls(records)}
    assert out["A"] == (["S639F"], [])
    assert out["B"] == ([], ["R1354H"])
    assert out["C"] == (["S639F"], ["R1354H"])
    assert out["D"] == ([], [])          # yielded, but both lists empty -> not a carrier


def test_release_sort_key_orders_by_release_number_not_lexically():
    tags = ["PDG000000067.12", "PDG000000067.709", "PDG000000067.710", "snapshot-x"]
    assert sorted(tags, key=ew._release_sort_key) == [
        "PDG000000067.12", "PDG000000067.709", "PDG000000067.710",  # numeric, not "12">"709"
        "snapshot-x",                                               # non-release tag sorts last
    ]


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


# --- geography backfill + historical breakdown (issue #170) ------------------
# The `country` column was added to the call-table schema after releases 709/710 were written,
# so their tables carry it empty; `backfill_fks1_calls_country` fills it from a snapshot's
# geography and `fks1_calls_geo_breakdown` reads it back, panel-gated, for the page's
# per-release historical geography.

def _geo_rec(key, call, source="sra-fks1-recaller:called", country="", geo_loc_name=""):
    r = _rec(key, source, call=call)
    r["country"] = country
    r["geo_loc_name"] = geo_loc_name
    return r


def test_backfill_fills_country_from_map_and_leaves_unknown_blank(tmp_path):
    path = tmp_path / "PDG000000067.709.tsv"
    ew.write_fks1_call_table(str(path), [
        _rec("A", "sra-fks1-recaller:called", call="S639F"),
        _rec("B", "sra-fks1-recaller:called", call="F635I"),
    ])
    sha, n, filled, missing = ew.backfill_fks1_calls_country(
        str(path), {"A": "USA", "B": ""})   # B maps to blank -> unlocated, never invented
    assert n == 2 and filled == 1 and missing == ["B"]
    table = ew.read_fks1_call_table(str(path))
    assert table["A"]["country"] == "USA"
    assert table["B"]["country"] == ""
    # row set unchanged, and byte-stable: re-reading + re-serializing is identical
    sha2, _ = ew.write_fks1_call_table(str(path), list(table.values()))
    assert sha2 == sha


def test_backfill_preserves_row_count_and_other_columns(tmp_path):
    path = tmp_path / "PDG000000067.709.tsv"
    ew.write_fks1_call_table(str(path), [_rec("A", "sra-fks1-recaller:called", call="S639F")])
    ew.backfill_fks1_calls_country(str(path), {"A": "China"})
    table = ew.read_fks1_call_table(str(path))
    assert set(table) == {"A"}                       # no rows added or dropped
    assert table["A"]["fks1_call"] == "S639F"         # existing columns intact
    assert table["A"]["run_acc"] == "SRR1"


def test_geo_breakdown_counts_panel_carriers_by_country(tmp_path):
    path = tmp_path / "PDG000000067.710.tsv"
    ew.write_fks1_call_table(str(path), [
        _geo_rec("A", "S639F", country="USA"),      # panel -> counted
        _geo_rec("B", "S639F", country="USA"),      # panel -> counted
        _geo_rec("C", "W691L", country="China"),    # panel -> counted
        _geo_rec("D", "F635I", country="USA"),      # non-panel VUS -> NOT counted
        _geo_rec("E", "S639F", country=""),         # located nowhere -> omitted
    ])
    assert ew.fks1_calls_geo_breakdown(str(path)) == {"USA": 2, "China": 1}


def test_geo_breakdown_empty_when_country_absent(tmp_path):
    # a pre-backfill table (country cells blank) yields {} -- honest degradation, and a
    # missing table reads as {} rather than raising
    path = tmp_path / "PDG000000067.700.tsv"
    ew.write_fks1_call_table(str(path), [_rec("A", "sra-fks1-recaller:called", call="S639F")])
    assert ew.fks1_calls_geo_breakdown(str(path)) == {}
    assert ew.fks1_calls_geo_breakdown(str(tmp_path / "missing.tsv")) == {}


def test_backfill_fills_geo_loc_name_from_dict_map(tmp_path):
    # a dict-valued map restores the finest place (geo_loc_name) alongside country; a blank
    # country leaves geo_loc_name blank too (never a region without its country)
    path = tmp_path / "PDG000000067.710.tsv"
    ew.write_fks1_call_table(str(path), [
        _rec("A", "sra-fks1-recaller:called", call="S639F"),
        _rec("B", "sra-fks1-recaller:called", call="W691L"),
    ])
    sha, n, filled, missing = ew.backfill_fks1_calls_country(str(path), {
        "A": {"country": "USA", "geo_loc_name": "USA: California"},
        "B": {"country": "", "geo_loc_name": "Somewhere"},   # unlocated country -> both blank
    })
    assert n == 2 and filled == 1 and missing == ["B"]
    table = ew.read_fks1_call_table(str(path))
    assert table["A"]["country"] == "USA"
    assert table["A"]["geo_loc_name"] == "USA: California"
    assert table["B"]["country"] == "" and table["B"]["geo_loc_name"] == ""
    # column-add-only + byte-stable
    sha2, _ = ew.write_fks1_call_table(str(path), list(table.values()))
    assert sha2 == sha


def test_backfill_string_map_stays_country_only(tmp_path):
    # legacy string-valued map still works and leaves geo_loc_name blank
    path = tmp_path / "PDG000000067.709.tsv"
    ew.write_fks1_call_table(str(path), [_rec("A", "sra-fks1-recaller:called", call="S639F")])
    ew.backfill_fks1_calls_country(str(path), {"A": "USA"})
    table = ew.read_fks1_call_table(str(path))
    assert table["A"]["country"] == "USA" and table["A"]["geo_loc_name"] == ""


def test_geo_breakdown_finest_prefers_region_else_country(tmp_path):
    path = tmp_path / "PDG000000067.710.tsv"
    ew.write_fks1_call_table(str(path), [
        _geo_rec("A", "S639F", country="USA", geo_loc_name="USA: California"),  # region kept
        _geo_rec("B", "S639F", country="USA", geo_loc_name="USA: California"),  # region kept
        _geo_rec("C", "W691L", country="USA", geo_loc_name="USA: Texas"),       # distinct region
        _geo_rec("D", "S639F", country="China", geo_loc_name=""),               # falls back to country
        _geo_rec("E", "F635I", country="USA", geo_loc_name="USA: California"),  # non-panel -> omitted
    ])
    assert ew.fks1_calls_geo_breakdown(str(path), finest=True) == {
        "USA: California": 2, "USA: Texas": 1, "China": 1}
    # the default (country-level) view collapses the regions
    assert ew.fks1_calls_geo_breakdown(str(path)) == {"USA": 3, "China": 1}
