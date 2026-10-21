"""Tests for the weather page's window anchoring (scripts/render_weather.py + openafr/fill_status).

The live FKS1 page shows a frozen, git-tracked fill. A trailing-`today` display window marches
forward and slowly ages the fill's isolates out, so the page silently decays toward "awaiting
re-call" even though the fill is the freshest data there is. The fix anchors the DISPLAY window to
the committed fill's own covered window, so the page stays lit at full coverage regardless of how
much wall-clock time has passed since the fill -- with the billable refill decision left to the
today-relative `fill_status.fill_due` watchdog, which these tests confirm is untouched.

No network, no cloud: the real committed release-710 table stands in for a live pull via a
monkeypatched `_load_records`, and the census fixtures are tiny in-tmp tables.
"""
import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
import render_weather as rw  # noqa: E402

from openafr import earlywarning as ew  # noqa: E402
from openafr import emergence, fill_status  # noqa: E402

CALLS_DIR = str(pathlib.Path(__file__).resolve().parent.parent
                / "data" / "earlywarning" / "fks1_calls")


def _args(**over):
    """A render_weather args namespace with the anchoring-relevant fields; sentinels default None."""
    base = dict(gene="fks1", snapshot=None, tsv=None, fks1_calls=None, as_of=None, window=None)
    base.update(over)
    return argparse.Namespace(**base)


# --- _resolve_window / display_window: the anchoring seam ------------------------------------

def test_display_window_anchors_to_committed_fill():
    # The committed release-710 fill is as_of 2026-04-14, window 180; display anchors to exactly it.
    as_of, window = fill_status.display_window(CALLS_DIR, "2099-01-01", 180)
    assert as_of == "2026-04-14"
    assert window == 180


def test_display_window_falls_back_when_no_fill():
    assert fill_status.display_window("/tmp/no_such_fill_dir", "2099-01-01", 99) == ("2099-01-01", 99)


def test_fks1_arm_self_anchors_independent_of_today():
    # The whole point: no --as-of/--window, fks1 live arm -> the fill's window, NOT today-180.
    args = _args()
    rw._resolve_window(args)
    assert args.as_of == "2026-04-14"
    assert args.window == 180


def test_explicit_flags_win_over_anchor():
    args = _args(as_of="2025-01-01", window=90)
    rw._resolve_window(args)
    assert args.as_of == "2025-01-01"
    assert args.window == 90


def test_snapshot_and_tsv_do_not_anchor():
    # A written snapshot / raw TSV carries its own calls, not the committed series -> plain today/180.
    import datetime
    today = datetime.date.today().isoformat()
    for src in ("snapshot", "tsv"):
        args = _args(**{src: "some.tsv"})
        rw._resolve_window(args)
        assert args.as_of == today
        assert args.window == 180


def test_erg11_arm_does_not_anchor():
    import datetime
    args = _args(gene="erg11")
    rw._resolve_window(args)
    assert args.as_of == datetime.date.today().isoformat()
    assert args.window == 180


# --- end-to-end: the page stays LIT past the fill's trailing-today lifetime -------------------

def test_page_renders_lit_without_as_of_regardless_of_today(monkeypatch, tmp_path, capsys):
    """Drive the real CLI with NO --as-of: it must anchor to the committed fill and come out lit
    (confident coverage clause, never 'awaiting re-call') -- the regression the fix prevents."""
    records = list(ew.read_fks1_call_table(ew.latest_fks1_call_table(CALLS_DIR)).values())
    monkeypatch.setattr(rw, "_load_records", lambda args: records)
    out = tmp_path / "index.html"
    monkeypatch.setattr(sys, "argv", ["render_weather.py", "--gene", "fks1", "--out", str(out)])
    rw.main()
    page = out.read_text()
    assert "window isolates re-called" in page      # a confident coverage %, i.e. lit
    assert "awaiting re-call" not in page
    assert "calls current through" in page and "2026-04-14" in page


# --- coverage denominator under a sliding today: anchored stays full, trailing decays to dark ---

def test_anchored_coverage_stays_full_while_trailing_window_goes_dark():
    import datetime
    records = list(ew.read_fks1_call_table(ew.latest_fks1_call_table(CALLS_DIR)).values())

    def coverage_at(as_of, window=180):
        _base, recent, _und = emergence.split_by_creation(records, as_of, window)
        resolved = len(emergence.resolved(recent, "fks1_call"))
        return resolved, len(recent)

    # Anchored window: the committed fill's own window -> every called isolate is in-window.
    anchor_as_of, anchor_window = fill_status.display_window(CALLS_DIR, "x", 180)
    a_called, a_total = coverage_at(anchor_as_of, anchor_window)
    assert a_called == a_total == len(records)      # 100% coverage, denominator == resolved count

    # Trailing-today window far past the fill (old behaviour): the fill has fully aged out -> dark.
    future = (datetime.date.fromisoformat("2027-05-01") - datetime.timedelta(days=180)).isoformat()
    t_called, t_total = coverage_at(future, 180)
    assert t_total == 0 and t_called == 0           # recent window no longer overlaps the fill


# --- VUS census: a token the geometry-free tool cannot place is DROPPED -----------------------

def test_vus_census_drops_unplaceable_token(tmp_path):
    """_vus_census lists only non-panel variants it can honestly situate: a token outside every
    hot-spot window (position_context -> None) is dropped; a placeable novel letter is kept."""
    table = tmp_path / "PDG000000067.999.tsv"
    header = "isolate_key\tfks1_call\tfks1_resistance_source\trun_acc\ttarget_creation_date\tcountry"
    rows = [
        # S639L: novel letter at a panel residue -> placeable, kept.
        "iso1\tS639L\tsra-fks1-recaller\tSRR1\t2026-05-01\tUSA",
        # A100V: outside every hot-spot window -> position_context None -> dropped.
        "iso2\tA100V\tsra-fks1-recaller\tSRR2\t2026-05-02\tUSA",
    ]
    table.write_text(header + "\n" + "\n".join(rows) + "\n")

    args = _args(fks1_calls=str(table))
    census = rw._vus_census(args)
    tokens = {e["token"] for e in (census or [])}
    assert "S639L" in tokens
    assert "A100V" not in tokens


# --- _geo_breakdown: reproducible current-release map from the committed table ----------------

def _panel_rec(key, geo_loc_name):
    """A minimal live-pull record carrying a panel-positive FKS1 call and finest geography."""
    return {"isolate_key": key, "fks1_call": "S639F", "geo_loc_name": geo_loc_name}


def test_geo_breakdown_reads_region_from_committed_table(tmp_path):
    # The committed table carries geo_loc_name -> the map is region-level and reproducible,
    # independent of whatever the (stand-in) live records say.
    table = tmp_path / "PDG000000067.999.tsv"
    header = ("isolate_key\tfks1_call\tfks1_resistance_source\trun_acc\t"
              "target_creation_date\tcountry\tgeo_loc_name")
    rows = [
        "iso1\tS639F\tsra-fks1-recaller:called\tSRR1\t2026-05-01\tUSA\tUSA: California",
        "iso2\tS639F\tsra-fks1-recaller:called\tSRR2\t2026-05-02\tUSA\tUSA: California",
        "iso3\tW691L\tsra-fks1-recaller:called\tSRR3\t2026-05-03\tIndia\tIndia",
        "iso4\tF635I\tsra-fks1-recaller:called\tSRR4\t2026-05-04\tUSA\tUSA: Texas",  # non-panel
    ]
    table.write_text(header + "\n" + "\n".join(rows) + "\n")
    args = _args(fks1_calls=str(table))
    # live records disagree; the committed table must win
    records = [_panel_rec("iso9", "Narnia")]
    assert rw._geo_breakdown(records, args) == {"USA: California": 2, "India": 1}


def test_geo_breakdown_falls_back_to_live_when_table_has_no_geography(tmp_path):
    # a pre-backfill table (no country/geo_loc_name) -> honest degradation: use the finest
    # geography the live records themselves carry, so the map never regresses to blank.
    table = tmp_path / "PDG000000067.998.tsv"
    header = ("isolate_key\tfks1_call\tfks1_resistance_source\trun_acc\t"
              "target_creation_date\tcountry\tgeo_loc_name")
    table.write_text(header + "\n"
                     + "iso1\tS639F\tsra-fks1-recaller:called\tSRR1\t2026-05-01\t\t\n")
    args = _args(fks1_calls=str(table))
    records = [_panel_rec("iso1", "USA: Oregon"), _panel_rec("iso2", "USA: Oregon")]
    assert rw._geo_breakdown(records, args) == {"USA: Oregon": 2}


def test_geo_breakdown_none_for_erg11_and_no_records(tmp_path):
    assert rw._geo_breakdown([], _args()) is None
    assert rw._geo_breakdown([_panel_rec("iso1", "USA")], _args(gene="erg11")) is None
