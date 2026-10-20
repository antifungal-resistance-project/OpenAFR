"""Render the C. auris resistance weather report to a static HTML page (issue: weather).

The page is a thin VIEW over the same early-warning output deliver_alert.py acts on --
it composes the current alert result and hands it to openafr.weather.render_page. It
pulls nothing new of its own beyond the record source the composer needs, and it NEVER
crashes the build: on a failed live pull or an empty snapshot it renders the honest
day-0 "watching" page instead of erroring, because a page that 500s on a quiet week is
worse than a page that says "quiet week".

The pocket image is not rendered here (PyMOL is out of the offline path). The page
references pockets/<token>.png; the workflow copies committed assets next to index.html.

The arm defaults to FKS1/echinocandin -- the surveillance product. Pass --gene erg11 to
render the azole/ERG11 research arm instead.

Usage
-----
  # Live pull, render the echinocandin page to the default site dir:
  python scripts/render_weather.py --as-of 2025-01-01 --out site/index.html

  python scripts/render_weather.py --snapshot data/.../PDG..tsv --as-of ... --out site/index.html
  python scripts/render_weather.py --gene erg11 --snapshot data/.../PDG..tsv --as-of ... --out ...
  python scripts/render_weather.py --tsv saved.tsv --as-of ... --out site/index.html
  python scripts/render_weather.py --watching-only --out site/index.html   # skip data, day-0 page
"""
import argparse
import csv
import datetime
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from openafr import alert, weather, weather_feed  # noqa: E402
from openafr import earlywarning as ew  # noqa: E402
from openafr import fill_status  # noqa: E402

DIGEST_DIR = ROOT / "data" / "earlywarning" / "digest"
DEFAULT_STATE = DIGEST_DIR / "STATE.json"
FKS1_CALLS_DIR = ROOT / "data" / "earlywarning" / "fks1_calls"
CONCORDANCE_RUNLOG = ROOT / "data" / "earlywarning" / "runlog" / "concordance-fks1.jsonl"
# The public page's canonical URL, for the shareable OpenGraph/Twitter card.
PAGE_URL = "https://antifungal-resistance-project.github.io/OpenAFR/"


def _overlay_live_calls(records, args):
    """Replay the durable, git-tracked FKS1 calls onto a freshly pulled set (issue #171).

    A live NCBI pull (or a raw --tsv) carries no echinocandin call -- every row is the honest
    `pending:` state. The billable re-caller output lives in a committed call table; overlay it
    so the page shows the real echinocandin picture instead of day-0. No-op for the erg11 arm
    or when no table is found. A written --snapshot is left untouched (it carries its own calls).
    """
    if args.gene != "fks1":
        return
    table = args.fks1_calls or ew.latest_fks1_call_table(str(FKS1_CALLS_DIR))
    if not table or not pathlib.Path(table).exists():
        return
    n = ew.overlay_fks1_calls(records, str(table))
    print(f"overlaid {n} FKS1 call(s) from {table}", file=sys.stderr)


def _calls_provenance(args):
    """Provenance of the FKS1 call table the page overlays (`as_of`, `n_called`), for the
    freshness clause -- or None for the erg11 arm or when no fill exists. Mirrors the table
    auto-discovery in `_overlay_live_calls`, so what the clause reports is what was overlaid.
    A written --snapshot carries its own calls (no overlay), so its freshness is its own
    as_of, not the committed table's -- return None there rather than claim the fill date."""
    if args.gene != "fks1" or args.snapshot:
        return None
    return ew.latest_fks1_calls_index_entry(_calls_dir(args))


def _calls_dir(args):
    """The FKS1-calls directory whose INDEX/tables back the trend + freshness, mirroring the
    overlay auto-discovery so what the page reports is what it overlaid."""
    if args.fks1_calls:
        return str(pathlib.Path(args.fks1_calls).resolve().parent)
    return str(FKS1_CALLS_DIR)


def _release_trend(args):
    """The per-release resistance trajectory the trend section renders, or None.

    Fully reproducible from committed data: one entry per committed fill, read from the calls
    INDEX (`release_tag`, `as_of`, resolved `n_called`) joined to that release's committed
    table for the resistant-carrier count. `n_resistant` is the KNOWN-resistance panel count
    (not every FKS1 variant token): non-panel hot-spot variants are of uncertain significance
    and are carried separately as `n_vus`, so the headline never overstates resistance (the
    distinction the real-data backtest surfaced). Prevalence is panel-resistant / resolved --
    the same resolved denominator the coverage fix uses. None for the erg11 arm or a written
    --snapshot (which carries its own calls, not the committed series)."""
    if args.gene != "fks1" or args.snapshot:
        return None
    calls_dir = _calls_dir(args)
    trend = []
    for row in ew.read_fks1_calls_index(calls_dir):
        tag = row.get("release_tag") or ""
        table = pathlib.Path(calls_dir) / f"{tag}.tsv"
        n_called = int(row.get("n_called") or 0)
        n_resistant, n_vus = ew.count_fks1_panel_positive(str(table))
        trend.append({
            "release_tag": tag,
            "as_of": row.get("as_of") or "",
            "n_called": n_called,
            "n_resistant": n_resistant,
            "n_vus": n_vus,
            "prevalence": (n_resistant / n_called) if n_called else None,
        })
    return trend or None


def _forecast(args):
    """The prevalence trajectory + near-term forecast bundle the page/feed surface, or None.

    Reproducible from the SAME latest committed call table the page overlays: build the
    resolved-denominator, panel-gated per-deposit-month series (openafr.forecast), run the
    pre-registered walk-forward skill test, and only mark `skill=True` when the trend model
    beat both naive baselines on held-out months. The page renders a projection only on skill;
    otherwise it shows the trajectory and says no skill was demonstrated -- the honest
    degradation the backtest enforces. None for the erg11 arm or a written --snapshot (which
    carries its own calls, not the committed series)."""
    if args.gene != "fks1" or args.snapshot:
        return None
    from openafr import forecast
    table = args.fks1_calls or ew.latest_fks1_call_table(_calls_dir(args))
    if not table or not pathlib.Path(table).exists():
        return None
    records = list(ew.read_fks1_call_table(str(table)).values())
    series = forecast.monthly_prevalence_series(records)
    if len(series) < 2:
        return None
    bt = forecast.backtest(series, min_train=3, horizon=1)
    return {
        "series": series,
        "skill": forecast.has_skill(bt),
        "n_folds": bt["n_folds"],
        "model_mae": bt["model_mae"],
        "naive_last_mae": bt["naive_last_mae"],
        "naive_pooled_mae": bt["naive_pooled_mae"],
        "projection": forecast.forecast_next(series, horizon=1),
    }


def _vus_census(args):
    """The per-variant census of non-panel hot-spot FKS1 variants (uncertain significance) the
    page/feed surface, each enriched with geometry-free position context, or None.

    Reads the SAME latest committed table the page overlays (`fks1_vus_census`), so the census
    is reproducible from committed data alone -- then attaches `position_context` + a
    pre-rendered geography sentence (openafr/fks1_position) exactly as the alert composer does
    (openafr/alert.py), keeping the view decoupled from the geometry derivation. A token the
    tool cannot place (outside every hot-spot window, or malformed) is dropped, so the census
    only lists variants it can honestly situate. None for the erg11 arm or a written --snapshot
    (which carries its own calls, not the committed series)."""
    if args.gene != "fks1" or args.snapshot:
        return None
    from openafr import fks1_position
    table = args.fks1_calls or ew.latest_fks1_call_table(_calls_dir(args))
    if not table or not pathlib.Path(table).exists():
        return None
    out = []
    for entry in ew.fks1_vus_census(table):
        ctx = fks1_position.position_context(entry["token"])
        if ctx is None:
            continue
        out.append({
            "token": entry["token"],
            "carriers": entry["carriers"],
            "position_context": ctx,
            "position_note": f"It {fks1_position.describe(ctx)}",
        })
    return out or None


def _geo_breakdown(records, args):
    """Where the current release's resolved carriers are: {location: carrier_count}, or None.

    Zero-cost -- it reads the finest geography the live NCBI pull already carries on each
    record (ew.display_geo prefers the raw geo_loc_name's region over the bare country);
    counts only resolved carriers of a KNOWN echinocandin-resistance panel substitution (not
    every FKS1 variant token), matching the panel-gated headline so the map never overstates
    resistance. Reflects the current release only, since the committed table does not yet carry
    geography -- the section labels it as such. None for the erg11 arm, no records, or no located
    panel carriers."""
    if args.gene != "fks1" or not records:
        return None
    from openafr import emergence, fks1_caller
    counts = {}
    for rec in records:
        toks = emergence.parse_call(rec.get("fks1_call"))
        if not any(fks1_caller.is_panel_token(t) for t in toks):
            continue
        loc = ew.display_geo(rec)
        if loc:
            counts[loc] = counts.get(loc, 0) + 1
    return counts or None


def _geo_history(args):
    """Per-release historical geography the page/feed surface, oldest release first, or None.

    The historical counterpart to `_geo_breakdown` (current release only): walks the committed
    calls INDEX (same source as `_release_trend`) and, for each release table, reads the
    panel-positive resistant carriers by `country` straight from the committed table
    (`ew.fks1_calls_geo_breakdown` -- same panel gate as the headline, so geography never
    overstates resistance). A release whose table predates the `country` column (709/710 before
    backfill) contributes an empty breakdown and is dropped, so the section only shows releases
    it can honestly place -- honest degradation, same as the live geo. None for the erg11 arm or
    a written --snapshot (which carries its own calls, not the committed series)."""
    if args.gene != "fks1" or args.snapshot:
        return None
    calls_dir = _calls_dir(args)
    out = []
    for row in ew.read_fks1_calls_index(calls_dir):
        tag = row.get("release_tag") or ""
        table = pathlib.Path(calls_dir) / f"{tag}.tsv"
        breakdown = ew.fks1_calls_geo_breakdown(str(table))
        if not breakdown:
            continue
        out.append({"release_tag": tag, "as_of": row.get("as_of") or "",
                    "geography": breakdown})
    return out or None


def _certification(runlog_path=CONCORDANCE_RUNLOG):
    """Certification verdict for the FKS1 caller, scored live from the committed concordance
    runlog through the SAME scorer (`openafr.concordance.evaluate` + `apply_bar`) the CI
    replay gate asserts against on every PR (#164/#179). The badge is therefore a computed,
    self-updating fact -- never a hardcoded claim. Returns {verdict, n, n_resolved,
    sensitivity, specificity, identity, benchmark} or None on any failure (fail-open: the
    page just omits the badge, the same page-never-fails discipline as `_watching`)."""
    try:
        from openafr import concordance
        lines = [ln for ln in pathlib.Path(runlog_path).read_text().splitlines() if ln.strip()]
        record = json.loads(lines[-1])
        items = record.get("items") or []
        if not items:
            return None
        summary = concordance.evaluate(items)
        verdict, _reasons = concordance.apply_bar(summary)

        def _pt(metric):
            p = summary[metric]["point"]
            return None if p is None else round(p, 4)

        return {
            "verdict": verdict,
            "n": summary["n_total"],
            "n_resolved": summary["n_resolved"],
            "sensitivity": _pt("sensitivity"),
            "specificity": _pt("specificity"),
            "identity": _pt("identity"),
            "benchmark": record.get("fixture") or "concordance panel",
        }
    except Exception as exc:  # noqa: BLE001 -- the badge must never fail the build
        print(f"certification unavailable ({exc!r}); omitting the badge.", file=sys.stderr)
        return None


def _watching_since(state_path, fallback):
    """The date to frame 'quiet since' by: the delivery state's `updated`, else fallback.

    STATE.json is `{"delivered": {...}, "updated": "<iso>"}`; an empty string (never
    delivered) means we are genuinely at day 0, so we fall back to the configured start.
    """
    try:
        state = json.loads(pathlib.Path(state_path).read_text())
        return state.get("updated") or fallback
    except (OSError, ValueError):
        return fallback


def _load_records(args):
    """Resolve the record source, mirroring deliver_alert.py. Raises on failure so the
    caller can fall back to the watching page."""
    if args.snapshot:
        return ew.read_snapshot(args.snapshot)
    if args.tsv:
        with open(args.tsv, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        records, _ = ew.normalize(rows)
        _overlay_live_calls(records, args)
        return records
    rows, _tag, _url = ew.fetch_release_rows(args.release)
    records, _ = ew.normalize(rows)
    _overlay_live_calls(records, args)
    return records


def _watching(args):
    """An empty, gene-marked result so the calm page names the right arm even when no
    data is loaded (watching-only, or a data hiccup) -- render_page reads `gene` for copy."""
    return {"gene": "FKS1" if args.gene == "fks1" else "ERG11", "alerts": []}


def _resolve_window(args):
    """Fill in the --as-of / --window sentinels IN PLACE.

    The live fks1 page shows a frozen, git-tracked fill; a trailing-`today` window would slowly age
    the fill's isolates out and silently decay the page toward "awaiting re-call". So when the user
    did NOT pin --as-of/--window, the fks1 arm anchors the displayed window to the committed fill's
    own window (fill_status.display_window), keeping the page lit at full coverage with an honest
    "calls current through <fill date>" clause. A written --snapshot / raw --tsv carries its own
    calls (not the committed series), so those are left on the plain today/180 default -- mirroring
    the snapshot guards in _calls_provenance / _release_trend. An explicit flag always wins."""
    today = datetime.date.today().isoformat()
    anchor = (args.gene == "fks1" and not args.snapshot and not args.tsv)
    if anchor and (args.as_of is None or args.window is None):
        as_of, window = fill_status.display_window(_calls_dir(args), today, 180)
        if args.as_of is None:
            args.as_of = as_of
        if args.window is None:
            args.window = window
    if args.as_of is None:
        args.as_of = today
    if args.window is None:
        args.window = 180


def _compose(args):
    """Compose the current alert result and return (result, records); on no/failed data, an
    empty gene-marked result with records=None (-> the calm watching page for the selected
    arm). The records ride along so the geo breakdown can aggregate them without a second load.
    """
    if args.watching_only:
        return _watching(args), None
    try:
        records = _load_records(args)
    except Exception as exc:  # noqa: BLE001 -- a page must never fail on a data hiccup
        print(f"data unavailable ({exc!r}); rendering the watching page.", file=sys.stderr)
        return _watching(args), None
    composer = (alert.compose_fks1_alerts if args.gene == "fks1"
                else alert.compose_alerts)
    result = composer(
        records, args.as_of, window_days=args.window,
        min_count=args.min_count, min_delta=args.min_delta,
        backlog_frac=args.backlog_frac)
    return result, records


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--snapshot", help="read an existing written snapshot TSV")
    src.add_argument("--tsv", help="offline raw NCBI metadata TSV")
    src.add_argument("--watching-only", action="store_true",
                     help="skip data entirely; render the day-0 watching page")
    ap.add_argument("--gene", choices=("fks1", "erg11"), default="fks1",
                    help="which arm to render: fks1/echinocandin (the product, default) "
                         "or erg11/azole (the research arm)")
    ap.add_argument("--release", type=int, default=None,
                    help="pin an NCBI release to pull (default: latest)")
    ap.add_argument("--fks1-calls", default=None,
                    help="committed FKS1 call table to overlay onto a live pull / --tsv "
                         "(default: newest in data/earlywarning/fks1_calls; fks1 arm only)")
    ap.add_argument("--as-of", default=None,
                    help="baseline cut date, ISO (default: the fks1 arm anchors to the latest "
                         "committed fill's window; otherwise today)")
    ap.add_argument("--window", type=int, default=None,
                    help="recent-window days (default: the fks1 arm anchors to the committed "
                         "fill's window; otherwise 180)")
    ap.add_argument("--min-count", type=int, default=3, help="min recent carriers to flag")
    ap.add_argument("--min-delta", type=float, default=0.05, help="min share increase for RISING")
    ap.add_argument("--backlog-frac", type=float, default=0.5,
                    help="archival-carrier fraction that trips the backlog guard")
    ap.add_argument("--state", default=str(DEFAULT_STATE),
                    help="committed delivery-state path (for the 'quiet since' framing)")
    ap.add_argument("--watching-start", default="2026-09-17",
                    help="fallback watch-start date when nothing has been delivered yet")
    ap.add_argument("--out", default="site/index.html", help="output HTML path")
    args = ap.parse_args()
    _resolve_window(args)

    result, records = _compose(args)
    watching_since = _watching_since(args.state, args.watching_start)
    # Certification is the caller's, not the fill's -- the azole research arm has no such
    # concordance record, so only the echinocandin product surfaces the badge.
    certification = _certification() if args.gene == "fks1" else None
    shared = dict(calls_provenance=_calls_provenance(args),
                  certification=certification,
                  release_trend=_release_trend(args),
                  geo_breakdown=_geo_breakdown(records, args),
                  geo_history=_geo_history(args),
                  vus_census=_vus_census(args),
                  forecast=_forecast(args),
                  page_url=PAGE_URL)
    # One clock reading so the page and both feeds agree on "generated".
    now = datetime.datetime.now(datetime.timezone.utc)
    page = weather.render_page(result, now=now, watching_since=watching_since,
                               feeds=True, **shared)

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page)
    # Sibling machine-readable feeds (same inputs as the page, so they cannot disagree).
    (out.parent / "data.json").write_text(json.dumps(
        weather_feed.build_data(result, now=now, **shared), indent=2, sort_keys=True) + "\n")
    (out.parent / "alerts.atom").write_text(
        weather_feed.render_atom(result, now=now, page_url=PAGE_URL))
    n = len((result or {}).get("alerts") or [])
    print(f"wrote {out} — {n} alert(s), watching since {watching_since}.")


if __name__ == "__main__":
    main()
