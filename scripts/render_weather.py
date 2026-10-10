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
from openafr import alert, weather  # noqa: E402
from openafr import earlywarning as ew  # noqa: E402

DIGEST_DIR = ROOT / "data" / "earlywarning" / "digest"
DEFAULT_STATE = DIGEST_DIR / "STATE.json"
FKS1_CALLS_DIR = ROOT / "data" / "earlywarning" / "fks1_calls"


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
    if args.fks1_calls:
        # Explicit table: read provenance from the sibling INDEX in its directory.
        calls_dir = str(pathlib.Path(args.fks1_calls).resolve().parent)
    else:
        calls_dir = str(FKS1_CALLS_DIR)
    return ew.latest_fks1_calls_index_entry(calls_dir)


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


def _compose(args):
    """Compose the current alert result; on no/failed data, an empty gene-marked result
    (-> the calm watching page for the selected arm)."""
    if args.watching_only:
        return _watching(args)
    try:
        records = _load_records(args)
    except Exception as exc:  # noqa: BLE001 -- a page must never fail on a data hiccup
        print(f"data unavailable ({exc!r}); rendering the watching page.", file=sys.stderr)
        return _watching(args)
    composer = (alert.compose_fks1_alerts if args.gene == "fks1"
                else alert.compose_alerts)
    return composer(
        records, args.as_of, window_days=args.window,
        min_count=args.min_count, min_delta=args.min_delta,
        backlog_frac=args.backlog_frac)


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
    ap.add_argument("--as-of", default=datetime.date.today().isoformat(),
                    help="baseline cut date, ISO (default: today)")
    ap.add_argument("--window", type=int, default=180, help="recent-window days")
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

    result = _compose(args)
    watching_since = _watching_since(args.state, args.watching_start)
    page = weather.render_page(result, watching_since=watching_since,
                               calls_provenance=_calls_provenance(args))

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page)
    n = len((result or {}).get("alerts") or [])
    print(f"wrote {out} — {n} alert(s), watching since {watching_since}.")


if __name__ == "__main__":
    main()
