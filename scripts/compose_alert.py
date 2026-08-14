"""Alert-composition CLI for the C. auris azole early-warning track (issue #25).

Assembles the human-readable alert that closes the loop: emergence signal (#22) +
mutation->pocket mapping (#23) + structural fit so-what (#24), joined with the
regions and first-seen date the detector dropped, rendered as Markdown (for a human)
or JSON (for a machine). All method logic lives in openafr.alert; this is the CLI.

Three subcommands
-----------------
  # Compose alerts from a real snapshot in the store (#21):
  python scripts/compose_alert.py run \
      data/earlywarning/snapshots/PDG000000067.671.tsv --as-of 2025-01-01

  # Prove the whole pipeline end-to-end on the synthetic #22 fixture:
  python scripts/compose_alert.py demo
  python scripts/compose_alert.py demo --json     # machine-readable form

On real snapshots today the erg11_call column is empty for every isolate (calls are
pending the ERG11 re-caller, #23), so `run` honestly reports "no calls to assess"
rather than a false all-clear. `demo` runs the synthetic called snapshot so the
composed alert -- including the archival-backlog CONTEXT case -- is reproducible now.
"""
import argparse
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from openafr import alert  # noqa: E402
from openafr import earlywarning as ew  # noqa: E402


def _emit(result, as_json):
    print(alert.render_json(result) if as_json else alert.render_markdown(result))


def cmd_run(args):
    records = ew.read_snapshot(args.snapshot)
    result = alert.compose_alerts(
        records, args.as_of, window_days=args.window, min_count=args.min_count,
        min_delta=args.min_delta, backlog_frac=args.backlog_frac)
    _emit(result, args.json)


def cmd_demo(args):
    sys.path.insert(0, str(ROOT / "scripts"))
    from detect_emergence import _demo_records  # reuse the #22 fixture
    result = alert.compose_alerts(_demo_records(), "2024-12-31", window_days=180,
                                  min_count=3, min_delta=0.05)
    _emit(result, args.json)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("run", help="compose alerts from a real snapshot")
    p.add_argument("snapshot")
    p.add_argument("--as-of", required=True, help="baseline cut date, ISO")
    p.add_argument("--window", type=int, default=180,
                   help="recent-window length in days (default 180)")
    p.add_argument("--min-count", type=int, default=3,
                   help="min recent carriers to flag (default 3)")
    p.add_argument("--min-delta", type=float, default=0.05,
                   help="min recent-vs-baseline share increase to call RISING")
    p.add_argument("--backlog-frac", type=float, default=0.5,
                   help="archival-carrier fraction that trips the backlog guard")
    p.add_argument("--json", action="store_true", help="emit JSON instead of Markdown")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("demo", help="compose alerts from the synthetic called snapshot")
    p.add_argument("--json", action="store_true", help="emit JSON instead of Markdown")
    p.set_defaults(func=cmd_demo)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
