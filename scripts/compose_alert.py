"""Alert-composition CLI for the C. auris azole early-warning track (issue #25).

Assembles the human-readable alert that closes the loop: emergence signal (#22) +
mutation->pocket mapping (#23) + structural fit so-what (#24), joined with the
regions and first-seen date the detector dropped, rendered as Markdown (for a human)
or JSON (for a machine). All method logic lives in openafr.alert; this is the CLI.

Two genes
---------
  --gene erg11 (default) -- azole/ERG11, with the #24 structural fit so-what.
  --gene fks1            -- echinocandin/FKS1 (v2), DETECTION-ONLY: emergence over the
                           `fks1_call` column with no structural verdict, because
                           echinocandins coordinate no metal, so the CYP51/heme-iron
                           geometry does not transfer (openafr/fks1_caller.py).

Subcommands
-----------
  # Compose alerts from a real snapshot in the store (#21):
  python scripts/compose_alert.py run \
      data/earlywarning/snapshots/PDG000000067.671.tsv --as-of 2025-01-01
  python scripts/compose_alert.py run <snapshot> --gene fks1 --as-of 2025-01-01

  # Prove the whole pipeline end-to-end on the synthetic #22 fixture:
  python scripts/compose_alert.py demo
  python scripts/compose_alert.py demo --gene fks1
  python scripts/compose_alert.py demo --json     # machine-readable form

On real snapshots today the erg11_call / fks1_call columns are empty for every isolate
(calls are pending the respective re-caller), so `run` honestly reports "no calls to
assess" rather than a false all-clear. `demo` runs a synthetic called snapshot so the
composed alert -- including the archival-backlog CONTEXT case -- is reproducible now.
"""
import argparse
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from openafr import alert  # noqa: E402
from openafr import earlywarning as ew  # noqa: E402


def _emit(result, as_json, gene):
    if as_json:
        print(alert.render_json(result))
    elif gene == "fks1":
        print(alert.render_fks1_markdown(result))
    else:
        print(alert.render_markdown(result))


def cmd_run(args):
    records = ew.read_snapshot(args.snapshot)
    if args.gene == "fks1":
        result = alert.compose_fks1_alerts(
            records, args.as_of, window_days=args.window, min_count=args.min_count,
            min_delta=args.min_delta, backlog_frac=args.backlog_frac)
    else:
        result = alert.compose_alerts(
            records, args.as_of, window_days=args.window, min_count=args.min_count,
            min_delta=args.min_delta, backlog_frac=args.backlog_frac)
    _emit(result, args.json, args.gene)


def cmd_demo(args):
    if args.gene == "fks1":
        result = alert.compose_fks1_alerts(_demo_fks1_records(), "2024-12-31",
                                           window_days=180, min_count=3, min_delta=0.05)
    else:
        sys.path.insert(0, str(ROOT / "scripts"))
        from detect_emergence import _demo_records  # reuse the #22 fixture
        result = alert.compose_alerts(_demo_records(), "2024-12-31", window_days=180,
                                      min_count=3, min_delta=0.05)
    _emit(result, args.json, args.gene)


def _demo_fks1_records():
    """A synthetic FKS1 snapshot mirroring the ERG11 #22 fixture, on `fks1_call`.

    Baseline (created 2024): S639Y endemic at ~50% of called isolates. Recent (2025):
    a genuine first-appearance S639F on freshly-collected isolates, plus an archival
    backlog spike of S639P whose carriers were all collected back in 2019. Exercises
    the watch (genuine) and context (backlog) tiers of the detection-only path.
    """
    recs = []

    def add(key, created, collected, call, country="Pakistan"):
        recs.append({
            "isolate_key": key, "target_acc": key + ".1",
            "target_creation_date": created, "collection_date": collected,
            "country": country, "fks1_call": call,
            "run_acc": "SRR" + key[-3:], "fks1_resistance_source": "recalled",
        })

    for i in range(10):
        add(f"PDFB{i:03d}", "2024-03-01", "2023", "S639Y" if i < 5 else "")
    for i in range(4):  # genuine first-appearance: S639F, freshly collected
        add(f"PDFR1{i:02d}", "2025-02-01", "2024", "S639F,S639Y")
    for i in range(5):  # archival backlog: S639P, all collected in 2019
        add(f"PDFR3{i:02d}", "2025-03-01", "2019", "S639P", country="India")
    add("PDFR400", "2025-02-15", "2024", "")  # uncalled, not susceptible
    return recs


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("run", help="compose alerts from a real snapshot")
    p.add_argument("snapshot")
    p.add_argument("--gene", choices=("erg11", "fks1"), default="erg11",
                   help="which re-caller's calls to assess: erg11 (azole, structural "
                        "so-what) or fks1 (echinocandin, detection-only)")
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
    p.add_argument("--gene", choices=("erg11", "fks1"), default="erg11",
                   help="erg11 (azole/structural) or fks1 (echinocandin/detection-only)")
    p.add_argument("--json", action="store_true", help="emit JSON instead of Markdown")
    p.set_defaults(func=cmd_demo)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
