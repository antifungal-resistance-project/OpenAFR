"""Emergence-detection CLI for the C. auris azole early-warning track (issue #22).

Diffs a snapshot's recent window against its baseline and flags ERG11 azole
substitutions that are newly-appearing or rising. All method logic lives in
openafr.emergence; this is the CLI + a synthetic demo.

Two subcommands
---------------
  # Run the detector on a real snapshot from the store (#21):
  python scripts/detect_emergence.py run \
      data/earlywarning/snapshots/PDG000000067.671.tsv --as-of 2025-01-01

  # Prove the method fires, without waiting on the ERG11 re-caller:
  python scripts/detect_emergence.py demo

On real snapshots today the erg11_call column is empty for every isolate (calls
are pending), so `run` will honestly report zero callable isolates rather than a
false all-clear. The `demo` builds a synthetic snapshot that *does* carry calls --
including a deliberate backlog spike -- so the detection and the backlog guard are
reproducible now and re-usable as a fixture when real calls arrive.
"""
import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import earlywarning as ew  # noqa: E402
from openafr import emergence as em  # noqa: E402


def _fmt_signal(s):
    tags = []
    if s["first_appearance"]:
        tags.append("FIRST-APPEARANCE")
    if s["rising"]:
        tags.append("RISING")
    if s["possible_backlog"]:
        tags.append("possible-backlog")
    if not s["canonical"]:
        tags.append("non-canonical-token")
    head = f"  {s['mutation']:10s} [{', '.join(tags)}]"
    body = (f"      baseline {s['baseline_count']}/{s['baseline_called']} "
            f"({s['baseline_prop']:.1%})  ->  recent {s['recent_count']}/"
            f"{s['recent_called']} ({s['recent_prop']:.1%})  "
            f"delta {s['delta']:+.1%}")
    why = f"      backlog: {s['backlog_reason']}"
    return "\n".join([head, body, why])


def _report(result):
    c = result["coverage"]
    print(f"as-of {result['as_of']}  window {result['window_days']}d  "
          f"thresholds {result['thresholds']}")
    print(f"coverage: baseline {c['baseline_called']}/{c['baseline_total']} called, "
          f"recent {c['recent_called']}/{c['recent_total']} called, "
          f"{c['undated']} undated (unplaceable)")
    if "note" in result:
        print(f"\nNOTE: {result['note']}")
    if not result["signals"]:
        print("\nno emergence signals over the configured thresholds.")
        return
    genuine = [s for s in result["signals"] if not s["possible_backlog"]]
    print(f"\n{len(result['signals'])} signal(s) "
          f"({len(genuine)} genuine, {len(result['signals']) - len(genuine)} "
          f"possible-backlog):\n")
    for s in result["signals"]:
        print(_fmt_signal(s))


def cmd_run(args):
    records = ew.read_snapshot(args.snapshot)
    result = em.detect_emergence(
        records, args.as_of, window_days=args.window,
        min_count=args.min_count, min_delta=args.min_delta,
        backlog_frac=args.backlog_frac)
    print(f"snapshot {pathlib.Path(args.snapshot).name}: {len(records)} isolates")
    _report(result)


def _demo_records():
    """A synthetic snapshot exercising every branch of the detector.

    Baseline (created 2024): Y132F endemic at ~50% of called isolates, K143R rare.
    Recent (created 2025): a genuine new substitution F126L on freshly-collected
    isolates, K143R climbing, and a fake spike of TR34 whose carriers were all
    *collected* in 2019 -- an archival backlog the guard must catch.
    """
    recs = []

    def add(key, created, collected, call, country="Pakistan"):
        recs.append({
            "isolate_key": key, "target_acc": key + ".1",
            "target_creation_date": created, "collection_date": collected,
            "country": country, "erg11_call": call,
            "run_acc": "SRR" + key[-3:], "resistance_source": "recalled",
        })

    # --- baseline window (created 2024) : 10 called isolates ---
    for i in range(10):
        call = "Y132F" if i < 5 else ("K143R" if i == 5 else "")
        add(f"PDTB{i:03d}", "2024-03-01", "2023", call)

    # --- recent window (created 2025) ---
    # genuine first-appearance: F126L on 4 freshly-collected isolates
    for i in range(4):
        add(f"PDTR1{i:02d}", "2025-02-01", "2024", "F126L,Y132F")
    # genuine rising: K143R jumps from 1/6 baseline to 4 recent carriers
    for i in range(4):
        add(f"PDTR2{i:02d}", "2025-02-10", "2024", "K143R,Y132F")
    # archival backlog: TR34 spike, but every carrier was collected back in 2019
    for i in range(5):
        add(f"PDTR3{i:02d}", "2025-03-01", "2019", "TR34,Y132F", country="India")
    # a couple of uncalled recent isolates -- must not be read as susceptible
    add("PDTR400", "2025-02-15", "2024", "")
    add("PDTR401", "2025-02-16", "2024", "")
    return recs


def cmd_demo(args):
    records = _demo_records()
    print("SYNTHETIC demo snapshot (calls are pending on real data, #23):")
    print(f"  {len(records)} isolates, {len(em.called(records))} called\n")
    result = em.detect_emergence(records, "2024-12-31", window_days=180,
                                 min_count=3, min_delta=0.05)
    _report(result)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("run", help="detect emergence in a real snapshot")
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
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("demo", help="run the method on a synthetic called snapshot")
    p.set_defaults(func=cmd_demo)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
