"""Backtest / honest-validation CLI for the C. auris azole early-warning track (#27).

The credibility test for the whole track. All logic lives in openafr.backtest and
openafr.emergence; this is the CLI + a labelled (synthetic) mechanism replay.

Three subcommands, one per pre-registered part (work/PREREGISTRATION_backtest.md):

  # H1-arrival -- the REAL lead-time budget (deposit lag) on a real snapshot:
  python scripts/backtest_earlywarning.py lag \
      data/earlywarning/snapshots/PDG000000067.671.tsv

  # H2-null -- the honest walk-forward null on the REAL snapshot (calls are empty):
  python scripts/backtest_earlywarning.py replay \
      data/earlywarning/snapshots/PDG000000067.671.tsv

  # H3-mechanism -- lead time on a LABELLED (synthetic) replay, proving the detector
  # turns a first-appearance into a dated, positive-lead-time flag once calls exist:
  python scripts/backtest_earlywarning.py mechanism

  # Azole event frequency -- the deliverable number once `fill` has populated calls:
  python scripts/backtest_earlywarning.py prevalence \
      data/earlywarning/snapshots/PDG000000067.673.tsv

On real snapshots `erg11_call` is empty for every isolate (the ERG11 re-caller was
never built), so `replay` honestly reports "cannot detect yet" at every step rather
than a false all-clear. `lag` needs no calls -- it measures the arrival budget the
warning would have to work within, which is real and measurable today.
"""
import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import backtest as bt  # noqa: E402
from openafr import earlywarning as ew  # noqa: E402


# ------------------------------- H1: lag --------------------------------------

def cmd_lag(args):
    records = ew.read_snapshot(args.snapshot)
    s = bt.lag_summary(records, within_days=args.within)
    print(f"snapshot {pathlib.Path(args.snapshot).name}: {len(records)} isolates")
    print("\nH1-arrival -- deposit lag (collection_date -> target_creation_date)")
    print(f"  datable {s['n_datable']}/{s['n_total']} "
          f"({s['n_undatable']} undatable, excluded not guessed)")
    if s["median"] is None:
        print("  no isolate has both a collection and a visibility day -- cannot decide.")
        return
    print(f"  lag days: median {s['median']:.0f}  "
          f"p25 {s['p25']:.0f}  p75 {s['p75']:.0f}  p90 {s['p90']:.0f}  "
          f"[min {s['min']}, max {s['max']}]")
    print(f"  visible within {s['within_days']}d of collection: "
          f"{s['n_within']}/{s['n_datable']} ({s['frac_within']:.0%})")
    if s["n_negative"]:
        print(f"  note: {s['n_negative']} isolate(s) recorded as collected AFTER "
              f"becoming visible (surfaced, not clamped)")
    passed, reason = bt.arrival_verdict(s)
    verdict = {True: "PASS-arrival", False: "FAIL-arrival", None: "UNDECIDED"}[passed]
    print(f"\n  pre-registered verdict: {verdict}\n    {reason}")


# ------------------------------- H2: replay -----------------------------------

def cmd_replay(args):
    records = ew.read_snapshot(args.snapshot)
    start, stop = bt.creation_date_range(records)
    steps = bt.month_steps(start, stop)
    print(f"snapshot {pathlib.Path(args.snapshot).name}: {len(records)} isolates")
    print(f"\nH2-null -- walk-forward over real creation dates "
          f"{start} .. {stop} ({len(steps)} monthly steps), target {args.target}")
    wf = bt.walk_forward(records, steps, target=args.target,
                         window_days=args.window, min_count=args.min_count,
                         min_delta=args.min_delta, backlog_frac=args.backlog_frac)
    n_callable = sum(1 for st in wf["steps"] if st["callable"])
    print(f"  steps with >=1 called isolate in the recent window: "
          f"{n_callable}/{len(wf['steps'])}")
    print(f"  steps that flagged {args.target}: "
          f"{sum(1 for st in wf['steps'] if st['target_flagged'])}")
    if not wf["any_callable"]:
        print("\n  HONEST NULL: no step had a single called isolate, because "
              "erg11_call is empty on real data (the ERG11 re-caller was never "
              "built). Every step returned 'cannot detect yet' -- the pipeline "
              "refusing to emit a false all-clear, NOT a validated warning.")
    if wf["first_flag"]:
        print(f"  first flag of {args.target}: {wf['first_flag']}")


# -------------------------- azole event frequency -----------------------------

def cmd_prevalence(args):
    """The deliverable number: azole event frequency on a re-called snapshot (TODOS step 2).

    On a snapshot the re-caller has NOT touched, every row is pending, so this reports 0
    resolved and an undefined frequency -- the honest 'not measured yet', not a zero signal.
    After `recall_erg11.py fill`, it becomes the measured number the FKS1/v2 decision gates on.
    """
    records = ew.read_snapshot(args.snapshot)
    p = bt.panel_prevalence(records)
    print(f"snapshot {pathlib.Path(args.snapshot).name}: {p['n_total']} isolates")
    print("\nAzole event frequency -- resolved ERG11 panel carrying a known resistance mutation")
    ex = p["excluded"]
    print(f"  resolved panel:   {p['n_resolved']}/{p['n_total']} "
          f"(excluded: {ex['partial']} partial, {ex['failed']} failed/refused, "
          f"{ex['pending']} pending -- not counted azole-negative)")
    if p["event_frequency"] is None:
        print("\n  NOT MEASURED YET: no isolate has a resolved panel (run "
              "`recall_erg11.py fill` first). This is 'cannot measure', NOT a zero signal.")
        return
    print(f"  panel-positive:   {p['n_panel_positive']}/{p['n_resolved']} "
          f"= {p['event_frequency']:.1%}  (azole event frequency)")
    if p["n_called_nonpanel"]:
        print(f"  called non-panel: {p['n_called_nonpanel']} resolved isolate(s) carry only "
              f"non-panel substitution(s) -- azole-negative, emergence.py judges novelty")
    if p["per_mutation"]:
        print("  by panel mutation (resolved isolates carrying it):")
        for tok, n in p["per_mutation"].items():
            print(f"    {tok}: {n}  ({n / p['n_resolved']:.1%})")


# ---------------------------- H3: mechanism -----------------------------------

def _labelled_records():
    """A LABELLED (synthetic) replay cohort in which Y132F first appears.

    Illustrative dates, exactly as detect_emergence's demo is synthetic -- the real
    timeline needs the ERG11 re-caller applied to historical SRA reads, which does
    not exist. Baseline (2022): a small called background, no Y132F. Then Y132F is
    first COLLECTED mid-2023 and becomes NCBI-VISIBLE after a ~5-month deposit lag,
    on three isolates across late-2023/early-2024 -- a first-appearance the detector
    must catch well before the reference recognition date.
    """
    recs = []

    def add(key, created, collected, call, country="Pakistan"):
        recs.append({
            "isolate_key": key, "target_acc": key + ".1",
            "target_creation_date": created, "collection_date": collected,
            "country": country, "erg11_call": call,
            "run_acc": "SRR" + key[-3:], "resistance_source": "recalled",
        })

    # baseline background, visible through 2022: endemic K143R, no Y132F
    for i in range(6):
        call = "K143R" if i < 3 else ""
        add(f"PDTB{i:03d}", f"2022-0{3 + i}-01", "2021", call)

    # Y132F first appearance: collected mid-2023, visible ~5 months later
    add("PDTY001", "2023-11-05", "2023-06", "Y132F", country="India")
    add("PDTY002", "2023-12-10", "2023-07", "Y132F", country="India")
    add("PDTY003", "2024-01-15", "2023-08", "Y132F", country="Kenya")
    return recs


def cmd_mechanism(args):
    records = _labelled_records()
    start, stop = bt.creation_date_range(records)
    steps = bt.month_steps(start, stop)
    print("LABELLED (synthetic) mechanism replay -- illustrative dates, NOT real "
          "calls (#27 H3):")
    print(f"  {len(records)} isolates, creation range {start} .. {stop}, "
          f"reference recognition date {args.reference}")
    wf = bt.walk_forward(records, steps, target="Y132F",
                         window_days=args.window, min_count=args.min_count,
                         min_delta=args.min_delta, backlog_frac=args.backlog_frac)
    print(f"\n  walk-forward: {len(wf['steps'])} monthly steps, "
          f"{sum(1 for st in wf['steps'] if st['target_flagged'])} flagged Y132F")
    if not wf["first_flag"]:
        print("  Y132F never flagged -- FAIL-mechanism (unexpected).")
        return
    lead = bt.lead_time_days(wf["first_flag"], args.reference)
    print(f"  first flag of Y132F: {wf['first_flag']}")
    print(f"  reference recognition: {args.reference}")
    passed = lead is not None and lead > 0
    verdict = "PASS-mechanism" if passed else "FAIL-mechanism"
    if lead is None:
        print(f"\n  {verdict}: lead time undefined (unparseable date).")
    else:
        print(f"  lead time: {lead} days ({lead / 30.0:.1f} months) "
              f"{'BEFORE' if lead > 0 else 'AFTER'} reference")
        print(f"\n  {verdict}: the detector converts a first-appearance into a dated, "
              f"positive-lead-time flag once calls exist. This is a MECHANISM proof, "
              f"not real-data evidence -- it isolates H2's null to the missing re-caller.")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("lag", help="H1: deposit-lag budget on a real snapshot")
    p.add_argument("snapshot")
    p.add_argument("--within", type=int, default=365,
                   help="lead-time window in days for the frac-within metric (default 365)")
    p.set_defaults(func=cmd_lag)

    p = sub.add_parser("replay", help="H2: honest walk-forward null on a real snapshot")
    p.add_argument("snapshot")
    p.add_argument("--target", default="Y132F", help="mutation to track (default Y132F)")
    p.add_argument("--window", type=int, default=180)
    p.add_argument("--min-count", type=int, default=3)
    p.add_argument("--min-delta", type=float, default=0.05)
    p.add_argument("--backlog-frac", type=float, default=0.5)
    p.set_defaults(func=cmd_replay)

    p = sub.add_parser("prevalence",
                       help="azole event frequency on a re-called snapshot (the deliverable)")
    p.add_argument("snapshot")
    p.set_defaults(func=cmd_prevalence)

    p = sub.add_parser("mechanism", help="H3: lead time on a labelled (synthetic) replay")
    p.add_argument("--reference", default="2024-09-01",
                   help="stand-in clinical/published recognition date (default 2024-09-01)")
    p.add_argument("--window", type=int, default=180)
    p.add_argument("--min-count", type=int, default=3)
    p.add_argument("--min-delta", type=float, default=0.05)
    p.add_argument("--backlog-frac", type=float, default=0.5)
    p.set_defaults(func=cmd_mechanism)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
