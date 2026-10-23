"""Is the billable FKS1 fill due? — the non-billable watchdog behind the weather cadence.

Pulls (or reads) a fresh *C. auris* record set, overlays the newest committed FKS1 call table,
and measures the in-window re-call gap via `openafr.fill_status.fill_due`. When the gap crosses a
threshold (or the last fill's whole window has rolled off the page) it writes a GitHub-issue body
and sets a `due=true` step-output so a scheduled workflow can ring the bell for a human to run the
billable GCP fill (`work/RUNBOOK_fks1_fill.md`). Otherwise it is a quiet no-op — the weekly page
publish keeps showing the last fill honestly, with its freshness clause, until a refill is worth
paying for.

This script does NO billable work: the only network it touches is the free NCBI metadata pull
(the same `fetch_release_rows` the early-warning schedule already runs in CI); it never downloads
SRA reads or runs the re-caller. It is the operational half that turns a one-off snapshot into a
self-driving service — up to the human-gated paid step, which it deliberately does not automate.

Mirrors scripts/deliver_alert.py's record-source + GitHub-output handshake so the workflow stays a
thin shell.
"""
import argparse
import csv
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from openafr import earlywarning as ew  # noqa: E402
from openafr import emergence            # noqa: E402
from openafr import fill_status          # noqa: E402

FKS1_CALLS_DIR = ROOT / "data" / "earlywarning" / "fks1_calls"
RUNBOOK = "work/RUNBOOK_fks1_fill.md"
# One stable label so a still-due month updates the open issue instead of spawning duplicates
# (the workflow keys its find-or-create on this).
ISSUE_LABEL = "fks1-fill-due"


def _load_records(args):
    """Resolve a FRESH record source: a written snapshot, an offline TSV, or a live NCBI pull.

    Records arrive in their raw `pending:` state (NCBI carries no AMR call); fill_status.fill_due
    overlays the committed call table itself, so this does not."""
    if args.snapshot:
        return ew.read_snapshot(args.snapshot), f"snapshot {args.snapshot}"
    if args.tsv:
        with open(args.tsv, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        records, _ = ew.normalize(rows)
        return records, f"offline tsv {args.tsv}"
    rows, tag, url = ew.fetch_release_rows(args.release)
    records, _ = ew.normalize(rows)
    return records, f"live pull {tag} ({url})"


def _issue_body(status, source, as_of, window_days):
    """The GitHub-issue body a workflow files when a refill is due -- the delta, how current the
    page is now, and the exact (non-billable) command to preview the workload before paying."""
    rel = status["newest_release"] or "(none yet)"
    through = status["committed_as_of"] or "—"
    cov = status["coverage_pct"]
    lines = [
        "## FKS1 fill is due — the resistance-weather page is drifting stale",
        "",
        f"**{status['reason']}**",
        "",
        f"- New in-window isolates uncalled (refill workload): **{status['uncalled_new']}**",
        f"- Window now: {status['resolved_in_window']}/{status['total_in_window']} isolates "
        f"re-called (**{cov:.0f}%** coverage) over a trailing {window_days}-day window "
        f"(as-of {as_of})",
        f"- Last committed fill: `{rel}`, calls current through **{through}**",
        f"- Threshold for this alert: {status['threshold']} new in-window isolate(s)",
        f"- Checked against: {source}",
        "",
        "### What to do",
        f"Run the billable GCP fill per [`{RUNBOOK}`]({RUNBOOK}). It is windowed, idempotent, and "
        "re-emits the union of calls-so-far, so it only pays for the new isolates. Preview the "
        "exact workload first (free, no cloud):",
        "",
        "```bash",
        "python scripts/snapshot_ncbi_auris.py pull        # free NCBI metadata pull",
        "python scripts/recall_fks1.py plan <snapshot> --json",
        "```",
        "",
        "Once the fill lands and the call table is committed, the weekly `weather.yml` publish "
        "re-lights the page automatically and this watchdog falls quiet.",
    ]
    return "\n".join(lines) + "\n"


def _write_github_output(path, status):
    """Emit workflow step-outputs (due=true/false, title, uncalled) if a path is given."""
    if not path:
        return
    title = (f"FKS1 fill due — {status['uncalled_new']} new in-window isolate(s) uncalled"
             if status["due"] else "")
    with open(path, "a") as fh:
        fh.write(f"due={'true' if status['due'] else 'false'}\n")
        fh.write(f"uncalled={status['uncalled_new']}\n")
        fh.write(f"label={ISSUE_LABEL}\n")
        if status["due"]:
            fh.write(f"title={title}\n")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--snapshot", help="read an existing written snapshot TSV")
    src.add_argument("--tsv", help="offline raw NCBI metadata TSV")
    ap.add_argument("--release", type=int, default=None,
                    help="pin an NCBI release to pull (default: latest)")
    ap.add_argument("--calls-dir", default=str(FKS1_CALLS_DIR),
                    help="committed FKS1 call-table dir to measure against")
    ap.add_argument("--window", type=int, default=180,
                    help="recent-window days (default 180, matches the weather page)")
    ap.add_argument("--as-of", default=None,
                    help="baseline cut date, ISO (default: today - window)")
    ap.add_argument("--threshold", type=int, default=fill_status.DEFAULT_THRESHOLD,
                    help=f"new in-window isolates that make a refill due "
                         f"(default {fill_status.DEFAULT_THRESHOLD})")
    ap.add_argument("--issue-out", default=None,
                    help="write the GitHub-issue body here when due (workflow files it)")
    ap.add_argument("--github-output", default=None,
                    help="append due/title/uncalled step-outputs here ($GITHUB_OUTPUT)")
    args = ap.parse_args()

    as_of = args.as_of or emergence.default_as_of(args.window)
    records, source = _load_records(args)
    status = fill_status.fill_due(records, args.calls_dir, as_of,
                                  window_days=args.window, threshold=args.threshold)

    verdict = "DUE" if status["due"] else "quiet"
    print(f"{verdict}: {status['reason']} (checked {source})")
    print(f"  in-window {status['resolved_in_window']}/{status['total_in_window']} re-called "
          f"({status['coverage_pct']:.0f}%); {status['uncalled_new']} new uncalled; "
          f"last fill {status['newest_release']} through {status['committed_as_of']}")

    if status["due"] and args.issue_out:
        pathlib.Path(args.issue_out).parent.mkdir(parents=True, exist_ok=True)
        pathlib.Path(args.issue_out).write_text(_issue_body(status, source, as_of, args.window))
        print(f"  wrote {args.issue_out} (issue body)")
    _write_github_output(args.github_output, status)
    return 0


if __name__ == "__main__":
    sys.exit(main())
