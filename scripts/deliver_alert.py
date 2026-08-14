"""Delivery CLI for the C. auris azole early-warning track (issue #26).

Closes the loop end to end: pull -> detect -> map -> estimate -> compose -> DELIVER,
minimally and on a cadence. The scheduled runner (a GitHub Action, see
.github/workflows/earlywarning.yml) invokes `run`; a `demo` proves the whole loop --
including the quiet no-op and an escalation -- against the synthetic #22 fixture with
no network.

What "deliver" means here
-------------------------
Two artefacts, produced ONLY when there is genuine new news (see openafr.delivery):

  * a digest committed to the repo (data/earlywarning/digest/DIGEST.md), plus the
    small STATE.json that remembers what we have already alerted on, and
  * the body of a GitHub issue (written to --issue-out) that the workflow files.

A no-new-emergence run writes nothing, files nothing, and says so on one line. That
quiet is the design, not a degenerate case: an early-warning that pages weekly for the
same endemic mutation is an early-warning nobody reads.

The honest today-state
-----------------------
On real NCBI snapshots every erg11_call is still empty (calls are pending the ERG11
re-caller, #23), so the composer returns "no calls to assess" and `run` correctly
no-ops -- the schedule stays silent until there is a real call to stand behind, and
lights up automatically when #23 lands. `demo` runs the synthetic *called* snapshot so
delivery, escalation, and the no-op are all reproducible right now.

Usage
-----
  # Pull the latest release and deliver if there's new news:
  python scripts/deliver_alert.py run --as-of 2025-01-01

  python scripts/deliver_alert.py run --release 671 --as-of 2025-01-01   # pin a release
  python scripts/deliver_alert.py run --snapshot data/.../PDG000000067.671.tsv --as-of ...
  python scripts/deliver_alert.py run --tsv saved.tsv --as-of ...        # fully offline
  python scripts/deliver_alert.py run --as-of ... --dry-run              # decide, write nothing

  # Prove the loop (fixture -> deliver -> re-run no-op -> escalation):
  python scripts/deliver_alert.py demo
"""
import argparse
import csv
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from openafr import alert, delivery  # noqa: E402
from openafr import earlywarning as ew  # noqa: E402

DIGEST_DIR = ROOT / "data" / "earlywarning" / "digest"
DEFAULT_DIGEST = DIGEST_DIR / "DIGEST.md"
DEFAULT_STATE = DIGEST_DIR / "STATE.json"


def _load_records(args):
    """Resolve the run's record source: written snapshot, offline TSV, or a live pull."""
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


def _read_text(path):
    p = pathlib.Path(path)
    return p.read_text() if p.exists() else ""


def _write_github_output(path, delivered, title):
    """Emit workflow step outputs (delivered=true/false, title=...) if $GITHUB_OUTPUT set."""
    if not path:
        return
    with open(path, "a") as fh:
        fh.write(f"delivered={'true' if delivered else 'false'}\n")
        if delivered:
            fh.write(f"title={title}\n")


def _deliver(result, state_text, digest_path, state_path, issue_out,
             github_output, as_of, dry_run, source):
    """The shared decide-then-emit core used by both `run` and `demo`.

    Returns True when it delivered (new news), False on a quiet no-op.
    """
    state = delivery.load_state(state_text)
    fresh = delivery.new_alerts(result, state)

    if not fresh:
        why = result.get("note", "no new or escalated alerts over the thresholds")
        print(f"no-op: {why} — nothing delivered ({source}).")
        _write_github_output(github_output, False, "")
        return False

    title = delivery.issue_title(result, fresh)
    digest = delivery.render_digest(result, fresh, state)
    body = delivery.issue_body(result, fresh, state)
    new_state = delivery.next_state(result, state, as_of=as_of)

    print(f"DELIVER: {len(fresh)} new/escalated alert(s) — {title}")
    for a in fresh:
        print(f"  - {a['mutation']} [{alert._BADGE[a['priority']]}]")

    if dry_run:
        print("(dry-run: nothing written)")
        _write_github_output(github_output, False, "")
        return True

    pathlib.Path(digest_path).parent.mkdir(parents=True, exist_ok=True)
    pathlib.Path(digest_path).write_text(digest)
    pathlib.Path(state_path).write_text(delivery.dumps_state(new_state))
    print(f"wrote {digest_path}")
    print(f"wrote {state_path}")
    if issue_out:
        pathlib.Path(issue_out).parent.mkdir(parents=True, exist_ok=True)
        pathlib.Path(issue_out).write_text(body)
        print(f"wrote {issue_out} (issue body)")
    _write_github_output(github_output, True, title)
    return True


def cmd_run(args):
    records, source = _load_records(args)
    result = alert.compose_alerts(
        records, args.as_of, window_days=args.window, min_count=args.min_count,
        min_delta=args.min_delta, backlog_frac=args.backlog_frac)
    _deliver(result, _read_text(args.state), args.digest, args.state, args.issue_out,
             args.github_output, args.as_of, args.dry_run, source)


def cmd_demo(args):
    """End-to-end proof on the synthetic fixture: deliver, re-run no-op, escalation."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from detect_emergence import _demo_records  # the shared #22 fixture

    records = _demo_records()
    result = alert.compose_alerts(records, "2024-12-31", window_days=180,
                                  min_count=3, min_delta=0.05)

    print("=== first run (empty state): everything is new ===")
    state_text = ""  # no prior deliveries
    delivered = delivery.new_alerts(result, delivery.load_state(state_text))
    for a in delivered:
        print(f"  DELIVER {a['mutation']} [{alert._BADGE[a['priority']]}]")
    state_text = delivery.dumps_state(
        delivery.next_state(result, delivery.load_state(state_text)))

    print("\n=== second run (same snapshot): quiet no-op ===")
    if delivery.is_noop(result, delivery.load_state(state_text)):
        print("  no-op: nothing new — nothing delivered (as designed).")
    else:
        print("  ERROR: expected a no-op on an unchanged re-run")

    print("\n=== third run: TR34 backlog clears, context -> watch (escalation) ===")
    escalated = _clear_tr34_backlog(records)
    result2 = alert.compose_alerts(escalated, "2024-12-31", window_days=180,
                                   min_count=3, min_delta=0.05)
    fresh = delivery.new_alerts(result2, delivery.load_state(state_text))
    for a in fresh:
        print(f"  DELIVER {a['mutation']} [{alert._BADGE[a['priority']]}] (escalation)")
    if not fresh:
        print("  (no escalation surfaced — fixture may not force one)")

    if args.markdown and fresh:
        print("\n--- digest that would be committed ---\n")
        print(delivery.render_digest(result2, fresh, delivery.load_state(state_text)))


def _clear_tr34_backlog(records):
    """Copy the fixture with TR34's carriers freshly collected instead of archival.

    In the base fixture TR34 spikes but every carrier was *collected* in 2019, so the
    backlog guard files it as CONTEXT. Here we move those collection dates into the
    recent window: the spike is now genuine, so the composer re-tiers TR34 up out of
    CONTEXT -- a real escalation the delivery rule must re-fire on, driven entirely by
    the real detector, not a hand-forged verdict.
    """
    out = []
    for r in records:
        r = dict(r)
        if "TR34" in r["erg11_call"]:
            r["collection_date"] = "2024"
        out.append(r)
    return out


def _add_source_args(p):
    src = p.add_mutually_exclusive_group()
    src.add_argument("--snapshot", help="read an existing written snapshot TSV")
    src.add_argument("--tsv", help="offline raw NCBI metadata TSV")
    p.add_argument("--release", type=int, default=None,
                   help="pin an NCBI release to pull (default: latest)")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("run", help="pull/compose and deliver if there's new news")
    _add_source_args(p)
    p.add_argument("--as-of", required=True, help="baseline cut date, ISO")
    p.add_argument("--window", type=int, default=180, help="recent-window days (default 180)")
    p.add_argument("--min-count", type=int, default=3, help="min recent carriers to flag")
    p.add_argument("--min-delta", type=float, default=0.05, help="min share increase for RISING")
    p.add_argument("--backlog-frac", type=float, default=0.5,
                   help="archival-carrier fraction that trips the backlog guard")
    p.add_argument("--digest", default=str(DEFAULT_DIGEST), help="committed digest path")
    p.add_argument("--state", default=str(DEFAULT_STATE), help="committed delivery-state path")
    p.add_argument("--issue-out", default=None,
                   help="write the GitHub-issue body here (workflow files it)")
    p.add_argument("--github-output", default=None,
                   help="append delivered/title step-outputs here ($GITHUB_OUTPUT)")
    p.add_argument("--dry-run", action="store_true", help="decide + report, write nothing")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("demo", help="prove the loop on the synthetic fixture")
    p.add_argument("--markdown", action="store_true", help="print the digest that would ship")
    p.set_defaults(func=cmd_demo)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
