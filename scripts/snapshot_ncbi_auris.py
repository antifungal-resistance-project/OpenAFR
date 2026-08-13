"""Snapshot store CLI for the C. auris azole early-warning track (issue #21).

Persists each NCBI Pathogen Detection pull as a normalized, byte-stable snapshot
so emergence detection (#22) has a well-defined "what we saw before", and exposes
the two derived views the track needs: a baseline as of a date, and a diff between
two pulls. All heavy lifting lives in openafr.earlywarning; this is the CLI.

Usage
-----
  # Pull the latest release, write a snapshot, update the git-tracked INDEX:
  python scripts/snapshot_ncbi_auris.py pull

  python scripts/snapshot_ncbi_auris.py pull --release 671   # pin a release
  python scripts/snapshot_ncbi_auris.py pull --tsv saved.tsv --release 671  # offline
  python scripts/snapshot_ncbi_auris.py pull --dry-run       # fetch + report, no write

  # The baseline: isolates NCBI had as of a date.
  python scripts/snapshot_ncbi_auris.py baseline data/earlywarning/snapshots/PDG000000067.671.tsv --as-of 2025-01-01

  # What isolates are new between two snapshots (the emergence-detection input).
  python scripts/snapshot_ncbi_auris.py diff old.tsv new.tsv

Snapshot files are multi-MB and gitignored; the INDEX manifest (small) is
committed and is the durable, diffable memory of what has been pulled.
"""
import argparse
import collections
import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import earlywarning as ew  # noqa: E402

SNAP_DIR = pathlib.Path(__file__).resolve().parent.parent / "data" / "earlywarning" / "snapshots"
INDEX = SNAP_DIR / "INDEX.tsv"


def _load_local_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def cmd_pull(args):
    if args.tsv:
        rows = _load_local_tsv(args.tsv)
        if args.release is None:
            sys.exit("--release is required with --tsv (to name the snapshot)")
        tag = f"{ew.ORG_PREFIX}.{args.release}"
        url = "local:" + args.tsv
    else:
        rows, tag, url = ew.fetch_release_rows(args.release)

    records, dropped = ew.normalize(rows)
    n = len(records)
    sha = ew.snapshot_sha256(records)

    print(f"release {tag}: {len(rows)} raw rows -> {n} isolates "
          f"({dropped} unkeyable dropped)")
    print(f"snapshot sha256: {sha}")

    run_ok = sum(1 for r in records if r["run_acc"])
    dated = sum(1 for r in records if len(r["target_creation_date"][:10]) == 10)
    print(f"  with SRA reads (ERG11-callable, #23): {run_ok}/{n} = {100.0*run_ok/n:.1f}%")
    print(f"  with a creation date (baseline-placeable): {dated}/{n} = {100.0*dated/n:.1f}%")
    top = collections.Counter(r["country"] or "NULL" for r in records).most_common(5)
    print("  top countries: " + ", ".join(f"{c}={k}" for c, k in top))

    if args.dry_run:
        print("(dry-run: nothing written)")
        return

    SNAP_DIR.mkdir(parents=True, exist_ok=True)
    out = SNAP_DIR / f"{tag}.tsv"
    written = ew.write_snapshot(out, records)
    assert written == sha
    ew.update_index(INDEX, {
        "release_tag": tag,
        "pulled_at_utc": ew.utcnow_iso(),
        "n_isolates": str(n),
        "n_dropped": str(dropped),
        "sha256": sha,
        "source_url": url,
        "tool_version": ew.TOOL_VERSION,
    })
    print(f"wrote {out}")
    print(f"updated {INDEX}")


def cmd_baseline(args):
    records = ew.read_snapshot(args.snapshot)
    in_base, undated = ew.baseline_as_of(records, args.as_of)
    print(f"baseline as of {args.as_of}: {len(in_base)} isolates "
          f"(of {len(records)} in {pathlib.Path(args.snapshot).name})")
    print(f"  undated (cannot be placed in time, excluded): {len(undated)}")
    countries = collections.Counter(r["country"] or "NULL" for r in in_base)
    print(f"  countries represented: {len(countries)}")


def cmd_diff(args):
    old = ew.read_snapshot(args.old)
    new = ew.read_snapshot(args.new)
    d = ew.diff_snapshots(old, new)
    print(f"diff {pathlib.Path(args.old).name} -> {pathlib.Path(args.new).name}")
    print(f"  new isolates:     {len(d['added'])}")
    print(f"  removed isolates: {len(d['removed'])}")
    print(f"  revised (version bumped): {len(d['revised'])}")
    for rec in d["added"][:args.show]:
        print(f"    + {rec['isolate_key']}  {rec['country'] or '?':20s} "
              f"created {rec['target_creation_date'] or '?'}  run {rec['run_acc'] or '-'}")
    if len(d["added"]) > args.show:
        print(f"    ... and {len(d['added']) - args.show} more")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("pull", help="fetch a release, write a snapshot, update INDEX")
    p.add_argument("--release", type=int, help="pin a PDG release number (default: latest)")
    p.add_argument("--tsv", help="read a locally saved metadata.tsv instead of fetching")
    p.add_argument("--dry-run", action="store_true", help="report but do not write")
    p.set_defaults(func=cmd_pull)

    p = sub.add_parser("baseline", help="isolates seen as of a date")
    p.add_argument("snapshot")
    p.add_argument("--as-of", required=True, help="ISO date, e.g. 2025-01-01")
    p.set_defaults(func=cmd_baseline)

    p = sub.add_parser("diff", help="new/removed isolates between two snapshots")
    p.add_argument("old")
    p.add_argument("new")
    p.add_argument("--show", type=int, default=10, help="how many new isolates to list")
    p.set_defaults(func=cmd_diff)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
