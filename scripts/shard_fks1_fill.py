#!/usr/bin/env python3
"""Shard a windowed FKS1 fill across N hosts, then merge the partial call tables (issue #171).

`recall_fks1.py fill` has no accession filter and writes ONE release-named call table from its
snapshot's resolved rows, so N hosts filling the same snapshot would clobber each other. This
splits the fill's bounded workload into N disjoint sub-snapshots (stable hash of run_acc, so the
assignment is reproducible and a re-split lands identically), each of which a host fills
independently with `recall_fks1.py fill <sub-snapshot> --no-window`; `merge` then unions the N
partial call tables by isolate_key into the final `PDG…<release>.tsv` + its INDEX provenance row.

The split/merge are pure TSV manipulation (no tools, no network) and run anywhere, including this
osx-arm64 host; only the per-shard `fill` is the billable linux-64 GCP step.

    # locally (non-billable):
    python scripts/shard_fks1_fill.py split <snapshot> --shards 8 --out-dir shards/
    # on host i (billable): run its one sub-snapshot
    python scripts/recall_fks1.py fill shards/PDG000000067.710.shard3-of-8.tsv --no-window
    # back locally after collecting the N partial call tables into a dir (non-billable):
    python scripts/shard_fks1_fill.py merge --calls-dir collected/ \
        --release PDG000000067.710 --shards 8 --as-of 2026-04-14 --window-days 180
"""
import argparse
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from openafr import earlywarning as ew  # noqa: E402

# Reuse the fill's exact window + eligibility definitions so a shard split covers precisely the
# rows `recall_fks1.py fill` would attempt -- no divergence between the two tools.
import importlib.util  # noqa: E402

_rf_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "recall_fks1.py")
_spec = importlib.util.spec_from_file_location("recall_fks1", _rf_path)
recall_fks1 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(recall_fks1)


def _shard_of(run_acc, n):
    """Stable shard index for a run accession: sha1 mod N. Deterministic across hosts/runs, so a
    re-split reproduces the same partition and shards stay disjoint."""
    first = run_acc.split(",")[0].strip()
    h = hashlib.sha1(first.encode("utf-8")).hexdigest()
    return int(h, 16) % n


def _shard_name(release_tag, i, n):
    return f"{release_tag}.shard{i}-of-{n}.tsv"


def cmd_split(args):
    records = ew.read_snapshot(args.snapshot)
    as_of, window_days = recall_fks1._resolve_window(args)
    window_keys = recall_fks1._in_window(records, as_of, window_days)
    pending = recall_fks1._pending_with_run(records, window_keys)
    release_tag = recall_fks1._release_tag_from_snapshot(args.snapshot)

    buckets = [[] for _ in range(args.shards)]
    for rec in pending:
        buckets[_shard_of(rec["run_acc"], args.shards)].append(rec)

    os.makedirs(args.out_dir, exist_ok=True)
    scope = (f"recent window as-of {as_of}, {window_days}d"
             if window_keys is not None else "whole snapshot (--no-window)")
    print(f"snapshot:    {args.snapshot}", file=sys.stderr)
    print(f"scope:       {scope} -> {len(pending)} pending-with-reads isolate(s)",
          file=sys.stderr)
    print(f"shards:      {args.shards}", file=sys.stderr)
    for i, bucket in enumerate(buckets):
        path = os.path.join(args.out_dir, _shard_name(release_tag, i, args.shards))
        ew.write_snapshot(path, bucket)
        print(f"  shard {i}: {len(bucket):5d} isolate(s) -> {path}", file=sys.stderr)
    sizes = [len(b) for b in buckets]
    print(f"balance:     min={min(sizes)} max={max(sizes)} "
          f"(each host runs: recall_fks1.py fill <shard> --no-window)", file=sys.stderr)
    return 0


def cmd_merge(args):
    merged = {}
    found = 0
    for i in range(args.shards):
        path = os.path.join(args.calls_dir, _shard_name(args.release, i, args.shards))
        if not os.path.exists(path):
            print(f"  shard {i}: MISSING {path}", file=sys.stderr)
            continue
        table = ew.read_fks1_call_table(path)
        found += 1
        overlap = merged.keys() & table.keys()
        if overlap:
            print(f"  WARNING shard {i}: {len(overlap)} isolate_key(s) already seen "
                  f"(shards should be disjoint): {sorted(overlap)[:3]}...", file=sys.stderr)
        merged.update(table)
        print(f"  shard {i}: {len(table):5d} call(s) <- {path}", file=sys.stderr)

    if found != args.shards:
        print(f"ERROR: only {found}/{args.shards} shard tables present; refusing a partial "
              f"merge (re-run after every host's call table is collected).", file=sys.stderr)
        return 1

    out_dir = args.out_dir or str(ew_default_calls_dir())
    os.makedirs(out_dir, exist_ok=True)
    final_path = os.path.join(out_dir, f"{args.release}.tsv")
    # merged values are already compact call-table rows (FKS1_CALLS_COLUMNS); serialize_fks1_calls
    # keeps only resolved rows and sorts by isolate_key, so the output is byte-stable.
    sha, n_called = ew.write_fks1_call_table(final_path, list(merged.values()))
    ew.update_fks1_calls_index(os.path.join(out_dir, "INDEX.tsv"), {
        "release_tag": args.release,
        "as_of": args.as_of,
        "window_days": args.window_days,
        "n_called": n_called,
        "sha256": sha,
        "written_at_utc": ew.utcnow_iso(),
    })
    print(f"merged {found} shard(s) -> {final_path} ({sha[:12]}..., {n_called} resolved call(s))",
          file=sys.stderr)
    print(f"INDEX row: release={args.release} as_of={args.as_of} window={args.window_days}",
          file=sys.stderr)
    return 0


def ew_default_calls_dir():
    return recall_fks1.FKS1_CALLS_DIR


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("split", help="split a snapshot's windowed fill workload into N "
                                     "disjoint sub-snapshots")
    s.add_argument("snapshot")
    s.add_argument("--shards", type=int, required=True)
    s.add_argument("--out-dir", required=True)
    recall_fks1._add_window_args(s)
    s.set_defaults(func=cmd_split)

    m = sub.add_parser("merge", help="union the N partial call tables into the final release "
                                     "call table + INDEX row")
    m.add_argument("--calls-dir", required=True,
                   help="dir holding the collected PDG....shardI-of-N.tsv partial call tables")
    m.add_argument("--release", required=True, help="release tag, e.g. PDG000000067.710")
    m.add_argument("--shards", type=int, required=True)
    m.add_argument("--as-of", required=True, help="window as-of for the INDEX provenance row")
    m.add_argument("--window-days", type=int, default=recall_fks1.DEFAULT_WINDOW_DAYS)
    m.add_argument("--out-dir", default=None,
                   help="where to write the final table (default: data/earlywarning/fks1_calls)")
    m.set_defaults(func=cmd_merge)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
