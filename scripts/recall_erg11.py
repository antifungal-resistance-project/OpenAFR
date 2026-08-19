"""ERG11 re-caller CLI -- the reads/consensus -> azole-resistance call (TODOS top-priority).

Fills the `erg11_call` column the snapshot store leaves `pending:sra-recaller`, so the
early-warning pipeline (emergence -> mapping -> structural so-what -> alert) can fire on
real isolates instead of an honest null. The deterministic call logic lives in
openafr.recaller (unit-tested against the pinned reference); this is the CLI, including
the tool-and-network orchestration that turns an isolate's SRA reads into the consensus
CDS the core consumes.

Two boundaries, deliberately separated
--------------------------------------
`call`  -- deterministic, no network, no external tools. Given a consensus ERG11 CDS
           FASTA (in-frame, same length as the reference), emit the call. This is the
           reproducible path the tests exercise and the one to trust.

`recall`/`fill` -- the messy half: SRA reads -> aligned -> consensus CDS -> `call`.
           Shells out to `fasterq-dump` (sra-tools), `minimap2`, and `samtools`
           (>=1.13, for `samtools consensus`). GATED on those binaries the same way the
           docking scripts gate on `vina`/`obabel`: if they are absent this exits with
           an install message rather than pretending. This half needs network + multi-GB
           downloads, so it is intentionally NOT in the tested core -- it is honestly
           labelled as the un-reproducible-here layer.

Usage
-----
  # Deterministic: call from an already-built consensus CDS (no network).
  python scripts/recall_erg11.py call --consensus-fasta isolate_erg11.fasta

  # Orchestrated: one isolate's SRA run -> call (needs the external tools + network).
  python scripts/recall_erg11.py recall SRR1234567

  # Offline preflight (no tools/network): what would a fill download and re-call?
  python scripts/recall_erg11.py plan data/earlywarning/snapshots/PDG000000067.671.tsv

  # Batch: fill every pending:sra-recaller row of a snapshot in place.
  python scripts/recall_erg11.py fill data/earlywarning/snapshots/PDG000000067.671.tsv

Honesty
-------
An isolate whose panel residues could not be covered is written
`resistance_source = sra-recaller:partial(...)` with an EMPTY call -- never a silent
"wild type". An isolate that fails to fetch/align is left `pending:sra-recaller` (or
marked `sra-recaller:failed(<reason>)`), not dropped and not guessed.
"""
import argparse
import os
import pathlib
import random
import re
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import recaller as R          # noqa: E402
from openafr import earlywarning as ew      # noqa: E402
from openafr import runlog                  # noqa: E402  (append-only per-run provenance)

# External binaries the orchestration layer needs. Kept in one place so the gate
# message and the actual calls cannot drift.
TOOLS = {
    "fasterq-dump": "sra-tools (SRA read download)",
    "minimap2": "minimap2 (short-read alignment)",
    "samtools": "samtools >=1.13 (samtools consensus)",
}
MIN_DEPTH = 10   # positions with fewer supporting reads become N -> uncalled, not wild type

# The consensus pipeline's LAST step, `samtools consensus`, only exists in samtools
# >=1.13. Presence on PATH is not enough: an older samtools passes `shutil.which` and then
# dies on `unrecognized command: consensus` -- but only AFTER fasterq-dump has pulled the
# isolate's multi-GB reads. We refuse up front so that failure is instant and actionable.
SAMTOOLS_MIN = (1, 13)


def _read_consensus_fasta(path):
    seqs = []
    cur = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    if not seqs:
        sys.exit(f"no sequence in {path}")
    if len(seqs) > 1:
        sys.exit(f"expected one sequence in {path}, found {len(seqs)}")
    return seqs[0]


# --------------------------- deterministic path -----------------------------

def cmd_call(args):
    """Consensus CDS FASTA -> (erg11_call, resistance_source). No network/tools."""
    reference = R.load_reference()
    consensus = _read_consensus_fasta(args.consensus_fasta)
    try:
        result = R.call_substitutions(consensus, reference=reference)
    except R.ConsensusError as e:
        # A frameshift/length mismatch is a real, honest outcome -- report it as the
        # source, do not crash the caller into pretending nothing was sequenced.
        print(f"\t{args.consensus_fasta}\tsra-recaller:refused({e})", file=sys.stderr)
        print(f"\tsra-recaller:refused")
        return 2
    call, source = R.format_call(result)
    print(f"call:   {call or '(none -- no substitution)'}")
    print(f"source: {source}")
    print(f"panel:  {', '.join(result['panel_hits']) or '(none)'}")
    print(f"covered: {result['n_covered']}/{result['n_residues']} residues"
          + (f"; uncalled panel positions: {result['uncalled_panel']}"
             if result["uncalled_panel"] else ""))
    return 0


# --------------------------- orchestration path -----------------------------

def _require_tools():
    missing = [(t, why) for t, why in TOOLS.items() if shutil.which(t) is None]
    if missing:
        lines = "\n".join(f"  - {t}: {why}" for t, why in missing)
        sys.exit(
            "the reads->consensus orchestration needs external tools that are not on "
            f"PATH:\n{lines}\n"
            "Install them (conda-forge/bioconda), or run the deterministic `call` "
            "subcommand on a consensus FASTA you built elsewhere."
        )
    _require_samtools_consensus()


def _run_text(cmd):
    """Capture a probe command's combined output; '' if it cannot be run. Never raises."""
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return ""
    return (p.stdout or "") + (p.stderr or "")


def _parse_version(text):
    """First 'MAJOR.MINOR' in text -> (major, minor); None if none present."""
    m = re.search(r"(\d+)\.(\d+)", text)
    return (int(m.group(1)), int(m.group(2))) if m else None


def samtools_capabilities():
    """(version_tuple_or_None, has_consensus) for the samtools on PATH.

    `consensus` is listed by `samtools --help` exactly on the versions that ship it, so
    that listing -- not the parsed version -- is the ground-truth capability check; the
    version is parsed only to make the error message concrete.
    """
    ver = _parse_version(_run_text(["samtools", "--version"]))
    has_consensus = re.search(r"\bconsensus\b", _run_text(["samtools", "--help"])) is not None
    return ver, has_consensus


def _require_samtools_consensus():
    ver, has_consensus = samtools_capabilities()
    if has_consensus and (ver is None or ver >= SAMTOOLS_MIN):
        return
    vstr = ".".join(map(str, ver)) if ver else "unknown"
    need = ".".join(map(str, SAMTOOLS_MIN))
    sys.exit(
        f"samtools {vstr} lacks the `consensus` subcommand (needs >= {need}). The "
        "reads->consensus orchestration would fail at its FINAL step -- after the "
        "multi-GB read download -- so it is refused now instead. Upgrade via "
        f"conda-forge/bioconda (`conda install -c bioconda 'samtools>={need}'`) and retry."
    )


def _run(cmd, **kw):
    return subprocess.run(cmd, check=True, **kw)


def reads_to_consensus(run_acc, workdir, reference_fasta, min_depth=MIN_DEPTH):
    """SRA run accession -> ERG11 consensus CDS FASTA path. Shells out; not tested here.

    Aligns the isolate's reads to the ERG11 CDS reference, then emits a consensus in
    reference coordinates with N wherever depth < min_depth -- so uncovered hotspots
    arrive at the core as uncalled, exactly the honesty the core relies on.
    """
    wd = pathlib.Path(workdir)
    # 1) fetch reads
    _run(["fasterq-dump", "--split-spot", "--skip-technical", "-O", str(wd), run_acc])
    fastqs = sorted(str(p) for p in wd.glob(f"{run_acc}*.fastq"))
    if not fastqs:
        raise RuntimeError(f"fasterq-dump produced no FASTQ for {run_acc}")
    # 2) align to the ERG11 CDS (short-read preset), sort, index
    sam = wd / "aln.sam"
    with open(sam, "w") as fh:
        _run(["minimap2", "-ax", "sr", reference_fasta, *fastqs], stdout=fh)
    bam = wd / "aln.sorted.bam"
    _run(["samtools", "sort", "-o", str(bam), str(sam)])
    _run(["samtools", "index", str(bam)])
    # 3) consensus in reference coordinates; low coverage -> N (uncalled)
    cons = wd / "erg11_consensus.fasta"
    with open(cons, "w") as fh:
        _run(["samtools", "consensus", "-f", "fasta", "--min-depth", str(min_depth),
              "--call-fract", "0.5", str(bam)], stdout=fh)
    return str(cons)


def _recall_one(run_acc, reference, reference_fasta, keep_tmp=False):
    """Orchestrate one isolate: run_acc -> (call, source). Returns ('', failed-source)
    on any tool/fetch error rather than raising -- a failed isolate is honestly marked,
    never dropped."""
    tmp = tempfile.mkdtemp(prefix=f"recall_{run_acc}_")
    try:
        cons_path = reads_to_consensus(run_acc, tmp, reference_fasta)
        consensus = _read_consensus_fasta(cons_path)
        result = R.call_substitutions(consensus, reference=reference)
        return R.format_call(result)
    except R.ConsensusError as e:
        return "", f"sra-recaller:refused({_short(e)})"
    except (subprocess.CalledProcessError, RuntimeError, OSError) as e:
        return "", f"sra-recaller:failed({_short(e)})"
    finally:
        if not keep_tmp:
            shutil.rmtree(tmp, ignore_errors=True)


def _short(e):
    return str(e).splitlines()[0][:80].replace(",", ";")


def cmd_recall(args):
    _require_tools()
    reference = R.load_reference()
    with runlog.record_run("recall", reference_path=R.DEFAULT_REFERENCE,
                           tools=list(TOOLS),
                           extra={"run_acc": args.run_acc}) as run:
        call, source = _recall_one(args.run_acc, reference, R.DEFAULT_REFERENCE,
                                   keep_tmp=args.keep_tmp)
        ok = source.startswith("sra-recaller:") and "failed" not in source
        run.add_item(run_acc=args.run_acc, call=call, source=source)
        run.set(status=("ok" if ok else "fail"), call=call, source=source)
    print(f"run:    {args.run_acc}")
    print(f"call:   {call or '(none)'}")
    print(f"source: {source}")
    return 0 if ok else 1


def cmd_fill(args):
    """Fill every pending:sra-recaller row of a snapshot in place (orchestration)."""
    _require_tools()
    reference = R.load_reference()
    records = ew.read_snapshot(args.snapshot)
    pending = _fill_order(records, args.limit, args.random, args.seed)
    how = f"random sample, seed={args.seed}" if args.random else "snapshot order"
    print(f"{len(pending)} isolate(s) to re-call (pending:sra-recaller with a run_acc; {how})",
          file=sys.stderr)
    # A `fill` overwrites the snapshot's calls in place, so without a log a re-run erases
    # its predecessor with no trace. Log every fill (append-only) with the before/after
    # snapshot digest, the per-isolate transcript, and full tool/code provenance. Written
    # on exit even if a mid-batch isolate crashes -- the digests then show it never wrote.
    digest_before = ew.snapshot_sha256(records)
    n_called = n_partial = n_failed = 0
    with runlog.record_run("fill", reference_path=R.DEFAULT_REFERENCE,
                           tools=list(TOOLS),
                           extra={"snapshot": os.path.basename(args.snapshot),
                                  "snapshot_sha256_before": digest_before,
                                  "n_pending": len(pending),
                                  "limit": args.limit or None,
                                  "random": bool(args.random),
                                  "seed": args.seed if args.random else None,
                                  "dry_run": bool(args.dry_run)}) as run:
        for i, rec in enumerate(pending, 1):
            run_acc = rec["run_acc"].split(",")[0].strip()   # first linked run
            call, source = _recall_one(run_acc, reference, R.DEFAULT_REFERENCE)
            rec["erg11_call"] = call
            rec["resistance_source"] = source
            if "failed" in source or "refused" in source:
                n_failed += 1
            elif "partial" in source:
                n_partial += 1
            elif call:
                n_called += 1
            run.add_item(isolate_key=rec["isolate_key"], run_acc=run_acc,
                         call=call, source=source)
            print(f"[{i}/{len(pending)}] {rec['isolate_key']} {run_acc} -> "
                  f"{call or '(no sub)'} [{source}]", file=sys.stderr)
        digest_after = ew.snapshot_sha256(records)
        if not args.dry_run:
            digest = ew.write_snapshot(args.snapshot, records)
            print(f"wrote {args.snapshot} ({digest[:12]}...)", file=sys.stderr)
        run.set(snapshot_sha256_after=digest_after, wrote=not args.dry_run,
                summary={"called": n_called, "partial": n_partial, "failed": n_failed,
                         "attempted": len(pending)})
    print(f"summary: {n_called} with substitution(s), {n_partial} partial, "
          f"{n_failed} failed/refused, of {len(pending)} attempted", file=sys.stderr)
    return 0


# ------------------------------ offline preflight ----------------------------

def _pending_with_run(records):
    """The rows `fill` would attempt, in fill's order: pending:sra-recaller + a run_acc."""
    return [r for r in records
            if (r.get("resistance_source") or "") == "pending:sra-recaller"
            and (r.get("run_acc") or "").strip()]


def _fill_order(records, limit=0, shuffle=False, seed=0):
    """The exact rows a `fill` will attempt, in the order it will attempt them.

    Shared by cmd_fill and cmd_plan so the offline preflight matches the real run.

    Snapshot order is not random -- it tracks NCBI accession/submission order, which
    correlates with outbreak/study, so the head of the list is a biased sample (an early
    `--limit 20` came back all one resistance study). With `shuffle`, deterministically
    permute (seeded) BEFORE slicing, so a `--limit N` is a random *representative* sample
    of the pending pool rather than its head. The seed keeps it reproducible: the same
    snapshot + seed always yield the same sample -- the same rebuild contract as the pinned
    environment and hash-frozen pre-registrations."""
    rows = _pending_with_run(records)
    if shuffle:
        random.Random(seed).shuffle(rows)
    return rows[:limit] if limit else rows


def cmd_plan(args):
    """Read-only preflight: what a `fill` WOULD fetch and re-call. No tools, no network.

    Runnable anywhere (incl. osx-arm64, where the orchestration itself cannot run), so a
    snapshot can be inspected -- pending count, run_acc coverage, the exact accessions the
    next `fill` will download -- before committing to the multi-GB downloads. Mirrors
    cmd_fill's row selection exactly (same source/run_acc filter, same first-linked-run
    pick, same --limit slice), so its counts ARE the next fill's workload."""
    records = ew.read_snapshot(args.snapshot)
    pending = [r for r in records
               if (r.get("resistance_source") or "") == "pending:sra-recaller"]
    with_run = _pending_with_run(records)
    without_run = len(pending) - len(with_run)
    attempted = _fill_order(records, args.limit, args.random, args.seed)
    first_runs = [r["run_acc"].split(",")[0].strip() for r in attempted]
    distinct = sorted(set(first_runs))

    print(f"snapshot:            {args.snapshot}")
    print(f"total isolates:      {len(records)}")
    print(f"pending:sra-recaller {len(pending)}"
          f"  ({without_run} with no run_acc -> skipped, stay pending)")
    print(f"fill would re-call:  {len(attempted)} isolate(s)"
          + (f"  (--limit {args.limit})" if args.limit else "")
          + (f"  [random sample, seed={args.seed}]" if args.random else "")
          + f"  == {len(attempted)} SRA download(s)")
    print(f"distinct run acc:    {len(distinct)}"
          + ("" if len(distinct) == len(attempted)
             else "  (some isolates share a run; each row still downloads once)"))
    if args.list_runs:
        for acc in distinct:
            print(acc)
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("call", help="deterministic: consensus CDS FASTA -> call (no network)")
    p.add_argument("--consensus-fasta", required=True,
                   help="in-frame ERG11 consensus CDS, same length as the reference")
    p.set_defaults(func=cmd_call)

    p = sub.add_parser("recall", help="orchestrated: one SRA run accession -> call")
    p.add_argument("run_acc", help="SRA Run accession (e.g. SRR1234567)")
    p.add_argument("--keep-tmp", action="store_true", help="keep the alignment workdir")
    p.set_defaults(func=cmd_recall)

    p = sub.add_parser("fill", help="orchestrated: fill a snapshot's pending rows in place")
    p.add_argument("snapshot", help="snapshot TSV to update in place")
    p.add_argument("--limit", type=int, default=0, help="only re-call N pending isolates")
    p.add_argument("--random", action="store_true",
                   help="sample the --limit slice randomly (seeded), not from snapshot order "
                        "-- for a representative frequency, since snapshot order is study-biased")
    p.add_argument("--seed", type=int, default=0,
                   help="RNG seed for --random (default 0); same snapshot+seed -> same sample")
    p.add_argument("--dry-run", action="store_true", help="re-call but do not write the snapshot")
    p.set_defaults(func=cmd_fill)

    p = sub.add_parser("plan", help="offline preflight: what a fill would fetch (no tools/network)")
    p.add_argument("snapshot", help="snapshot TSV to inspect (not modified)")
    p.add_argument("--limit", type=int, default=0, help="preview fill's N-isolate slice")
    p.add_argument("--random", action="store_true",
                   help="preview the seeded random sample fill --random would take")
    p.add_argument("--seed", type=int, default=0, help="RNG seed for --random (default 0)")
    p.add_argument("--list-runs", action="store_true", help="print the distinct SRA accessions")
    p.set_defaults(func=cmd_plan)

    args = ap.parse_args()
    return args.func(args) or 0


if __name__ == "__main__":
    sys.exit(main())
