"""FKS1 re-caller CLI -- the reads/consensus -> echinocandin-resistance call (v2, T2).

The FKS1 analog of scripts/recall_erg11.py. The azole/ERG11 track measured an ~80% azole
event frequency (work/RESULTS_prevalence.md), which saturates the azole *emergence* signal;
echinocandin (FKS1) resistance is the dynamic axis with headroom. The detection core lives in
openafr.fks1_caller (unit-tested against the pinned B8441 GSC1/FKS1 reference); this is the
CLI, including the tool-and-network orchestration that turns an isolate's SRA reads into the
per-window consensus the core consumes.

DETECTION ONLY. FKS1 is a 1,3-beta-glucan synthase; echinocandins do not coordinate a metal,
so the CYP51/heme-iron geometry moat does NOT transfer and no structural "so-what" is claimed
for an FKS1 call. This CLI emits substitution tokens for a downstream emergence/alert pipeline;
it never borrows the azole pocket verdict.

Two boundaries, deliberately separated (same split as the ERG11 re-caller)
--------------------------------------------------------------------------
`call`   -- deterministic, no network, no external tools. Given an isolate's consensus over
            the hot-spot windows (either a full-length CDS FASTA, or one FASTA per window),
            emit the call. This is the reproducible path the tests exercise.

`recall` -- the messy half: SRA reads -> aligned -> per-window consensus -> `call`. Shells out
            to `fasterq-dump` (sra-tools), `minimap2`, and `samtools` (>=1.13, for
            `samtools consensus`). GATED on those binaries exactly like recall_erg11.py: if
            they are absent this exits with an install message rather than pretending. Needs
            network + multi-GB downloads, so it is intentionally NOT in the tested core.

The one departure from the ERG11 orchestration: WINDOWED consensus
------------------------------------------------------------------
ERG11 builds ONE whole-CDS consensus and refuses anything whose length != reference. FKS1 is
~5.6 kb (1888 aa) and all echinocandin resistance lives in two short hot-spots (HS1 ~S639,
HS2 ~R1354). Building one 5.6 kb consensus makes the length-exact frame contract fragile --
a single low-coverage indel anywhere would discard an isolate whose hot-spots were perfectly
covered. So `recall` emits ONE consensus PER hot-spot region (`samtools consensus -r`), in
reference coordinates: an indel in one window's reads is refused for that window alone and
never corrupts the other. The core (fks1_caller.call_windows) then calls each window
independently.

Batch over an NCBI snapshot (T3) -- the ERG11 batch path, now for FKS1
---------------------------------------------------------------------
`plan`/`fill`/`prevalence`/`migrate` mirror scripts/recall_erg11.py, over the FKS1 columns
the snapshot schema now carries (openafr.earlywarning.SNAPSHOT_COLUMNS gained
`fks1_call`/`fks1_resistance_source` in T3; a pre-v2 snapshot is migrated on read, seeding
them to an honest pending state from run_acc):

  plan       -- offline preflight: what a fill WOULD fetch/re-call (no tools/network)
  fill       -- orchestrated: fill every `pending:sra-fks1-recaller` row in place
  prevalence -- offline: the echinocandin event frequency over RESOLVED isolates
                (undefined -- 'NOT MEASURED YET' -- until a fill populates fks1_call)
  migrate    -- rewrite a pre-v2 snapshot to the current schema in place

`plan`/`prevalence`/`migrate` run anywhere (incl. osx-arm64); only `fill` needs the external
tools + network, gated exactly like `recall`.

Honesty
-------
A window whose hot-spot residues could not be covered is reported per-window as
`sra-fks1-recaller:HS1=partial(uncalled:639),...` with those residues left OUT of the call --
never a silent "wild type". A window with an in-window indel is `refused(...)` for that window
alone. An isolate that fails to fetch/align yields an empty call with a `failed(<reason>)`
source, not a dropped or guessed one.
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
from openafr import backtest as bt           # noqa: E402  (fks1_panel_prevalence)
from openafr import earlywarning as ew        # noqa: E402  (snapshot store + migration)
from openafr import fks1_caller as F         # noqa: E402
from openafr import runlog                   # noqa: E402  (append-only per-run provenance)

# External binaries the orchestration layer needs. Same toolchain as the ERG11 re-caller;
# kept in one place so the gate message and the actual calls cannot drift.
TOOLS = {
    "fasterq-dump": "sra-tools (SRA read download)",
    "minimap2": "minimap2 (short-read alignment)",
    "samtools": "samtools >=1.13 (samtools consensus)",
}
MIN_DEPTH = 10   # positions with fewer supporting reads become N -> uncalled, not wild type

# `samtools consensus` only exists in samtools >=1.13. Presence on PATH is not enough: an
# older samtools passes `shutil.which` and then dies on `unrecognized command: consensus` --
# but only AFTER fasterq-dump has pulled the isolate's multi-GB reads. We refuse up front.
SAMTOOLS_MIN = (1, 13)


def _read_consensus_fasta(path):
    """Read a single-record consensus FASTA -> the bare sequence string. Refuses empty /
    header-only / multi-record input (the same bridge the ERG11 CLI uses), so an upstream
    tool hiccup can never read as a silently wild-type empty sequence."""
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

def _print_result(call, source, result):
    print(f"call:   {call or '(none -- no substitution)'}")
    print(f"source: {source}")
    print(f"panel:  {', '.join(result['panel_hits']) or '(none)'}")
    for name in F.FKS1_WINDOWS:
        w = result["windows"][name]
        cov = f"{w['n_covered']}/{w['n_residues']}"
        detail = ""
        if w["status"] == "refused":
            detail = f"  refused: {w['refused']}"
        elif w["uncalled_positions"]:
            detail = f"  uncalled: {'+'.join(map(str, w['uncalled_positions']))}"
        print(f"  {name}: {w['status']:9s} covered {cov}{detail}")


def cmd_call(args):
    """Consensus -> (fks1_call, resistance_source). No network/tools.

    Accepts EITHER a full-length CDS consensus (--consensus-fasta; sliced to windows, the
    fragile length-exact path) OR one consensus FASTA per window (--window NAME=PATH; the
    robust path `recall` produces). A window not supplied is reported `missing`, never
    assumed wild-type."""
    reference = F.load_reference()
    if args.window:
        window_seqs = {}
        for spec in args.window:
            name, sep, path = spec.partition("=")
            name = name.strip()
            if not sep or not path:
                sys.exit(f"--window expects NAME=PATH, got {spec!r}")
            if name not in F.FKS1_WINDOWS:
                sys.exit(f"unknown window {name!r}; known: {', '.join(F.FKS1_WINDOWS)}")
            window_seqs[name] = _read_consensus_fasta(path)
    else:
        full = _read_consensus_fasta(args.consensus_fasta)
        try:
            window_seqs = F.slice_windows_from_cds(full)
        except F.WindowError as e:
            # A whole-consensus that cannot be anchored is a real, honest outcome -- report
            # it, do not crash into pretending the windows were covered.
            print(f"source: sra-fks1-recaller:refused({e})", file=sys.stderr)
            print("refused")
            return 2
    result = F.call_windows(window_seqs, reference=reference)
    call, source = F.format_call(result)
    _print_result(call, source, result)
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

    `consensus` is listed by `samtools --help` exactly on the versions that ship it, so that
    listing -- not the parsed version -- is the ground-truth capability check; the version is
    parsed only to make the error message concrete."""
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


def _reference_contig(reference_fasta):
    """The reference sequence name minimap2 will emit into the BAM (the FASTA header's first
    token) -- what `samtools consensus -r` must address the hot-spot regions by."""
    with open(reference_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                return line[1:].split()[0]
    raise RuntimeError(f"no FASTA header in {reference_fasta}")


def _window_region(contig, window):
    """`samtools consensus -r` region string for a hot-spot window, in 1-based inclusive
    reference-CDS nucleotide coordinates. Derived from the window's residue span so it can
    never drift from the core's numbering (residue r -> nt (r-1)*3+1 .. r*3)."""
    start_nt = (window.start_res - 1) * 3 + 1
    end_nt = window.end_res * 3
    return f"{contig}:{start_nt}-{end_nt}"


def reads_to_window_consensus(run_acc, workdir, reference_fasta, min_depth=MIN_DEPTH):
    """SRA run accession -> {window_name: consensus_nt}. Shells out; not tested here.

    The windowed departure from ERG11: align the isolate's reads to the full FKS1 CDS once,
    then emit ONE consensus per hot-spot region in reference coordinates. `-a` fills every
    region position (N where depth < min_depth), so an uncovered hot-spot residue arrives at
    the core as uncalled -- exactly the honesty the core relies on -- and an indel in one
    window is confined to that window's consensus, never discarding the other."""
    wd = pathlib.Path(workdir)
    contig = _reference_contig(reference_fasta)
    # 1) fetch reads
    _run(["fasterq-dump", "--split-spot", "--skip-technical", "-O", str(wd), run_acc])
    fastqs = sorted(str(p) for p in wd.glob(f"{run_acc}*.fastq"))
    if not fastqs:
        raise RuntimeError(f"fasterq-dump produced no FASTQ for {run_acc}")
    # 2) align to the full FKS1 CDS (short-read preset), sort, index
    sam = wd / "aln.sam"
    with open(sam, "w") as fh:
        _run(["minimap2", "-ax", "sr", reference_fasta, *fastqs], stdout=fh)
    bam = wd / "aln.sorted.bam"
    _run(["samtools", "sort", "-o", str(bam), str(sam)])
    _run(["samtools", "index", str(bam)])
    # 3) one consensus per hot-spot window, in reference coordinates; low cov -> N (uncalled)
    window_seqs = {}
    for name, w in F.FKS1_WINDOWS.items():
        region = _window_region(contig, w)
        cons = wd / f"fks1_{name}_consensus.fasta"
        with open(cons, "w") as fh:
            _run(["samtools", "consensus", "-f", "fasta", "-a", "-r", region,
                  "--min-depth", str(min_depth), "--call-fract", "0.5", str(bam)],
                 stdout=fh)
        window_seqs[name] = _read_consensus_fasta(str(cons))
    return window_seqs


def _short(e):
    return str(e).splitlines()[0][:80].replace(",", ";")


def _recall_one(run_acc, reference, reference_fasta, keep_tmp=False):
    """Orchestrate one isolate: run_acc -> (call, source). Returns ('', failed-source) on any
    tool/fetch error rather than raising -- a failed isolate is honestly marked, never dropped.

    Note: unlike the ERG11 core, per-window frame refusals are NOT exceptions here; they are
    carried inside the call result (window status 'refused') and rendered into the source by
    format_call. Only tool/fetch/reference failures land in the failed(...) branch."""
    tmp = tempfile.mkdtemp(prefix=f"fks1_recall_{run_acc}_")
    try:
        window_seqs = reads_to_window_consensus(run_acc, tmp, reference_fasta)
        result = F.call_windows(window_seqs, reference=reference)
        return F.format_call(result)
    except (F.ReferenceError, subprocess.CalledProcessError, RuntimeError, OSError) as e:
        return "", f"sra-fks1-recaller:failed({_short(e)})"
    finally:
        if not keep_tmp:
            shutil.rmtree(tmp, ignore_errors=True)


def cmd_recall(args):
    _require_tools()
    reference = F.load_reference()
    with runlog.record_run("recall-fks1", reference_path=F.DEFAULT_REFERENCE,
                           tools=list(TOOLS),
                           extra={"run_acc": args.run_acc}) as run:
        call, source = _recall_one(args.run_acc, reference, F.DEFAULT_REFERENCE,
                                   keep_tmp=args.keep_tmp)
        ok = source.startswith("sra-fks1-recaller:") and "failed" not in source
        run.add_item(run_acc=args.run_acc, call=call, source=source)
        run.set(status=("ok" if ok else "fail"), call=call, source=source)
    print(f"run:    {args.run_acc}")
    print(f"call:   {call or '(none)'}")
    print(f"source: {source}")
    return 0 if ok else 1


# ------------------------------- batch (snapshot) ----------------------------
# The FKS1 analog of recall_erg11.py's fill/plan/prevalence, over fks1_call /
# fks1_resistance_source. `pending:sra-fks1-recaller` is the fill-eligible pending state
# (a run_acc + no call yet); rows without a run stay `pending:no-reads`.

FKS1_PENDING = "pending:sra-fks1-recaller"


def _pending_with_run(records):
    """The rows `fill` would attempt, in fill's order: pending:sra-fks1-recaller + a run_acc."""
    return [r for r in records
            if (r.get("fks1_resistance_source") or "") == FKS1_PENDING
            and (r.get("run_acc") or "").strip()]


def _fill_order(records, limit=0, shuffle=False, seed=0):
    """The exact rows a `fill` will attempt, in the order it will attempt them (shared by
    cmd_fill and cmd_plan so the offline preflight matches the real run).

    Snapshot order tracks NCBI accession/submission order, which correlates with
    outbreak/study, so the head is a biased sample. With `shuffle`, deterministically permute
    (seeded) BEFORE slicing, so a `--limit N` is a representative sample of the pending pool --
    the same discipline as recall_erg11.py and the pinned pre-registrations."""
    rows = _pending_with_run(records)
    if shuffle:
        random.Random(seed).shuffle(rows)
    return rows[:limit] if limit else rows


def cmd_fill(args):
    """Fill every pending:sra-fks1-recaller row of a snapshot in place (orchestration)."""
    _require_tools()
    reference = F.load_reference()
    records = ew.read_snapshot(args.snapshot)   # migrates a pre-v2 snapshot on load
    pending = _fill_order(records, args.limit, args.random, args.seed)
    how = f"random sample, seed={args.seed}" if args.random else "snapshot order"
    print(f"{len(pending)} isolate(s) to re-call (pending:sra-fks1-recaller with a run_acc; "
          f"{how})", file=sys.stderr)
    # A `fill` overwrites calls in place, so without a log a re-run erases its predecessor with
    # no trace. Log every fill (append-only) with before/after snapshot digests, the
    # per-isolate transcript, and full tool/code provenance -- written on exit even if a
    # mid-batch isolate crashes, so the digests then show it never wrote.
    digest_before = ew.snapshot_sha256(records)
    n_called = n_partial = n_failed = 0
    with runlog.record_run("fill-fks1", reference_path=F.DEFAULT_REFERENCE,
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
            call, source = _recall_one(run_acc, reference, F.DEFAULT_REFERENCE)
            rec["fks1_call"] = call
            rec["fks1_resistance_source"] = source
            # 'failed' is a whole-isolate failure; per-window 'refused'/'partial' live inside
            # the source string, so classify on the presence of those window tokens.
            if "failed(" in source:
                n_failed += 1
            elif "partial(" in source or "refused(" in source:
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
    print(f"summary: {n_called} with substitution(s), {n_partial} partial/refused, "
          f"{n_failed} failed, of {len(pending)} attempted", file=sys.stderr)
    return 0


def cmd_plan(args):
    """Read-only preflight: what a `fill` WOULD fetch and re-call. No tools, no network.

    Runnable anywhere (incl. osx-arm64, where the orchestration cannot run), so a snapshot can
    be inspected before committing to the multi-GB downloads. Mirrors cmd_fill's row selection
    exactly, so its counts ARE the next fill's workload."""
    records = ew.read_snapshot(args.snapshot)
    pending = [r for r in records
               if (r.get("fks1_resistance_source") or "") == FKS1_PENDING]
    with_run = _pending_with_run(records)
    without_run = len(pending) - len(with_run)
    attempted = _fill_order(records, args.limit, args.random, args.seed)
    first_runs = [r["run_acc"].split(",")[0].strip() for r in attempted]
    distinct = sorted(set(first_runs))

    print(f"snapshot:                {args.snapshot}")
    print(f"total isolates:          {len(records)}")
    print(f"pending:sra-fks1-recaller {len(pending)}"
          f"  ({without_run} with no run_acc -> skipped, stay pending)")
    print(f"fill would re-call:      {len(attempted)} isolate(s)"
          + (f"  (--limit {args.limit})" if args.limit else "")
          + (f"  [random sample, seed={args.seed}]" if args.random else "")
          + f"  == {len(attempted)} SRA download(s)")
    print(f"distinct run acc:        {len(distinct)}"
          + ("" if len(distinct) == len(attempted)
             else "  (some isolates share a run; each row still downloads once)"))
    if args.list_runs:
        for acc in distinct:
            print(acc)
    return 0


def cmd_prevalence(args):
    """The echinocandin event frequency over a (re-called) snapshot. No tools/network.

    Of the isolates whose FKS1 panel the caller could RESOLVE, the fraction carrying a known
    S639F/P/Y substitution, with a Wilson 95% CI. Undefined ('NOT MEASURED YET') until a
    `fill` populates fks1_call -- an empty snapshot reports 0 resolved, never a false 0%."""
    records = ew.read_snapshot(args.snapshot)
    p = bt.fks1_panel_prevalence(records)
    print(f"snapshot:          {args.snapshot}")
    print(f"total isolates:    {p['n_total']}")
    print(f"resolved panel:    {p['n_resolved']}"
          f"  (excluded: {p['excluded']['partial']} partial, "
          f"{p['excluded']['failed']} failed, {p['excluded']['pending']} pending)")
    if not p["n_resolved"]:
        print("echinocandin event frequency: NOT MEASURED YET "
              "(no isolate's FKS1 panel is resolved -- run a fill first)")
        return 0
    lo, hi = p["ci95"]
    print(f"panel-positive:    {p['n_panel_positive']} "
          f"(+{p['n_called_nonpanel']} resolved but non-panel/clade only)")
    print(f"event frequency:   {p['event_frequency']:.1%}  "
          f"Wilson 95% CI [{lo:.1%}, {hi:.1%}]")
    if p["per_mutation"]:
        print("by substitution:   "
              + ", ".join(f"{m} x{n}" for m, n in p["per_mutation"].items()))
    return 0


def cmd_migrate(args):
    """Rewrite a pre-v2 snapshot to the current schema in place (adds the FKS1 columns,
    seeded to an honest pending state from run_acc). Idempotent: a current snapshot is
    reported unchanged rather than needlessly rewritten. No tools/network."""
    sha, changed = ew.migrate_snapshot_file(args.snapshot)
    if changed:
        print(f"migrated {args.snapshot} -> current schema ({sha[:12]}...)")
    else:
        print(f"{args.snapshot} already current ({sha[:12]}...); no change")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("call",
                       help="deterministic: consensus -> windowed call (no network)")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--consensus-fasta",
                   help="full-length FKS1 CDS consensus, same length as the reference "
                        "(sliced to windows -- the length-exact fragile path)")
    g.add_argument("--window", action="append", metavar="NAME=PATH",
                   help="per-window consensus FASTA, e.g. --window HS1=hs1.fasta "
                        "(repeatable; a window not given is reported 'missing')")
    p.set_defaults(func=cmd_call)

    p = sub.add_parser("recall",
                       help="orchestrated: one SRA run accession -> windowed call")
    p.add_argument("run_acc", help="SRA Run accession (e.g. SRR1234567)")
    p.add_argument("--keep-tmp", action="store_true", help="keep the alignment workdir")
    p.set_defaults(func=cmd_recall)

    p = sub.add_parser("fill",
                       help="orchestrated: fill a snapshot's pending FKS1 rows in place")
    p.add_argument("snapshot", help="snapshot TSV to update in place")
    p.add_argument("--limit", type=int, default=0, help="only re-call N pending isolates")
    p.add_argument("--random", action="store_true",
                   help="sample the --limit slice randomly (seeded), not from snapshot order "
                        "-- for a representative frequency, since snapshot order is study-biased")
    p.add_argument("--seed", type=int, default=0,
                   help="RNG seed for --random (default 0); same snapshot+seed -> same sample")
    p.add_argument("--dry-run", action="store_true",
                   help="re-call but do not write the snapshot")
    p.set_defaults(func=cmd_fill)

    p = sub.add_parser("plan",
                       help="offline preflight: what a fill would fetch (no tools/network)")
    p.add_argument("snapshot", help="snapshot TSV to inspect (not modified)")
    p.add_argument("--limit", type=int, default=0, help="preview fill's N-isolate slice")
    p.add_argument("--random", action="store_true",
                   help="preview the seeded random sample fill --random would take")
    p.add_argument("--seed", type=int, default=0, help="RNG seed for --random (default 0)")
    p.add_argument("--list-runs", action="store_true",
                   help="print the distinct SRA accessions")
    p.set_defaults(func=cmd_plan)

    p = sub.add_parser("prevalence",
                       help="offline: echinocandin event frequency over resolved isolates")
    p.add_argument("snapshot", help="snapshot TSV to read (not modified)")
    p.set_defaults(func=cmd_prevalence)

    p = sub.add_parser("migrate",
                       help="rewrite a pre-v2 snapshot to the current schema in place")
    p.add_argument("snapshot", help="snapshot TSV to migrate in place")
    p.set_defaults(func=cmd_migrate)

    args = ap.parse_args()
    return args.func(args) or 0


if __name__ == "__main__":
    sys.exit(main())
