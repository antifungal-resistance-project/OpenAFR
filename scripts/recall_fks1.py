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

Not implemented here -- deferred, not silently dropped
------------------------------------------------------
`fill`/`plan` over an NCBI snapshot (the ERG11 batch path) are NOT provided. The snapshot
schema (openafr.earlywarning.SNAPSHOT_COLUMNS) carries `erg11_call`/`resistance_source` only,
and `_serialize` writes a FIXED fieldname set with `extrasaction="ignore"` -- so an `fks1_call`
put on a record would be silently discarded on write. Populating a snapshot with FKS1 calls is
a separate increment (T3: add fks1 columns to the schema, seed `pending:sra-fks1-recaller` in
the snapshot builder, migrate existing snapshots, wire the emergence/diff surface). Until then
this CLI operates per isolate (`recall`) and per consensus (`call`).

Honesty
-------
A window whose hot-spot residues could not be covered is reported per-window as
`sra-fks1-recaller:HS1=partial(uncalled:639),...` with those residues left OUT of the call --
never a silent "wild type". A window with an in-window indel is `refused(...)` for that window
alone. An isolate that fails to fetch/align yields an empty call with a `failed(<reason>)`
source, not a dropped or guessed one.
"""
import argparse
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
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

    args = ap.parse_args()
    return args.func(args) or 0


if __name__ == "__main__":
    sys.exit(main())
