#!/usr/bin/env python3
"""Interpret one *C. auris* isolate's genotype -> an RUO resistance-marker verdict (issue #163).

The command-line front door to the genotype->verdict engine. The engine itself already exists
as library code -- openafr/interpret.py turns a genotype into a #133 verdict, openafr/report.py
renders one isolate's per-drug-class verdicts into a human brief or machine JSON. This script is
the thin glue a lab user actually runs: it collects the handful of raw input shapes a WGS /
targeted-panel pipeline produces, calls interpret() once per drug-class gene, and renders the
combined report. It re-implements NO calling or interpretation logic -- every base/codon/pocket
decision stays in the pinned-reference callers behind interpret().

A lab user holds one of three input shapes per gene; pick whichever their pipeline emits:

  * variant tokens -- the easiest path, no files, no external tools:
        interpret_isolate.py --isolate B12345 --gene FKS1 --variants S639P
  * a full-length consensus CDS FASTA (ERG11: in-frame; FKS1: the fragile length-exact bridge):
        interpret_isolate.py --isolate B12345 --gene ERG11 --cds erg11_consensus.fasta
  * one consensus FASTA per FKS1 hot-spot window (the robust path `recall_fks1.py recall` emits):
        interpret_isolate.py --isolate B12345 --gene FKS1 --window HS1=hs1.fa --window HS2=hs2.fa

Type an isolate for BOTH drug-classes in one call by repeating --gene with its own input:
    interpret_isolate.py --isolate B12345 --gene FKS1 --variants S639P --gene ERG11 --variants Y132F

Everything here is RUO (research use only), not a clinical determination, and checks only the
named ERG11/FKS1 markers -- see docs/QUICKSTART_LAB.md for how to read each verdict and
docs/DIAGNOSTIC_VERDICT_CONTRACT.md for the full contract. Stdlib only, offline; the
variant-token path needs no conda env and no external binaries.
"""
import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import report as _report       # noqa: E402  (build_report / render_*)
from openafr.interpret import interpret     # noqa: E402  (genotype -> verdict)


def _read_consensus_fasta(path):
    """Read a single-record consensus FASTA -> the bare sequence string. Refuses empty /
    header-only / multi-record input, matching scripts/recall_fks1.py's reader so an upstream
    hiccup can never read as a silently wild-type empty sequence."""
    seqs, cur = [], []
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


def _parse_windows(specs):
    """--window NAME=PATH ... -> {NAME: consensus_seq}. interpret() validates the window names
    (FKS1-only, known hot-spots); here we only parse the NAME=PATH shape and read each file."""
    windows = {}
    for spec in specs:
        name, sep, path = spec.partition("=")
        name = name.strip()
        if not sep or not path:
            sys.exit(f"--window expects NAME=PATH, got {spec!r}")
        windows[name] = _read_consensus_fasta(path)
    return windows


def _interpret_gene(args):
    """Turn one parsed --gene block into a verdict object, delegating every call/validation
    decision to interpret(). Exactly-one-input and token/gene validation live in interpret();
    we surface its ValueError as a clean exit rather than a traceback."""
    kwargs = {}
    if args.variants is not None:
        kwargs["variants"] = args.variants
    if args.cds is not None:
        kwargs["cds"] = _read_consensus_fasta(args.cds)
    if args.window:
        kwargs["windows"] = _parse_windows(args.window)
    if args.uncalled:
        kwargs["uncalled"] = args.uncalled
    try:
        return interpret(args.gene, **kwargs)
    except ValueError as e:
        sys.exit(f"{args.gene}: {e}")


def _build_gene_parser():
    """A sub-parser for one --gene block's inputs, reused for every repeat of --gene so each
    drug-class carries its own one-of {variants, cds, windows}."""
    gp = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    gp.add_argument("--gene", required=True,
                    help="ERG11 (azole) or FKS1 (echinocandin); CYP51A is refused (no caller yet)")
    gp.add_argument("--variants", nargs="*",
                    help="protein substitution tokens from a variant caller, e.g. S639P Y132F "
                         "(an empty list is a valid 'nothing found' input)")
    gp.add_argument("--cds", help="full-length consensus CDS FASTA (single record)")
    gp.add_argument("--window", action="append", default=[], metavar="NAME=PATH",
                    help="FKS1 hot-spot consensus FASTA, e.g. HS1=hs1.fa (repeatable)")
    gp.add_argument("--uncalled", nargs="*", type=int, default=[], metavar="POS",
                    help="residue positions NOT covered; a panel residue here forces UNRESOLVED")
    return gp


def _split_gene_blocks(argv, gene_parser):
    """Partition argv (already stripped of the top-level flags) into one block per --gene so
    --variants/--cds/--window after each --gene bind to that gene. A single pass of argparse
    cannot express 'repeat this whole group', so we slice on --gene ourselves. Anything before
    the first --gene is a stray token the top parser did not recognize -- reported, not ignored."""
    idxs = [i for i, a in enumerate(argv) if a == "--gene"]
    if not idxs:
        sys.exit("at least one --gene block is required (e.g. --gene FKS1 --variants S639P)")
    if idxs[0] != 0:
        sys.exit(f"unrecognized argument(s) before the first --gene: {argv[:idxs[0]]}")
    bounds = idxs + [len(argv)]
    blocks = [argv[bounds[i]:bounds[i + 1]] for i in range(len(idxs))]
    return [gene_parser.parse_args(b) for b in blocks]


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    top = argparse.ArgumentParser(
        prog="interpret_isolate.py", description=__doc__, allow_abbrev=False,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    top.add_argument("--isolate", required=True, help="isolate label / accession for the report")
    top.add_argument("--json", action="store_true",
                     help="emit machine-readable JSON instead of the human markdown brief")
    top.add_argument("--operator", help="optional operator/run tag, surfaced in report metadata")

    # Pull the top-level flags out wherever they sit, so --json/--isolate aren't swallowed by a
    # preceding --gene block; the rest is purely gene blocks, sliced on --gene.
    args, rest = top.parse_known_args(argv)
    gene_parser = _build_gene_parser()
    gene_blocks = _split_gene_blocks(rest, gene_parser)

    verdicts = [_interpret_gene(g) for g in gene_blocks]

    metadata = {"operator": args.operator} if args.operator else None
    report = _report.build_report(args.isolate, verdicts, metadata=metadata)
    print(_report.render_json(report) if args.json else _report.render_markdown(report))
    return 0


if __name__ == "__main__":
    sys.exit(main())
