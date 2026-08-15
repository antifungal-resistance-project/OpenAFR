"""Sanity-run the ERG11 re-caller against a fixed, published known-genotype truth set.

The de-risking slice named in TODOS.md (OPEN step 1): before running the reads->consensus
orchestration on a 29k-isolate snapshot, prove it on a handful of strains whose ERG11
azole-resistance genotype is already published. For each strain in
`data/earlywarning/recaller_sanity/known_genotypes.tsv`, this:

  1. downloads its SRA reads and builds an ERG11 consensus CDS  (reuses the exact
     orchestration in scripts/recall_erg11.py -- same tools, same alignment, same
     consensus params, so a PASS here vouches for the real `fill` path), then
  2. calls substitutions with the tested deterministic core (openafr.recaller), and
  3. checks the emitted PANEL hits against the strain's expected clade genotype.

The strongest row is B8441: the caller's reference CDS (XM_085597798.1) is derived from
B8441, so B8441's own reads must call WILD-TYPE -- if that row shows spurious panel hits,
the orchestration (alignment/consensus), not the science, is wrong.

Gated on the same bioinformatics tools as `recall`/`fill` (fasterq-dump, minimap2,
samtools) and needs network + multi-GB downloads, so like that path it is intentionally
outside the tested core. Start here, then widen to `scripts/recall_erg11.py fill`.

Outcomes per strain:
  PASS          expected panel token(s) present (extra panel hits, e.g. a co-occurring
                V125A, are allowed and reported); wild-type strains emit no panel hit.
  FAIL          expected token missing, or an unexpected panel hit on a wild-type strain.
  INCONCLUSIVE  a panel residue was uncalled (coverage gap) -- honestly not a pass or a
                fail; the caller correctly refused to guess. Try more reads / lower depth.

Usage:
  python scripts/recaller_sanity.py                 # all rows
  python scripts/recaller_sanity.py --only B8441     # one strain (repeatable)
  python scripts/recaller_sanity.py --limit 1        # first row only (cheap smoke test)
  python scripts/recaller_sanity.py --keep-tmp       # keep alignment workdirs to inspect

Exit code is nonzero if any row FAILs (INCONCLUSIVE does not fail the run).
"""
import argparse
import os
import pathlib
import shutil
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))          # repo root -> openafr
sys.path.insert(0, str(HERE))                 # scripts/ -> recall_erg11
from openafr import recaller as R              # noqa: E402
import recall_erg11 as orch                    # noqa: E402  (reads->consensus + tool gate)

DEFAULT_FIXTURE = (HERE.parent / "data" / "earlywarning" / "recaller_sanity"
                   / "known_genotypes.tsv")


def _read_fixture(path):
    rows = []
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                continue
            rows.append(dict(zip(header, fields)))
    return rows


def _expected_set(cell):
    cell = (cell or "").strip()
    if cell in ("", "-"):
        return set()
    return {t.strip() for t in cell.split(",") if t.strip()}


def _evaluate(expected, result):
    """(outcome, detail) for one strain given its call_substitutions() result."""
    got_panel = set(result["panel_hits"])
    if result["uncalled_panel"]:
        gap = "+".join(str(p) for p in result["uncalled_panel"])
        return "INCONCLUSIVE", f"panel residue(s) uncalled (coverage gap at {gap})"
    if expected:
        missing = expected - got_panel
        if missing:
            return "FAIL", f"expected panel {sorted(expected)} missing {sorted(missing)}"
        extra = got_panel - expected
        extra_note = f"; also saw panel {sorted(extra)}" if extra else ""
        return "PASS", f"panel {sorted(got_panel)} covers expected{extra_note}"
    # wild-type expectation: no panel hit allowed
    if got_panel:
        return "FAIL", f"expected wild-type but got panel hit(s) {sorted(got_panel)}"
    return "PASS", "no panel substitution (wild-type), as expected"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fixture", default=str(DEFAULT_FIXTURE),
                    help="known-genotype TSV (default: the checked-in sanity set)")
    ap.add_argument("--only", action="append", default=[],
                    help="restrict to these strain name(s); repeatable")
    ap.add_argument("--limit", type=int, default=0, help="only run the first N rows")
    ap.add_argument("--keep-tmp", action="store_true",
                    help="keep each strain's alignment workdir for inspection")
    args = ap.parse_args()

    orch._require_tools()                        # same gate as recall/fill
    reference = R.load_reference()

    rows = _read_fixture(args.fixture)
    if args.only:
        wanted = set(args.only)
        rows = [r for r in rows if r["strain"] in wanted]
    if args.limit:
        rows = rows[:args.limit]
    if not rows:
        sys.exit("no fixture rows selected")

    print(f"ERG11 re-caller sanity run -- {len(rows)} strain(s), "
          f"reference {os.path.basename(R.DEFAULT_REFERENCE)}\n", file=sys.stderr)

    n_pass = n_fail = n_incon = 0
    for i, row in enumerate(rows, 1):
        strain, run_acc = row["strain"], row["run_acc"].strip()
        expected = _expected_set(row.get("expected_panel"))
        exp_str = ",".join(sorted(expected)) or "wild-type"
        print(f"[{i}/{len(rows)}] {strain} (clade {row.get('clade','?')}) "
              f"{run_acc}  expect: {exp_str}", file=sys.stderr)

        tmp = tempfile.mkdtemp(prefix=f"sanity_{strain}_")
        try:
            cons_path = orch.reads_to_consensus(run_acc, tmp, R.DEFAULT_REFERENCE)
            consensus = orch._read_consensus_fasta(cons_path)
            result = R.call_substitutions(consensus, reference=reference)
        except R.ConsensusError as e:
            outcome, detail = "INCONCLUSIVE", f"consensus refused: {orch._short(e)}"
            result = None
        except Exception as e:                   # tool/fetch failure -> not a science FAIL
            outcome, detail = "INCONCLUSIVE", f"orchestration failed: {orch._short(e)}"
            result = None
        else:
            outcome, detail = _evaluate(expected, result)
        finally:
            if not args.keep_tmp:
                shutil.rmtree(tmp, ignore_errors=True)

        if result is not None:
            call = ",".join(result["tokens"]) or "(none)"
            cov = f"{result['n_covered']}/{result['n_residues']} residues"
            print(f"        call: {call}   covered: {cov}", file=sys.stderr)
        print(f"        {outcome}: {detail}\n", file=sys.stderr)

        n_pass += outcome == "PASS"
        n_fail += outcome == "FAIL"
        n_incon += outcome == "INCONCLUSIVE"

    print(f"summary: {n_pass} PASS, {n_fail} FAIL, {n_incon} INCONCLUSIVE "
          f"of {len(rows)}", file=sys.stderr)
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
