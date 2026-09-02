"""Sanity-run the FKS1 (echinocandin) re-caller against a fixed, published truth set.

The v2 analog of scripts/recaller_sanity.py (ERG11). The de-risking slice named in TODOS.md
(FKS1 v2, OPEN step 1) and in work/RUNBOOK_fks1_run.md (step 3): before running the FKS1
reads->window-consensus orchestration on a 29k-isolate snapshot, prove it on a handful of
strains whose FKS1 echinocandin-resistance genotype is already published. For each strain in
`data/earlywarning/recaller_sanity/fks1_known_genotypes.tsv`, this:

  1. downloads its SRA reads and builds a per-hot-spot-window consensus  (reuses the exact
     orchestration in scripts/recall_fks1.py -- same tools, same alignment, same per-window
     `samtools consensus -r`, so a PASS here vouches for the real `fill` path), then
  2. calls substitutions with the tested deterministic core (openafr.fks1_caller.call_windows),
     and
  3. checks the emitted PANEL hits (S639F/P/Y) against the strain's expected genotype.

Controls: the wild-type control is B11220 (clade II, echinocandin-susceptible) -- it must emit
NO S639 panel hit through the exact Illumina path, so a spurious hit there means the
orchestration (alignment/consensus), not the science, is wrong. The three positive controls
(B11211 S639F, B12136 S639P, B12663 S639Y) are drawn from the published FKS1 benchmark dataset
(PMC12323592) and cover all three panel tokens. Every fixture run is Illumina PAIRED WGS
(verified against ENA; see the fixture header), matching the real `fill` path's `-ax sr` preset.

Gated on the same bioinformatics tools as `recall`/`fill` (fasterq-dump, minimap2,
samtools>=1.13) and needs network + multi-GB downloads, so like that path it is intentionally
outside the tested core. Start here, then widen to `scripts/recall_fks1.py fill`.

Outcomes per strain:
  PASS          expected panel token(s) present (extra panel hits allowed and reported);
                wild-type strains emit no panel hit.
  FAIL          expected token missing, or an unexpected panel hit on a wild-type strain.
  INCONCLUSIVE  the panel window was uncalled (coverage gap), refused (an in-window indel), or
                missing (no consensus), or the orchestration/fetch failed -- honestly not a pass
                or a fail; the caller correctly refused to guess. Try more reads / lower depth.

Usage:
  python scripts/recaller_sanity_fks1.py                 # all rows
  python scripts/recaller_sanity_fks1.py --only B11220     # one strain (repeatable)
  python scripts/recaller_sanity_fks1.py --limit 1        # first row only (cheap smoke test)
  python scripts/recaller_sanity_fks1.py --keep-tmp       # keep alignment workdirs to inspect

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
sys.path.insert(0, str(HERE))                 # scripts/ -> recall_fks1
from openafr import fks1_caller as F            # noqa: E402
from openafr import runlog                       # noqa: E402  (append-only per-run provenance)
import recall_fks1 as orch                       # noqa: E402  (reads->window consensus + tool gate)

DEFAULT_FIXTURE = (HERE.parent / "data" / "earlywarning" / "recaller_sanity"
                   / "fks1_known_genotypes.tsv")

# The windows whose panel is non-empty (i.e. carry a known-resistance token) -- only these can
# make or break a positive expectation, so only their status gates INCONCLUSIVE. HS1 today.
PANEL_WINDOWS = [name for name, w in F.FKS1_WINDOWS.items() if w.panel]


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
    """(outcome, detail) for one strain given its call_windows() result.

    Mirrors the ERG11 sanity evaluator, but the FKS1 core carries per-window refusal/coverage
    INSIDE the result (not as an exception), so a panel window that was refused/missing/uncalled
    is INCONCLUSIVE here rather than a science FAIL -- the caller could not assess the panel."""
    got_panel = set(result["panel_hits"])
    # A panel residue explicitly uncalled (coverage gap) -> not a pass or a fail.
    if result["uncalled_panel"]:
        gap = "+".join(str(p) for p in result["uncalled_panel"])
        return "INCONCLUSIVE", f"panel residue(s) uncalled (coverage gap at {gap})"
    # A panel-bearing window that could not be placed at all (indel refusal / no consensus).
    blocked = [name for name in PANEL_WINDOWS
               if result["windows"][name]["status"] in ("refused", "missing")]
    if blocked:
        why = "; ".join(
            f"{name}={result['windows'][name]['status']}"
            f"({result['windows'][name]['refused']})"
            if result["windows"][name]["refused"]
            else f"{name}={result['windows'][name]['status']}"
            for name in blocked
        )
        return "INCONCLUSIVE", f"panel window not assessable ({why})"
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


def _covered(result):
    """Coverage across all windows, e.g. '18/18 residues' (sum over the two hot-spots)."""
    cov = sum(w["n_covered"] for w in result["windows"].values())
    tot = sum(w["n_residues"] for w in result["windows"].values())
    return f"{cov}/{tot} residues"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fixture", default=str(DEFAULT_FIXTURE),
                    help="known-genotype TSV (default: the checked-in FKS1 sanity set)")
    ap.add_argument("--only", action="append", default=[],
                    help="restrict to these strain name(s); repeatable")
    ap.add_argument("--limit", type=int, default=0, help="only run the first N rows")
    ap.add_argument("--keep-tmp", action="store_true",
                    help="keep each strain's alignment workdir for inspection")
    args = ap.parse_args()

    orch._require_tools()                        # same gate as recall/fill
    reference = F.load_reference()

    rows = _read_fixture(args.fixture)
    if args.only:
        wanted = set(args.only)
        rows = [r for r in rows if r["strain"] in wanted]
    if args.limit:
        rows = rows[:args.limit]
    if not rows:
        sys.exit("no fixture rows selected")

    print(f"FKS1 re-caller sanity run -- {len(rows)} strain(s), "
          f"reference {os.path.basename(F.DEFAULT_REFERENCE)}\n", file=sys.stderr)

    n_pass = n_fail = n_incon = 0
    # Every sanity run is logged (append-only) with its full provenance and per-strain outcomes,
    # same as the ERG11 harness -- a PASS/FAIL is a durable, git-diffable record. The record is
    # written on exit even if a strain crashes mid-loop.
    with runlog.record_run("sanity-fks1", reference_path=F.DEFAULT_REFERENCE,
                           tools=list(orch.TOOLS),
                           extra={"fixture": os.path.basename(args.fixture),
                                  "n_strains": len(rows)}) as run:
        for i, row in enumerate(rows, 1):
            strain, run_acc = row["strain"], row["run_acc"].strip()
            expected = _expected_set(row.get("expected_panel"))
            exp_str = ",".join(sorted(expected)) or "wild-type"
            print(f"[{i}/{len(rows)}] {strain} (clade {row.get('clade','?')}) "
                  f"{run_acc}  expect: {exp_str}", file=sys.stderr)

            tmp = tempfile.mkdtemp(prefix=f"sanity_fks1_{strain}_")
            try:
                window_seqs = orch.reads_to_window_consensus(
                    run_acc, tmp, F.DEFAULT_REFERENCE)
                result = F.call_windows(window_seqs, reference=reference)
            except Exception as e:                # tool/fetch failure -> not a science FAIL
                outcome, detail = "INCONCLUSIVE", f"orchestration failed: {orch._short(e)}"
                result = None
            else:
                outcome, detail = _evaluate(expected, result)
            finally:
                if not args.keep_tmp:
                    shutil.rmtree(tmp, ignore_errors=True)

            call = ",".join(result["tokens"]) if result is not None else ""
            if result is not None:
                print(f"        call: {call or '(none)'}   covered: {_covered(result)}",
                      file=sys.stderr)
            print(f"        {outcome}: {detail}\n", file=sys.stderr)

            run.add_item(strain=strain, clade=row.get("clade", ""), run_acc=run_acc,
                         expected=exp_str, call=call, outcome=outcome, detail=detail,
                         covered=(_covered(result) if result is not None else None))
            n_pass += outcome == "PASS"
            n_fail += outcome == "FAIL"
            n_incon += outcome == "INCONCLUSIVE"

        run.set(status=("fail" if n_fail else "ok"),
                summary={"pass": n_pass, "fail": n_fail, "inconclusive": n_incon,
                         "total": len(rows)})

    print(f"summary: {n_pass} PASS, {n_fail} FAIL, {n_incon} INCONCLUSIVE "
          f"of {len(rows)}", file=sys.stderr)
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
