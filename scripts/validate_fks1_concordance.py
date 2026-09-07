"""Grade the FKS1 re-caller against the full 100-isolate external benchmark (concordance run).

The graded version of scripts/recaller_sanity_fks1.py. Where the sanity harness runs a handful
of strains as PASS/FAIL controls, this runs the whole pre-registered truth set
(work/PREREGISTRATION_fks1_concordance.md, sha frozen in work/PREREG_fks1_concordance.sha256)
and computes the caller's measured sensitivity / specificity / per-token identity with Wilson
CIs via openafr.concordance, then applies the pre-committed pass bar.

It reuses the EXACT production orchestration (recall_fks1.reads_to_window_consensus ->
fks1_caller.call_windows), so a PASS certifies the real `fill` path, not a parallel one.

Integrity gate (refuses to certify if any of these moved from freeze):
  * this run's pre-registration document hash,
  * the pinned reference FKS1 CDS hash (load_reference already asserts this), and
  * the truth-set fixture hash.
Pins live in work/PREREG_fks1_concordance.sha256 (prereg + fixture) and are re-checked here.

Truth-set fixture (verified manual transcription of PMC12323592 Table 1; NOT scraped):
  data/earlywarning/recaller_sanity/fks1_benchmark_100.tsv
Columns (tab-separated, '#'-comments allowed):
  strain  clade  run_acc  expected_panel  susceptibility  note  citation
    expected_panel  HS1 panel token(s) the paper attributes (comma-sep), '-' if none/non-HS1
    susceptibility  'R' or 'S' (echinocandin) -- the phenotype anchor for resistance-detection

Gated on the same bioinformatics tools + network + multi-GB downloads as `recall`/`fill`, so
like those it is intentionally outside the tested core; openafr.concordance (the scoring) IS
unit-tested offline (tests/test_concordance.py). Meant for the Linux/x86 cloud host, per
work/RUNBOOK_fks1_run.md.

Usage:
  python scripts/validate_fks1_concordance.py                 # full benchmark, then verdict
  python scripts/validate_fks1_concordance.py --limit 5       # first N rows (cheap smoke test)
  python scripts/validate_fks1_concordance.py --keep-tmp      # keep alignment workdirs
"""
import argparse
import hashlib
import os
import pathlib
import shutil
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))
from openafr import fks1_caller as F        # noqa: E402
from openafr import concordance as cc        # noqa: E402
from openafr import runlog                    # noqa: E402
import recall_fks1 as orch                    # noqa: E402  (reads->window consensus + tool gate)

PREREG = ROOT / "work" / "PREREGISTRATION_fks1_concordance.md"
PINS = ROOT / "work" / "PREREG_fks1_concordance.sha256"
DEFAULT_FIXTURE = ROOT / "data" / "earlywarning" / "recaller_sanity" / "fks1_benchmark_100.tsv"

PANEL_WINDOWS = [name for name, w in F.FKS1_WINDOWS.items() if w.panel]

# Pre-committed pass bar (work/PREREGISTRATION_fks1_concordance.md, frozen 2026-09-07).
BAR = {
    "specificity_point": 0.95, "specificity_lo": 0.90,
    "sensitivity_point": 0.90, "sensitivity_lo": 0.75,
    "identity_point": 0.95,
    "resolved_frac": 0.80,
}


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _check_integrity(fixture):
    """Re-check the frozen prereg + fixture hashes. Returns (ok, messages)."""
    msgs, ok = [], True
    pinned = {}
    if PINS.exists():
        for line in PINS.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            digest, _, name = line.partition("  ")
            pinned[name.strip()] = digest.strip()
    for target in (PREREG, fixture):
        rel = os.path.relpath(target, ROOT)
        want = pinned.get(rel)
        got = _sha256(target) if target.exists() else None
        if got is None:
            ok = False; msgs.append(f"MISSING: {rel} does not exist")
        elif want is None:
            ok = False; msgs.append(f"UNPINNED: {rel} has no frozen hash in {PINS.name} "
                                    f"(freeze it: `shasum -a 256 {rel} >> {PINS.name}`)")
        elif want != got:
            ok = False; msgs.append(f"HASH MOVED: {rel}\n    frozen {want}\n    now    {got}")
        else:
            msgs.append(f"ok: {rel} matches frozen hash")
    return ok, msgs


def _read_fixture(path):
    rows, header = [], None
    with open(path) as fh:
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
    return set() if cell in ("", "-") else {t.strip() for t in cell.split(",") if t.strip()}


def _resolved(result):
    """True if the panel was ASSESSABLE (mirrors the sanity harness INCONCLUSIVE logic).

    Unresolved = a panel residue uncalled (coverage gap) OR a panel-bearing window
    refused/missing (indel refusal / no consensus). Such isolates are excluded from every
    concordance denominator, never scored as a negative."""
    if result is None:
        return False
    if result["uncalled_panel"]:
        return False
    if any(result["windows"][n]["status"] in ("refused", "missing") for n in PANEL_WINDOWS):
        return False
    return True


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fixture", default=str(DEFAULT_FIXTURE))
    ap.add_argument("--limit", type=int, default=0, help="only run the first N rows")
    ap.add_argument("--keep-tmp", action="store_true")
    ap.add_argument("--allow-hash-mismatch", action="store_true",
                    help="run anyway if integrity fails (prints UNCERTIFIED; for debugging)")
    args = ap.parse_args()

    fixture = pathlib.Path(args.fixture)
    ok, msgs = _check_integrity(fixture)
    print("integrity:", file=sys.stderr)
    for m in msgs:
        print("  " + m, file=sys.stderr)
    if not ok and not args.allow_hash_mismatch:
        sys.exit("\nREFUSING TO CERTIFY: integrity check failed (see above). The truth-set "
                 "fixture must exist and its hash + the prereg hash must match "
                 f"{PINS.name}. Freeze the fixture before running, or pass "
                 "--allow-hash-mismatch to dry-run uncertified.")

    orch._require_tools()
    reference = F.load_reference()
    rows = _read_fixture(fixture)
    if args.limit:
        rows = rows[:args.limit]
    if not rows:
        sys.exit("no fixture rows selected")

    print(f"\nFKS1 concordance run -- {len(rows)} isolate(s), reference "
          f"{os.path.basename(F.DEFAULT_REFERENCE)}\n", file=sys.stderr)

    evals = []
    with runlog.record_run("concordance-fks1", reference_path=F.DEFAULT_REFERENCE,
                           tools=list(orch.TOOLS),
                           extra={"fixture": os.path.basename(str(fixture)),
                                  "n_isolates": len(rows),
                                  "prereg_sha": _sha256(PREREG)}) as run:
        for i, row in enumerate(rows, 1):
            strain, run_acc = row["strain"], row["run_acc"].strip()
            expected = _expected_set(row.get("expected_panel"))
            resistant = (row.get("susceptibility", "").strip().upper() == "R")
            print(f"[{i}/{len(rows)}] {strain} {run_acc}  expect: "
                  f"{','.join(sorted(expected)) or 'wild-type'}  pheno: "
                  f"{'R' if resistant else 'S'}", file=sys.stderr)

            tmp = tempfile.mkdtemp(prefix=f"conc_fks1_{strain}_")
            try:
                window_seqs = orch.reads_to_window_consensus(
                    run_acc, tmp, F.DEFAULT_REFERENCE)
                result = F.call_windows(window_seqs, reference=reference)
            except Exception as e:
                result = None
                print(f"        orchestration failed: {orch._short(e)}", file=sys.stderr)
            finally:
                if not args.keep_tmp:
                    shutil.rmtree(tmp, ignore_errors=True)

            resolved = _resolved(result)
            called = set(result["panel_hits"]) if result is not None else set()
            evals.append({"expected": expected, "called": called,
                          "resolved": resolved, "resistant": resistant})
            run.add_item(strain=strain, run_acc=run_acc,
                         expected=",".join(sorted(expected)) or "-",
                         called=",".join(sorted(called)) or "-",
                         resolved=resolved, resistant=resistant)
            print(f"        called: {','.join(sorted(called)) or '(none)'}  "
                  f"resolved: {resolved}\n", file=sys.stderr)

        summary = cc.evaluate(evals)
        verdict, reasons = _verdict(summary, certified=ok)
        run.set(status=("ok" if verdict == "PASS" else "fail"),
                summary={"verdict": verdict, "confusion": summary["confusion"],
                         "n_resolved": summary["n_resolved"],
                         "n_total": summary["n_total"]})

    print("\n" + cc.format_report(summary))
    print("\npre-registered verdict: " + verdict)
    for r in reasons:
        print("  " + r)
    return 0 if verdict == "PASS" else 1


def _verdict(summary, certified):
    """Apply the pre-committed pass bar. Returns (verdict, [reason lines])."""
    reasons = []
    n_total = summary["n_total"] or 1
    resolved_frac = summary["n_resolved"] / n_total
    if resolved_frac < BAR["resolved_frac"]:
        return "UNDERPOWERED", [
            f"resolved {summary['n_resolved']}/{summary['n_total']} = {resolved_frac:.0%} "
            f"< {BAR['resolved_frac']:.0%} required -- coverage too thin to certify; "
            f"try more reads / a lower min-depth. Neither a pass nor a fail."]
    if not certified:
        reasons.append("NOTE: uncertified (integrity check was overridden).")

    sens, spec, ident = summary["sensitivity"], summary["specificity"], summary["identity"]
    checks = [
        ("specificity", spec["point"], BAR["specificity_point"],
         spec["ci95"][0] if spec["ci95"] else None, BAR["specificity_lo"]),
        ("sensitivity", sens["point"], BAR["sensitivity_point"],
         sens["ci95"][0] if sens["ci95"] else None, BAR["sensitivity_lo"]),
        ("identity", ident["point"], BAR["identity_point"], None, None),
    ]
    passed = True
    for name, point, pbar, lo, lobar in checks:
        if point is None:
            passed = False; reasons.append(f"{name}: no resolved isolates to measure -- FAIL")
            continue
        pt_ok = point >= pbar
        lo_ok = (lobar is None) or (lo is not None and lo >= lobar)
        passed = passed and pt_ok and lo_ok
        detail = f"point {point:.1%} (bar {pbar:.0%})"
        if lobar is not None:
            detail += f", lower {lo:.1%} (bar {lobar:.0%})"
        reasons.append(f"{name}: {detail} -- {'ok' if pt_ok and lo_ok else 'BELOW BAR'}")
    return ("PASS" if passed else "FAIL"), reasons


if __name__ == "__main__":
    sys.exit(main())
