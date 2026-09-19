"""Score the diagnostics VERDICT against phenotypic AST truth, per drug class (#137).

The graded accuracy run governed by work/PREREGISTRATION_diagnostic_accuracy.md (sha frozen
in work/PREREG_diagnostic_accuracy.sha256). Where the concordance grader
(scripts/validate_fks1_concordance.py) asks "is the caller's TOKEN right?" against a genotype
truth set, this asks the clinically loaded question against a PHENOTYPE truth set: when the
engine emits a verdict, how often does it agree with the isolate's R/S, and how often does it
make the two dangerous errors -- a VERY MAJOR error (called NO_KNOWN_MARKER on a resistant
isolate) or a MAJOR error (called RESISTANCE_MARKER_DETECTED on a susceptible one). It reports
VME/ME/sensitivity/specificity/abstention per drug class (openafr.accuracy) and applies the
pre-committed pass bar.

  *** THIS RUN IS BLOCKED ON DATA AND DOES NOT CERTIFY ANYTHING TODAY. ***

Scoring a verdict against phenotype needs paired genotype + phenotypic AST per isolate. The
echinocandin arm reuses the 100-isolate PMC12323592 FKS1 benchmark
(data/earlywarning/recaller_sanity/fks1_benchmark_100.tsv) -- itself blocked on a verified
manual (non-scraped) transcription of the paywalled Table 1 + a GCP run. The azole arm has
only the 4 Lockhart clade strains and is declared underpowered until a larger paired ERG11
panel exists. Until a fixture is assembled and its hash pinned in
PREREG_diagnostic_accuracy.sha256, this grader refuses to run -- exactly as
validate_fks1_concordance.py / validate_calibration.py refuse without their fixtures. The
scoring math it would apply (openafr.accuracy) IS unit-tested offline (tests/test_accuracy.py).

Expected fixture schema (TSV, '#'-comments allowed) once a paired panel is frozen -- one row
per resolvable isolate:
  isolate  species  run_acc  gene  drug_class  called_verdict  susceptibility  note  citation
    called_verdict  the enum the production path (recall* -> caller -> build_verdict) emits;
                    filled by the orchestration, NOT hand-authored
    susceptibility  'R' | 'S' against a STATED CLSI/EUCAST breakpoint

Usage (once unblocked):
  python scripts/validate_diagnostic_accuracy.py --fixture <panel.tsv>
  python scripts/validate_diagnostic_accuracy.py --fixture <panel.tsv> --allow-hash-mismatch
"""
import argparse
import hashlib
import os
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from openafr import accuracy as ac        # noqa: E402

PREREG = ROOT / "work" / "PREREGISTRATION_diagnostic_accuracy.md"
PINS = ROOT / "work" / "PREREG_diagnostic_accuracy.sha256"

# Pre-committed pass bar (work/PREREGISTRATION_diagnostic_accuracy.md, frozen 2026-09-10).
BAR = {
    "vme_point": 0.03, "vme_upper": 0.15, "me_point": 0.05,
    "min_resistant": 15, "min_susceptible": 15, "max_abstention": 0.30,
}


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _check_integrity(fixture):
    """Re-check the frozen prereg hash and that the accuracy fixture exists + is pinned."""
    msgs, ok = [], True
    pinned = {}
    if PINS.exists():
        for line in PINS.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            digest, _, name = line.partition("  ")
            pinned[name.strip()] = digest.strip()
    rel = os.path.relpath(PREREG, ROOT)
    want, got = pinned.get(rel), (_sha256(PREREG) if PREREG.exists() else None)
    if got is None:
        ok = False; msgs.append(f"MISSING: {rel}")
    elif want != got:
        ok = False; msgs.append(f"HASH MOVED: {rel}\n    frozen {want}\n    now    {got}")
    else:
        msgs.append(f"ok: {rel} matches frozen hash")
    if fixture is None:
        ok = False
        msgs.append("BLOCKED: no accuracy fixture given. A paired genotype+phenotype panel "
                    "does not exist yet (see work/PREREGISTRATION_diagnostic_accuracy.md). "
                    "This run cannot certify until that fixture is assembled, transcribed "
                    "(verified, not scraped) and its hash pinned in "
                    "PREREG_diagnostic_accuracy.sha256.")
    elif not fixture.exists():
        ok = False; msgs.append(f"MISSING: fixture {fixture} does not exist")
    else:
        frel = os.path.relpath(fixture, ROOT)
        fwant, fgot = pinned.get(frel), _sha256(fixture)
        if fwant is None:
            ok = False; msgs.append(f"UNPINNED: {frel} has no frozen hash in {PINS.name} "
                                    f"(freeze it: `shasum -a 256 {frel} >> {PINS.name}`)")
        elif fwant != fgot:
            ok = False; msgs.append(f"HASH MOVED: {frel}\n    frozen {fwant}\n    now    {fgot}")
        else:
            msgs.append(f"ok: {frel} matches frozen hash")
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


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fixture", default=None, help="the frozen paired-panel TSV (required to "
                    "certify; absent -> BLOCKED)")
    ap.add_argument("--allow-hash-mismatch", action="store_true",
                    help="run anyway if integrity fails (prints UNCERTIFIED; for debugging)")
    args = ap.parse_args()

    fixture = pathlib.Path(args.fixture) if args.fixture else None
    ok, msgs = _check_integrity(fixture)
    print("integrity:", file=sys.stderr)
    for m in msgs:
        print("  " + m, file=sys.stderr)
    if not ok and not args.allow_hash_mismatch:
        sys.exit("\nREFUSING TO CERTIFY: integrity check failed (see above). This run is "
                 "BLOCKED until a paired genotype+phenotype panel exists and its hash is "
                 "pinned (work/PREREGISTRATION_diagnostic_accuracy.md). Pass "
                 "--allow-hash-mismatch to dry-run uncertified on a stand-in fixture.")

    rows = _read_fixture(fixture)
    if not rows:
        sys.exit("no fixture rows")

    records = [{"drug_class": r.get("drug_class", "").strip(),
                "verdict": r.get("called_verdict", "").strip(),
                "phenotype": r.get("susceptibility", "").strip()} for r in rows]
    summary = ac.evaluate(records)
    print(ac.format_report(summary))

    any_fail = False
    for cls in summary["drug_classes"]:
        outcome, reasons = ac.apply_bar(summary["by_drug_class"][cls], BAR)
        any_fail = any_fail or outcome == "FAIL"
        print(f"\n=== {cls} -- {outcome} ===")
        for rr in reasons:
            print("  " + rr)

    print("\n(uncertified dry-run)" if not ok else "")
    return 1 if any_fail else 0


if __name__ == "__main__":
    sys.exit(main())
