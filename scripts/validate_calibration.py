"""Calibrate the resistance probability and MEASURE its reliability on a held-out split (#136).

The graded calibration run governed by work/PREREGISTRATION_calibration.md (sha frozen in
work/PREREG_calibration.sha256). Where the concordance grader
(scripts/validate_fks1_concordance.py) asks "is the caller's token right?", this asks the
next question: when the engine stamps a probability, is it reliable? It fits a stratified
probability map on a 70% training split and reports Brier / reliability curve / ECE /
calibration-in-the-large on the disjoint 30% test split (openafr.calibration), per drug and
per regime, then applies the pre-committed pass bar.

  *** THIS RUN IS BLOCKED ON DATA AND DOES NOT CERTIFY ANYTHING TODAY. ***

Per the frozen #132 gate (work/PREREGISTRATION_diagnostics_panel.md, sha 1211e973...), no
obtainable public Candida collection clears the calibration bar; the verdict was FAIL ->
concordance-only, with the calibrated-probability product deferred to "option C" -- a pooled
C. albicans ERG11 panel that needs (1) a separate frozen prereg AMENDMENT for cross-study
MIC-method heterogeneity, (2) new C. albicans caller/structure work (the port is auris/5TZ1),
and (3) a verified manual (non-scraped) transcription of the paired genotype+MIC tables into a
hashed fixture. Until that fixture exists and its hash is pinned, this grader refuses to run --
exactly as validate_fks1_concordance.py refuses without its 100-isolate fixture. The scoring
math it would apply (openafr.calibration) IS unit-tested offline (tests/test_calibration.py).

Expected fixture schema (TSV, '#'-comments allowed) once the option-C panel is assembled and
frozen -- one row per resolvable isolate:
  isolate  species  source_study  drug  stratum  regime  susceptibility  note  citation
    drug            the azole/echinocandin this row's phenotype is against
    stratum         known regime: the panel token (Y132F, ...); novel regime: the structural
                    evidence tier (e.g. 'weaker-fit/estimate') the caller emits for the variant
    regime          'known' | 'novel'
    susceptibility  'R' | 'S' against a STATED CLSI/EUCAST breakpoint (the calibration target)

Usage (once unblocked):
  python scripts/validate_calibration.py --fixture <panel.tsv>
  python scripts/validate_calibration.py --fixture <panel.tsv> --allow-hash-mismatch   # dry-run
"""
import argparse
import hashlib
import os
import pathlib
import random
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from openafr import calibration as cal        # noqa: E402

PREREG = ROOT / "work" / "PREREGISTRATION_calibration.md"
PINS = ROOT / "work" / "PREREG_calibration.sha256"

# Pre-committed pass bar (work/PREREGISTRATION_calibration.md, frozen 2026-09-09).
BAR = {"ece": 0.10, "brier": 0.20, "in_large_gap": 0.10, "min_per_class": 15}
SPLIT_SEED = 20260909          # frozen in the pre-registration
TEST_FRACTION = 0.30


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _check_integrity(fixture):
    """Re-check the frozen prereg hash and that the panel fixture exists + is pinned."""
    msgs, ok = [], True
    pinned = {}
    if PINS.exists():
        for line in PINS.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            digest, _, name = line.partition("  ")
            pinned[name.strip()] = digest.strip()
    # prereg must match its frozen hash
    rel = os.path.relpath(PREREG, ROOT)
    want, got = pinned.get(rel), (_sha256(PREREG) if PREREG.exists() else None)
    if got is None:
        ok = False; msgs.append(f"MISSING: {rel}")
    elif want != got:
        ok = False; msgs.append(f"HASH MOVED: {rel}\n    frozen {want}\n    now    {got}")
    else:
        msgs.append(f"ok: {rel} matches frozen hash")
    # panel fixture must exist and be pinned (the data block lives here)
    if fixture is None:
        ok = False
        msgs.append("BLOCKED: no panel fixture given. The option-C C. albicans calibration "
                    "panel does not exist yet (see work/PREREGISTRATION_calibration.md). This "
                    "run cannot certify until that fixture is assembled, transcribed (verified, "
                    "not scraped) and its hash pinned in PREREG_calibration.sha256.")
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


def _stratified_split(rows, seed, test_fraction):
    """Deterministic 70/30 split stratified on the R/S label (frozen seed)."""
    rng = random.Random(seed)
    by_label = {}
    for r in rows:
        by_label.setdefault(r["susceptibility"].strip().upper(), []).append(r)
    train, test = [], []
    for label, group in sorted(by_label.items()):
        idx = list(range(len(group)))
        rng.shuffle(idx)
        cut = int(round(len(group) * (1 - test_fraction)))
        train += [group[i] for i in idx[:cut]]
        test += [group[i] for i in idx[cut:]]
    return train, test


def _cell_verdict(test_labels, probs, n_unseen):
    """Apply the pre-committed bar to one (drug, regime) cell's held-out predictions."""
    summary = cal.evaluate_calibration(probs, test_labels)
    cil = summary["calibration_in_the_large"]
    n_pos = sum(1 for y in cal._clean_labels(test_labels) if y == 1)
    n_neg = sum(1 for y in cal._clean_labels(test_labels) if y == 0)
    reasons = [f"test split: {n_pos} R / {n_neg} S; {n_unseen} unseen-stratum excluded"]
    if min(n_pos, n_neg) < BAR["min_per_class"]:
        return "UNDERPOWERED", summary, reasons + [
            f"< {BAR['min_per_class']}/class on the test split -- too thin to read a "
            f"reliability curve. Neither pass nor fail."]
    checks = [
        ("ECE", summary["ece"], BAR["ece"], "<="),
        ("Brier", summary["brier"], BAR["brier"], "<="),
        ("in-large gap", abs((cil["mean_predicted"] or 0) - (cil["observed"] or 0)),
         BAR["in_large_gap"], "<="),
    ]
    passed = True
    for name, val, bar, _op in checks:
        ok = val is not None and val <= bar
        passed = passed and ok
        reasons.append(f"{name}: {val:.3f} (bar <= {bar:.2f}) -- {'ok' if ok else 'BELOW BAR'}")
    return ("PASS" if passed else "FAIL"), summary, reasons


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fixture", default=None, help="the frozen option-C panel TSV (required "
                    "to certify; absent -> BLOCKED)")
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
                 "BLOCKED until the option-C calibration panel exists and its hash is pinned "
                 "(work/PREREGISTRATION_calibration.md). Pass --allow-hash-mismatch to dry-run "
                 "uncertified on a stand-in fixture.")

    rows = _read_fixture(fixture)
    if not rows:
        sys.exit("no fixture rows")

    # group by (drug, regime) -> fit per cell on train, measure on held-out test
    cells = {}
    for r in rows:
        cells.setdefault((r["drug"].strip(), r["regime"].strip().lower()), []).append(r)

    any_fail = False
    for (drug, regime), cell_rows in sorted(cells.items()):
        train, test = _stratified_split(cell_rows, SPLIT_SEED, TEST_FRACTION)
        rates = cal.stratified_rates([r["stratum"] for r in train],
                                     [r["susceptibility"] for r in train])
        probs, n_unseen = cal.apply_rates(rates, [r["stratum"] for r in test])
        test_labels = [r["susceptibility"] for r in test]
        verdict, summary, reasons = _cell_verdict(test_labels, probs, n_unseen)
        any_fail = any_fail or verdict == "FAIL"
        print(f"\n=== {drug} / {regime} regime -- {verdict} ===")
        print(cal.format_report(summary))
        for rr in reasons:
            print("  " + rr)

    print("\n(uncertified dry-run)" if not ok else "")
    return 1 if any_fail else 0


if __name__ == "__main__":
    sys.exit(main())
