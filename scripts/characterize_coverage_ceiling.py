#!/usr/bin/env python3
"""Characterise the echinocandin verdict-accuracy coverage ceiling (#137 follow-up).

The #137 echinocandin arm FAILed its Rung-A bar with a very-major-error (VME) rate of
4/40 = 10.0% (`work/RESULTS_diagnostic_accuracy.md`). That run reported the *rate* but not the
*mechanism* behind each miss. This script answers "which resistance mechanism did each missed /
abstained isolate carry, and is it a tractable FKS1-hotspot the caller could cover?" by joining
the two pinned, on-disk fixtures — no reads, no network, no GCP:

  * data/earlywarning/recaller_sanity/fks1_benchmark_100.tsv  (the PMC12323592 truth set; its
    `note` column carries the paper's per-isolate `paper_mut=<sub> (hsN)` mechanism annotation)
  * data/earlywarning/recaller_sanity/fks1_accuracy_98.tsv     (the harvested #137 fixture; its
    `called_verdict` column is the enum the live caller emitted per isolate)

It reproduces the measured confusion matrix (must match RESULTS_diagnostic_accuracy.md) and then
PROJECTS what an HS3 window (which the caller does not currently read) plus an HS1-panel widening
would recover. The projection is a projection — it is NOT a measured re-claim; that requires the
caller change + a GCP re-run under work/PREREGISTRATION_diagnostic_accuracy_v2.md. Every mechanism
label is copied verbatim from the benchmark `note`, never hand-authored.

Stdlib + openafr.backtest.wilson_interval (the same CI math the measured run uses). Run:
    python3 scripts/characterize_coverage_ceiling.py
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from openafr.backtest import wilson_interval  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SANITY = os.path.join(ROOT, "data", "earlywarning", "recaller_sanity")
BENCH = os.path.join(SANITY, "fks1_benchmark_100.tsv")
ACC = os.path.join(SANITY, "fks1_accuracy_98.tsv")

DET = "RESISTANCE_MARKER_DETECTED"
NKM = "NO_KNOWN_MARKER"
UNC = "UNCHARACTERIZED_VARIANT"

# Windows the caller READS today (openafr/fks1_caller.py FKS1_WINDOWS): HS1 635-643, HS2 1350-1358.
# HS3 (around W691) is NOT read. A residue's hotspot is taken from the benchmark's own `(hsN)` tag.
READ_WINDOWS = {"hs1", "hs2"}


def load_tsv(path):
    rows, header = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            rows.append(dict(zip(header, parts)))
    return rows


def paper_mut(note):
    """Pull `paper_mut=<value>` out of a benchmark note; '' if absent."""
    for tok in (note or "").split(";"):
        tok = tok.strip()
        if tok.startswith("paper_mut="):
            return tok[len("paper_mut="):]
    return ""


def hotspot_of(mut):
    """The (hsN) hotspot the benchmark tags on a mechanism string, lowercased, or 'none'/'?'."""
    m = re.search(r"\(hs(\d)\)", mut)
    if m:
        return "hs" + m.group(1)
    low = mut.strip().lower()
    if low in ("none", "undetermined", ""):
        return "none"
    return "?"


def w(k, n):
    lo, hi = wilson_interval(k, n)
    return f"{k}/{n} = {k / n:.1%}  [{lo:.1%}, {hi:.1%}]" if n else f"{k}/{n} = n/a"


def main():
    bench = {r["strain"]: r for r in load_tsv(BENCH)}
    acc = load_tsv(ACC)

    def mut_of(iso):
        return paper_mut(bench.get(iso, {}).get("note", ""))

    tp = fn = fp = tn = 0
    vme, abst = [], []
    for r in acc:
        ph, v, iso = r["susceptibility"], r["called_verdict"], r["isolate"]
        if v == UNC:
            abst.append(iso)
            continue
        if ph == "R" and v == DET:
            tp += 1
        elif ph == "R" and v == NKM:
            fn += 1
            vme.append(iso)
        elif ph == "S" and v == DET:
            fp += 1
        elif ph == "S" and v == NKM:
            tn += 1

    print("== MEASURED (must match work/RESULTS_diagnostic_accuracy.md) ==")
    print(f"  confusion TP={tp} FN={fn} FP={fp} TN={tn}  abstained={len(abst)}")
    print(f"  VME (missed R): {w(fn, tp + fn)}   bar: point<=3%, Wilson upper<=15%")
    print(f"  ME  (false R):  {w(fp, fp + tn)}")
    print(f"  abstention:     {w(len(abst), len(abst) + tp + fn + fp + tn)}")

    print("\n== VME isolates (scored misses) and mechanism ==")
    print("  isolate   paper_mut                 hotspot  caller reads it?")
    for iso in vme:
        mut = mut_of(iso)
        hs = hotspot_of(mut)
        reads = "yes" if hs in READ_WINDOWS else ("NO (no HS3 window)" if hs == "hs3" else "n/a")
        print(f"  {iso:8}  {mut:24}  {hs:6}   {reads}")

    print("\n== Abstained isolates (UNCHARACTERIZED_VARIANT, excluded from 2x2) and mechanism ==")
    for iso in abst:
        mut = mut_of(iso)
        print(f"  {iso:8}  {mut:24}  {hotspot_of(mut)}")

    # --- PROJECTION (not a measured claim) --------------------------------------------------
    # An HS3 window would let the caller SEE the W691 substitutions it is currently blind to.
    hs3_vme = [iso for iso in vme if hotspot_of(mut_of(iso)) == "hs3" and "W691" in mut_of(iso)]
    residual = [iso for iso in vme if iso not in hs3_vme]
    print("\n== PROJECTION: add an HS3 window tagging W691L/W691C (literature-pinned) ==")
    print(f"  HS3 recovers (VME -> detected): {hs3_vme}")
    print(f"  residual miss: {residual}  ({[mut_of(i) or 'none' for i in residual]})")
    tp2, fn2 = tp + len(hs3_vme), fn - len(hs3_vme)
    print(f"  projected confusion TP={tp2} FN={fn2} FP={fp} TN={tn}")
    print(f"  projected VME: {w(fn2, tp2 + fn2)}   -> clears bar (point<=3%, upper<=15%)")

    print("\n== PROJECTION: widen the HS1 panel (D642Y/F635 already READ but untagged) ==")
    hs1_abst = [iso for iso in abst if hotspot_of(mut_of(iso)) == "hs1"]
    hs2_abst = [iso for iso in abst if hotspot_of(mut_of(iso)) == "hs2"]
    print(f"  HS1 abstentions a widened HS1 panel would tag as detected: {len(hs1_abst)} "
          f"({sorted({mut_of(i) for i in hs1_abst})})")
    print(f"  HS2 abstention an HS2 panel tag would detect: {hs2_abst} "
          f"({[mut_of(i) for i in hs2_abst]})")
    print("  (These lower the abstention rate; they do not change the VME denominator arithmetic.)")

    print("\nNOTE: projections assume the panel mutant set is pinned from literature INDEPENDENT of")
    print("this benchmark (else recovery is circular). The measured re-claim requires the caller")
    print("change + a GCP re-run under work/PREREGISTRATION_diagnostic_accuracy_v2.md.")


if __name__ == "__main__":
    main()
