"""Verdict-level accuracy metrics against phenotypic AST truth (issue #137).

The credibility test for the *diagnostics MVP product* -- the verdict contract
(openafr/verdict.py, docs/DIAGNOSTIC_VERDICT_CONTRACT.md), as opposed to the caller
underneath it. openafr/concordance.py already grades the caller's TOKEN accuracy against a
genotype truth set ("did it emit S639F when S639F was present?"). This module asks the
next, clinically loaded question: when the engine emits a per-drug-class VERDICT, how often
does it agree with the isolate's phenotype, and -- above all -- how often does it get the
two *dangerous* directions wrong:

  * VERY MAJOR ERROR (VME): phenotype R, verdict emits no marker (NO_KNOWN_MARKER). The
    assay missed real resistance. This is the load-bearing error for a resistance test: it
    is precisely the failure the contract's rule 2 ("NO_KNOWN_MARKER is NOT susceptible")
    warns about, and this module MEASURES its rate instead of asserting it is rare.
  * MAJOR ERROR (ME): phenotype S, verdict RESISTANCE_MARKER_DETECTED. A false resistance
    call (the marker is real but the isolate is phenotypically susceptible -- e.g. a
    low-penetrance ERG11 change).

These are the CLSI-M23 error names an AST-device evaluation reports, denominated the CLSI
way: VME over the RESISTANT isolates, ME over the SUSCEPTIBLE ones. Sensitivity = 1 - VME,
specificity = 1 - ME; we report both framings because the error names are what a clinician
reads and the sens/spec framing is what a caller person reads.

The verdict enum -> R/S CALL mapping (fixed by the pre-registration, work/
PREREGISTRATION_diagnostic_accuracy.md)
--------------------------------------------------------------------------------------
  RESISTANCE_MARKER_DETECTED  -> positive call  ("R")
  NO_KNOWN_MARKER             -> negative call  ("not-R")  [NOT a susceptibility CLAIM;
                                 it is the assay's negative prediction, and the VME rate is
                                 exactly the honest price of reading it as one]
  UNCHARACTERIZED_VARIANT     -> ABSTENTION -- excluded from the 2x2, reported as an
                                 abstention rate. The contract's first-class "I don't know"
                                 (rule 3) is never scored as either a right or a wrong call.
  UNRESOLVED                  -> EXCLUDED -- reported as n_unresolved, never a negative
                                 (rule 4 / the "excluded, not guessed" discipline shared with
                                 concordance.evaluate and backtest.panel_prevalence).

An isolate with no phenotype truth is unscorable and counted as n_no_phenotype (excluded).
Abstentions and unresolved verdicts are held OUT of the confusion matrix, so a low VME can
never be bought by quietly reclassifying hard isolates as negatives.

Pure and offline: it consumes already-built verdict objects (or bare enum values) paired
with a phenotype label, so it needs no reads, no tools and no network, and is unit-tested
against synthetic verdicts. The orchestration that turns SRA reads into verdicts lives in
the (tool-gated, data-blocked) grader scripts/validate_diagnostic_accuracy.py.
"""
from openafr.backtest import wilson_interval
from openafr.verdict import (
    NO_KNOWN_MARKER,
    RESISTANCE_MARKER_DETECTED,
    UNCHARACTERIZED_VARIANT,
    UNRESOLVED,
)

# Verdict -> scoring disposition. "R"/"not_R" enter the 2x2; "abstain"/"exclude" do not.
_DISPOSITION = {
    RESISTANCE_MARKER_DETECTED: "R",
    NO_KNOWN_MARKER: "not_R",
    UNCHARACTERIZED_VARIANT: "abstain",
    UNRESOLVED: "exclude",
}


def _proportion(k, n):
    """A (k, n, point, ci) record; point/ci are None when n == 0 (nothing to divide)."""
    return {
        "k": k,
        "n": n,
        "point": (k / n if n else None),
        "ci95": wilson_interval(k, n),
    }


def _phenotype(value):
    """Normalise a phenotype label to 'R', 'S' or None (unscorable).

    Accepts 'R'/'S' (any case), True/False (True=resistant), or None/''/'-' -> None.
    """
    if value is None:
        return None
    if value is True:
        return "R"
    if value is False:
        return "S"
    s = str(value).strip().upper()
    if s in ("R", "RESISTANT"):
        return "R"
    if s in ("S", "SUSCEPTIBLE"):
        return "S"
    if s in ("", "-"):
        return None
    raise ValueError(f"unrecognised phenotype {value!r}; expected R/S, True/False or empty")


def _verdict_of(rec):
    """Pull the enum value out of a record that carries either a bare enum or a verdict obj."""
    v = rec.get("verdict")
    if isinstance(v, dict):          # a full verdict object from openafr.verdict.build_verdict
        return v.get("verdict")
    return v


def _drug_class_of(rec):
    """Prefer an explicit drug_class; else read it off an embedded verdict object."""
    if rec.get("drug_class"):
        return rec["drug_class"]
    v = rec.get("verdict")
    if isinstance(v, dict):
        return v.get("drug_class")
    return None


def record(verdict_obj, phenotype):
    """Build one accuracy record from a build_verdict() object + a phenotype label."""
    return {
        "drug_class": verdict_obj.get("drug_class"),
        "verdict": verdict_obj.get("verdict"),
        "phenotype": phenotype,
    }


def _evaluate_one_class(records):
    """Confusion matrix + VME/ME/sens/spec for a single drug class's records."""
    n_total = n_unresolved = n_no_phenotype = n_abstained = 0
    tp = fn = fp = tn = 0                       # R/not_R call x R/S phenotype

    for rec in records:
        n_total += 1
        verdict = _verdict_of(rec)
        disp = _DISPOSITION.get(verdict)
        if disp is None:
            raise ValueError(f"unrecognised verdict {verdict!r}")
        pheno = _phenotype(rec.get("phenotype"))

        if disp == "exclude":                  # UNRESOLVED -- never a negative
            n_unresolved += 1
            continue
        if pheno is None:                      # no AST truth -> unscorable
            n_no_phenotype += 1
            continue
        if disp == "abstain":                  # UNCHARACTERIZED -- held out of the 2x2
            n_abstained += 1
            continue

        # disp is "R" or "not_R" and we have a phenotype: this row enters the confusion matrix
        if pheno == "R" and disp == "R":
            tp += 1
        elif pheno == "R" and disp == "not_R":
            fn += 1                             # VERY MAJOR: missed real resistance
        elif pheno == "S" and disp == "R":
            fp += 1                             # MAJOR: false resistance call
        else:
            tn += 1

    n_scored = tp + fn + fp + tn
    return {
        "n_total": n_total,
        "n_unresolved": n_unresolved,
        "n_no_phenotype": n_no_phenotype,
        "n_abstained": n_abstained,
        "n_scored": n_scored,
        "confusion": {"tp": tp, "fn": fn, "fp": fp, "tn": tn},
        "sensitivity": _proportion(tp, tp + fn),
        "specificity": _proportion(tn, tn + fp),
        "categorical_agreement": _proportion(tp + tn, n_scored),
        # CLSI-named errors: VME over the resistant isolates, ME over the susceptible.
        "very_major_error": _proportion(fn, tp + fn),
        "major_error": _proportion(fp, fp + tn),
        # abstention rate over the isolates that HAD a phenotype and a categorical verdict
        # (scored + abstained); an honest denominator that excludes the unresolved/unphenotyped.
        "abstention": _proportion(n_abstained, n_abstained + n_scored),
    }


def evaluate(records):
    """Per-drug-class verdict accuracy over an iterable of records.

    Each record is a dict with:
      verdict     one of openafr.verdict.VERDICTS, OR a full verdict object (dict) from
                  build_verdict() -- in which case drug_class is read off it too.
      drug_class  'azole' | 'echinocandin' (optional if `verdict` is a verdict object).
      phenotype   the isolate's AST truth: 'R'/'S', True/False, or None/'' (unscorable).

    Returns {"by_drug_class": {cls: <metrics>}, "n_total": int, "drug_classes": [...]}.
    Metrics are grouped per drug class because VME/ME are per-drug quantities (an azole VME
    and an echinocandin VME are different clinical facts and are never pooled).
    """
    grouped = {}
    n_total = 0
    for rec in records:
        n_total += 1
        cls = _drug_class_of(rec)
        if not cls:
            raise ValueError("record has no drug_class (pass one, or a full verdict object)")
        grouped.setdefault(cls, []).append(rec)

    return {
        "n_total": n_total,
        "drug_classes": sorted(grouped),
        "by_drug_class": {cls: _evaluate_one_class(recs)
                          for cls, recs in sorted(grouped.items())},
    }


def apply_bar(class_summary, bar):
    """Grade one drug class's metrics against a pre-committed pass bar -> verdict + reasons.

    `bar` keys (all pre-registered in work/PREREGISTRATION_diagnostic_accuracy.md):
      vme_point / vme_upper   max very-major-error point estimate / Wilson UPPER bound
      me_point                max major-error point estimate
      min_resistant / min_susceptible
                              minimum scored isolates per phenotype arm to read an error rate
      max_abstention          max fraction of characterised isolates the engine may abstain on

    Returns one of 'PASS', 'FAIL', 'UNDERPOWERED' plus a list of human-readable reasons. The
    error bounds are gated on the Wilson UPPER bound (not just the point) because with small n
    a 0/20 VME still admits a double-digit true rate, and the pre-commitment must be honest
    about the power the truth set actually supplies.
    """
    c = class_summary["confusion"]
    n_res = c["tp"] + c["fn"]
    n_sus = c["fp"] + c["tn"]
    reasons = [f"scored {class_summary['n_scored']} isolates: {n_res} R / {n_sus} S; "
               f"{class_summary['n_abstained']} abstained, "
               f"{class_summary['n_unresolved']} unresolved, "
               f"{class_summary['n_no_phenotype']} no-phenotype (all excluded)"]

    if n_res < bar["min_resistant"] or n_sus < bar["min_susceptible"]:
        reasons.append(f"< {bar['min_resistant']} R or < {bar['min_susceptible']} S scored -- "
                       "too thin to read VME/ME. Neither pass nor fail.")
        return "UNDERPOWERED", reasons

    vme = class_summary["very_major_error"]
    me = class_summary["major_error"]
    abst = class_summary["abstention"]
    vme_upper = vme["ci95"][1] if vme["ci95"] else 1.0
    checks = [
        ("VME point", vme["point"], bar["vme_point"]),
        ("VME Wilson upper", vme_upper, bar["vme_upper"]),
        ("ME point", me["point"], bar["me_point"]),
        ("abstention", abst["point"], bar["max_abstention"]),
    ]
    passed = True
    for name, val, lim in checks:
        ok = val is not None and val <= lim
        passed = passed and ok
        shown = "n/a" if val is None else f"{val:.1%}"
        reasons.append(f"{name}: {shown} (bar <= {lim:.0%}) -- {'ok' if ok else 'BELOW BAR'}")
    return ("PASS" if passed else "FAIL"), reasons


def _fmt(p):
    if p["point"] is None:
        return "n/a (0 isolates)"
    lo, hi = p["ci95"]
    return f"{p['k']}/{p['n']} = {p['point']:.1%}  [{lo:.1%}, {hi:.1%}]"


def format_report(summary):
    """Human-readable report of an evaluate() result (no side effects)."""
    lines = [f"isolates: {summary['n_total']} total across "
             f"{len(summary['drug_classes'])} drug class(es)"]
    for cls in summary["drug_classes"]:
        s = summary["by_drug_class"][cls]
        c = s["confusion"]
        lines += [
            "",
            f"== {cls} ==",
            f"  scored {s['n_scored']}  (excluded: {s['n_abstained']} abstained, "
            f"{s['n_unresolved']} unresolved, {s['n_no_phenotype']} no-phenotype)",
            f"  confusion (R/S x call): TP={c['tp']} FN={c['fn']} FP={c['fp']} TN={c['tn']}",
            f"  very major error (missed R):  {_fmt(s['very_major_error'])}",
            f"  major error (false R):        {_fmt(s['major_error'])}",
            f"  sensitivity (= 1 - VME):      {_fmt(s['sensitivity'])}",
            f"  specificity (= 1 - ME):       {_fmt(s['specificity'])}",
            f"  categorical agreement:        {_fmt(s['categorical_agreement'])}",
            f"  abstention (uncharacterised): {_fmt(s['abstention'])}",
        ]
    return "\n".join(lines)
