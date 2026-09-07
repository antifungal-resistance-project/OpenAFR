"""Concordance metrics for a resistance re-caller against an EXTERNAL truth set.

The credibility test for the caller itself (as opposed to the emergence detector on top
of it). The re-caller sanity harnesses (scripts/recaller_sanity{,_fks1}.py) run a *handful*
of published-genotype strains as PASS/FAIL controls; that proves the orchestration is not
grossly broken, but a 3- or 4-strain control is not a measured error rate. This module turns
that control into a graded validation: given a per-isolate list of (expected panel token set,
called panel token set, whether the panel was resolvable, and the isolate's phenotype), it
computes the confusion matrix and the sensitivity/specificity/identity metrics -- each with a
Wilson 95% CI -- against a fixed external benchmark (for FKS1, the 100-isolate PMC12323592 set).

Design decisions, all matching the project's existing discipline:

- **Unresolved isolates are excluded, never counted as negative.** An isolate whose panel
  window was uncalled (coverage gap), refused (in-window indel) or failed to orchestrate is
  `resolved=False`; it is reported as a separate count and drops out of every denominator,
  the same "excluded not guessed" rule as backtest.panel_prevalence and the #22 undated
  discipline. A caller that honestly refuses to guess must not be scored as if it had.

- **Three honest framings, not one headline.** The caller's job is narrow (call the HS1 panel
  tokens); the early-warning question is broad (catch echinocandin resistance). We report:
    1. panel DETECTION -- expected-panel-positive vs any panel call (sensitivity + specificity);
    2. panel IDENTITY -- among detected positives, did the *exact* token match (S639F != S639P);
    3. resistance DETECTION -- of phenotypically resistant isolates, the fraction the HS1 caller
       flags. This is EXPECTED to sit below (1) because non-HS1 mechanisms (F635/D642/R1354/
       W691, non-FKS1) exist that an HS1 caller cannot and should not call. Reporting it is the
       point: it bounds what a detection-only pipeline can catch, honestly.

Pure and offline: it consumes already-evaluated token sets, so it needs no reads, no tools and
no network, and is unit-tested against synthetic evaluations. The orchestration that produces
those evaluations from SRA reads lives in the (tool-gated) validate script, exactly as the
sanity harness reuses recall_fks1.
"""
from openafr.backtest import wilson_interval


def _proportion(k, n):
    """A (k, n, point, ci) record; point/ci are None when n == 0 (nothing to divide)."""
    return {
        "k": k,
        "n": n,
        "point": (k / n if n else None),
        "ci95": wilson_interval(k, n),
    }


def _norm(cell):
    """Normalise a token-set cell into a set of upper-cased tokens; '' / '-' -> empty set."""
    if cell is None:
        return set()
    if isinstance(cell, (set, frozenset, list, tuple)):
        toks = cell
    else:
        s = str(cell).strip()
        toks = [] if s in ("", "-") else s.split(",")
    return {str(t).strip().upper() for t in toks if str(t).strip() not in ("", "-")}


def evaluate(records):
    """Confusion matrix + metrics over per-isolate evaluations against the truth set.

    `records` is an iterable of dicts, one per truth-set isolate:
      expected   token set the truth set attributes to this isolate ('' / '-' / empty = the
                 isolate is panel-negative, i.e. wild-type AT THE PANEL)
      called     panel token set the caller emitted (result['panel_hits'])
      resolved   bool -- False if the panel could not be assessed (uncalled / refused /
                 missing / orchestration failure). Unresolved rows are excluded from every
                 denominator and surfaced as `n_unresolved`.
      resistant  optional bool -- the isolate's echinocandin PHENOTYPE (True=resistant); used
                 only for the resistance-detection framing. None/absent -> omitted from it.

    Returns a dict of metric records (see module docstring for the three framings). Every
    proportion is (k, n, point, Wilson 95% ci).
    """
    n_total = 0
    n_unresolved = 0

    # panel-DETECTION confusion (any expected panel token vs any panel call)
    tp = tn = fp = fn = 0
    # panel-IDENTITY: of detection-TP isolates, exact-token agreement
    identity_match = identity_total = 0
    # per-token sensitivity: token -> [hits, expected-carriers]
    per_token = {}
    # resistance-DETECTION: of resolved *resistant*-phenotype isolates, any panel call
    res_flagged = res_total = 0

    for rec in records:
        n_total += 1
        if not rec.get("resolved", True):
            n_unresolved += 1
            continue
        exp = _norm(rec.get("expected"))
        got = _norm(rec.get("called"))
        exp_pos, got_pos = bool(exp), bool(got)

        if exp_pos and got_pos:
            tp += 1
            identity_total += 1
            if exp <= got:                     # every expected token was emitted
                identity_match += 1
        elif exp_pos and not got_pos:
            fn += 1
        elif not exp_pos and got_pos:
            fp += 1
        else:
            tn += 1

        for tok in exp:
            slot = per_token.setdefault(tok, [0, 0])
            slot[1] += 1
            if tok in got:
                slot[0] += 1

        # resistance-detection is a SENSITIVITY: its denominator is the resistant isolates
        # only (susceptibles belong to panel specificity, not here).
        if rec.get("resistant") is True:
            res_total += 1
            if got_pos:
                res_flagged += 1

    return {
        "n_total": n_total,
        "n_resolved": n_total - n_unresolved,
        "n_unresolved": n_unresolved,
        "confusion": {"tp": tp, "tn": tn, "fp": fp, "fn": fn},
        # panel detection
        "sensitivity": _proportion(tp, tp + fn),
        "specificity": _proportion(tn, tn + fp),
        "accuracy": _proportion(tp + tn, tp + tn + fp + fn),
        # panel identity (exact token, among detected positives)
        "identity": _proportion(identity_match, identity_total),
        # per-token sensitivity
        "per_token": {t: _proportion(k, n) for t, (k, n) in sorted(per_token.items())},
        # resistance detection (phenotype-anchored ceiling on the detection-only pipeline)
        "resistance_detection": _proportion(res_flagged, res_total),
    }


def _fmt(p):
    if p["point"] is None:
        return f"n/a (0 isolates)"
    lo, hi = p["ci95"]
    return f"{p['k']}/{p['n']} = {p['point']:.1%}  [{lo:.1%}, {hi:.1%}]"


def format_report(summary):
    """Human-readable report of an evaluate() result (stderr/stdout friendly, no side effects)."""
    c = summary["confusion"]
    lines = [
        f"isolates: {summary['n_total']} total, {summary['n_resolved']} resolved, "
        f"{summary['n_unresolved']} unresolved (excluded, not counted negative)",
        "",
        f"panel confusion (resolved): TP={c['tp']} TN={c['tn']} FP={c['fp']} FN={c['fn']}",
        f"  sensitivity (panel detection):   {_fmt(summary['sensitivity'])}",
        f"  specificity (wild-type at panel):{_fmt(summary['specificity'])}",
        f"  accuracy:                        {_fmt(summary['accuracy'])}",
        f"  identity (exact token | detected):{_fmt(summary['identity'])}",
    ]
    if summary["per_token"]:
        lines.append("  per-token sensitivity:")
        for tok, p in summary["per_token"].items():
            lines.append(f"    {tok}: {_fmt(p)}")
    rd = summary["resistance_detection"]
    if rd["n"]:
        lines += [
            "",
            f"resistance detection (of resolved RESISTANT-phenotype isolates, the HS1 caller "
            f"flags): {_fmt(rd)}",
            "  (expected BELOW panel sensitivity: non-HS1 mechanisms exist that an HS1 caller "
            "cannot call; this bounds the detection-only pipeline honestly)",
        ]
    return "\n".join(lines)
