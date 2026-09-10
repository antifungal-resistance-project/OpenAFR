"""RUO verdict report -- human- and machine-readable output (issue #138).

Turns the per-drug-class verdict object(s) built by openafr/verdict.py (issue #133) and
produced by the one entrypoint openafr/interpret.py (issue #134) into a single isolate-level
report a lab or a pipeline can consume: a structured JSON form and a human-readable Markdown
form, from the SAME assembled dict.

This module renders; it never re-decides. Every categorical judgement was already made and
frozen by verdict.py/interpret.py; report.py only assembles those verdict objects for one
isolate, attaches a plain-English rationale to each, and lays them out. It invents no new
field the contract does not already define -- above all, NO probability/score/MIC (the
day-one contract is concordance-only, docs/DIAGNOSTIC_VERDICT_CONTRACT.md #132).

Three things this report is built to do (issue #138)
----------------------------------------------------
1. RUO / NO-CLINICAL-CLAIM DISCLAIMER, BAKED IN. Every report -- JSON and Markdown -- carries
   the research-use-only scope and the explicit "not a clinical determination / do not treat
   from this" statement at the top level, not as fine print. The verdict objects already
   carry the per-verdict `scope` string; the report raises it to a standing banner and
   refuses to assemble verdicts whose scope disagrees (a scope drift is a bug, not a report).
2. THE UNCHARACTERIZED VERDICT IS SURFACED, NOT BURIED. A concordance tool's whole moat is
   that it says "I don't know" out loud (#135). So `UNCHARACTERIZED_VARIANT` gets its own
   prominent callout block near the top of the human report and its own top-level list in the
   JSON summary -- never left to be discovered by reading down a table.
3. NO_KNOWN_MARKER IS NEVER RENDERED AS "SUSCEPTIBLE." The rationale for a no-marker result
   spells out what was and was not checked (named ERG11/FKS1 markers only; not efflux/ERG3/
   non-target mechanisms), so the reader cannot mistake it for a clinical "S" (contract rule 2).

Confidence, honestly (issue #138 checklist "... rationale, confidence")
-----------------------------------------------------------------------
The day-one confidence model is CATEGORICAL -- there is no calibrated probability to report
(#132 found no data to calibrate one; deferred to #136). The report states that plainly, and
where an ERG11 verdict carries a structural best-guess (#135) it surfaces that guess's
own calibrated-LOW confidence tier -- clearly flagged as a mechanism inference, never a
resistance call. No numeric confidence is ever emitted.

Stdlib only, pure, offline. Consumes verdict dicts; produces a dict + two strings.
"""
import json

from openafr.verdict import (
    RUO_SCOPE,
    NO_KNOWN_MARKER,
    RESISTANCE_MARKER_DETECTED,
    UNCHARACTERIZED_VARIANT,
    UNRESOLVED,
    VERDICTS,
)

# The standing, top-level disclaimer. Baked into every report so it can never be stripped by
# rendering only part of a verdict object. Says the scope AND the concrete "do not act on this
# clinically" consequence -- an RUO label the reader can act on, not boilerplate.
DISCLAIMER = (
    "RESEARCH USE ONLY (RUO). This report is NOT a clinical determination and must not be "
    "used to guide treatment. It reports concordance with NAMED ERG11/FKS1 resistance "
    "markers only -- it does not measure an MIC, does not emit a calibrated probability, and "
    "does not cover non-target resistance mechanisms (efflux, ERG3, promoter/tandem-repeat). "
    "Absence of a known marker is NOT a susceptibility call."
)

# One-line human headline per verdict enum value -- the "so-what" a reader sees first.
_HEADLINE = {
    RESISTANCE_MARKER_DETECTED: "Known resistance marker detected",
    UNCHARACTERIZED_VARIANT: "Uncharacterized variant -- no known-resistance call",
    NO_KNOWN_MARKER: "No known resistance marker (NOT a susceptibility call)",
    UNRESOLVED: "Unresolved -- could not be called; excluded, not a negative",
}

# Emphasis badge per verdict -- mirrors alert.py's _BADGE idiom.
_BADGE = {
    RESISTANCE_MARKER_DETECTED: "RESISTANCE MARKER",
    UNCHARACTERIZED_VARIANT: "UNCHARACTERIZED",
    NO_KNOWN_MARKER: "NO KNOWN MARKER",
    UNRESOLVED: "UNRESOLVED",
}


def _known_tokens(verdict):
    return [t["token"] for t in verdict["called_tokens"] if t["class"] == "known_resistance"]


def _uncharacterized_tokens(verdict):
    return [t["token"] for t in verdict["called_tokens"] if t["class"] == "uncharacterized"]


def _rationale(verdict):
    """A plain-English rationale for one drug-class verdict -- the 'why', spelled out so a
    non-bioinformatician can read it without decoding the enum. Never asserts susceptibility,
    never implies a probability; states exactly what was and was not checked."""
    v = verdict["verdict"]
    gene = verdict["gene"]
    if v == RESISTANCE_MARKER_DETECTED:
        hits = ", ".join(_known_tokens(verdict))
        note = (f" A co-occurring uncharacterized variant "
                f"({', '.join(_uncharacterized_tokens(verdict))}) is also reported below."
                if _uncharacterized_tokens(verdict) else "")
        return (f"Known {gene} resistance-panel substitution(s) detected: {hits}. This is a "
                f"concordance hit against the named marker panel, not an MIC.{note}")
    if v == UNCHARACTERIZED_VARIANT:
        toks = ", ".join(_uncharacterized_tokens(verdict))
        base = (f"A non-synonymous {gene} change was observed ({toks}) but NONE is a "
                f"known-resistance panel marker. This is the explicit 'uncharacterized' "
                f"verdict -- we do not have a resistance call for it and do not guess one.")
        if verdict.get("structural"):
            base += (" A mechanism-based structural best-guess (calibrated-LOW, see the "
                     "structural block) is attached; it is an inference, not a resistance call.")
        return base
    if v == NO_KNOWN_MARKER:
        return (f"No known-resistance panel marker was found in the resolved {gene} "
                f"resistance windows. This does NOT establish susceptibility: the panel "
                f"covers only named {gene} markers, not efflux/ERG3/non-target mechanisms.")
    # UNRESOLVED
    return (f"The {gene} genotype could not be called honestly "
            f"({verdict.get('resolution_note') or 'unresolved'}). This isolate/drug-class is "
            f"EXCLUDED -- it is not scored as a negative and not read as wild-type.")


def build_report(isolate_id, verdicts, metadata=None):
    """Assemble one isolate's per-drug-class verdicts into a structured report dict.

    isolate_id  the isolate label/accession the verdicts belong to.
    verdicts    an iterable of verdict objects (from openafr.interpret / openafr.verdict), one
                per drug-class typed. Must be non-empty and share a single RUO scope string.
    metadata    optional free-form dict (run id, operator, timestamp) surfaced verbatim.

    Returns a JSON-serializable dict with a standing RUO disclaimer, the verdict objects each
    augmented with a `headline` and plain-English `rationale`, and a `summary` bucketing the
    drug-classes by verdict -- with `uncharacterized` surfaced as its own top-level list so it
    is never buried. Raises ValueError on no verdicts or a disagreeing scope (a scope drift is
    a bug, not something to paper over in a report). Does not mutate the inputs.
    """
    verdicts = list(verdicts)
    if not verdicts:
        raise ValueError("build_report needs at least one verdict object")

    scopes = {v.get("scope") for v in verdicts}
    if scopes != {RUO_SCOPE}:
        raise ValueError(
            f"every verdict must carry the RUO scope {RUO_SCOPE!r}; got {sorted(scopes)!r}"
        )
    for v in verdicts:
        if v.get("verdict") not in VERDICTS:
            raise ValueError(f"unknown verdict value {v.get('verdict')!r}")

    reported = []
    for v in verdicts:
        item = dict(v)                     # copy -- never mutate the caller's verdict object
        item["headline"] = _HEADLINE[v["verdict"]]
        item["rationale"] = _rationale(v)
        reported.append(item)

    def classes_where(pred):
        return [v["drug_class"] for v in verdicts if pred(v)]

    summary = {
        "drug_classes": [v["drug_class"] for v in verdicts],
        "resistance_markers": classes_where(
            lambda v: v["verdict"] == RESISTANCE_MARKER_DETECTED),
        "uncharacterized": classes_where(
            lambda v: v["verdict"] == UNCHARACTERIZED_VARIANT),
        "no_known_marker": classes_where(lambda v: v["verdict"] == NO_KNOWN_MARKER),
        "unresolved": classes_where(lambda v: v["verdict"] == UNRESOLVED),
    }

    return {
        "isolate_id": isolate_id,
        "scope": RUO_SCOPE,
        "disclaimer": DISCLAIMER,
        "confidence_model": (
            "categorical / concordance-only; no calibrated probability day one (#132)"),
        "metadata": dict(metadata or {}),
        "verdicts": reported,
        "summary": summary,
    }


def render_json(report, indent=2):
    """The machine-readable form: the build_report dict as stable JSON text."""
    return json.dumps(report, indent=indent, sort_keys=True)


def _render_structural(structural, lines):
    """Render the ERG11 structural best-guess block (#135) -- clearly a calibrated-LOW
    mechanism inference, mirroring alert.py's structural rendering. FKS1 never reaches here
    (its structural field is null by contract)."""
    lines.append(f"  - **Structural best-guess** [{structural['flag']}, "
                 f"confidence={structural['confidence']}]: {structural['summary']}")
    lines.append(f"    - _basis:_ {structural['basis']}")
    for var in structural["variants"]:
        lines.append(f"    - `{var['token']}` [{var['evidence_class']}, "
                     f"conf={var['confidence']}]: {var['fluconazole_fit_verdict']}")
        if var.get("lost_contacts"):
            lc = ", ".join(f"{c['atom']}@{c['dist_to_drug']:.2f}A"
                           for c in var["lost_contacts"])
            lines.append(f"      - _lost drug contacts (geometry):_ {lc}")


def render_markdown(report):
    """The human form: a one-screen RUO brief a lab can read.

    Standing RUO banner first (never fine print), then -- when present -- a prominent
    UNCHARACTERIZED callout so the honest 'I don't know' is seen before anything else, then
    one compact block per drug-class verdict: headline/badge, the plain-English rationale, the
    called tokens, and (ERG11 only) the calibrated-low structural best-guess.
    """
    lines = [f"# RUO azole/echinocandin resistance-marker report -- isolate `{report['isolate_id']}`",
             "",
             f"> **{report['disclaimer']}**",
             ""]

    meta = report.get("metadata") or {}
    if meta:
        lines.append("- **Run metadata:** "
                     + "; ".join(f"{k}={v}" for k, v in sorted(meta.items())))
    lines.append(f"- **Confidence model:** {report['confidence_model']}")
    lines.append("")

    # Prominent uncharacterized callout -- surfaced up top, never buried (issue #138 / #135).
    unchar = report["summary"]["uncharacterized"]
    if unchar:
        lines += ["> **⚠ Uncharacterized variant(s) present** in: "
                  + ", ".join(unchar)
                  + ". These carry NO known-resistance call -- see the per-class blocks below "
                    "for the explicit verdict and any structural best-guess.",
                  ""]

    for v in report["verdicts"]:
        lines.append(f"## {v['drug_class']} ({v['gene']}) -- {_BADGE[v['verdict']]}")
        lines.append("")
        lines.append(f"**{v['headline']}.** {v['rationale']}")
        lines.append("")
        if v["called_tokens"]:
            toks = ", ".join(f"`{t['token']}` ({t['class']})" for t in v["called_tokens"])
            lines.append(f"- **Called variants:** {toks}")
        else:
            lines.append("- **Called variants:** none")
        if v["resolution"] == "unresolved":
            lines.append(f"- **Resolution:** unresolved "
                         f"({v.get('resolution_note') or 'no reason given'}) -- excluded, "
                         f"not a negative")
        if v.get("structural"):
            _render_structural(v["structural"], lines)
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"
