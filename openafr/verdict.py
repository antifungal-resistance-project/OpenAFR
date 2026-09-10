"""Diagnostic verdict contract (issue #133) -- concordance-only, day one.

Maps a resistance re-caller's output (openafr/recaller.py for ERG11/azole,
openafr/fks1_caller.py for FKS1/echinocandin) into the per-isolate verdict object defined
in docs/DIAGNOSTIC_VERDICT_CONTRACT.md. The engine's whole product is this object and,
above all, the four things it refuses to claim.

Day-one scope is fixed by the #132 panel spike (work/RESULTS_diagnostics_panel_spike.md):
no obtainable public Candida panel can calibrate a resistance PROBABILITY, so the verdict
is CONCORDANCE-ONLY -- it detects *named resistance markers*, categorically. There is no
probability field; that is deferred to the pooled-C. albicans calibration track (#136) and
re-enters this contract only after that track produces a validated calibration.

Pure and offline: it consumes the callers' already-computed output dicts (the same
`tokens` / `panel_hits` / `uncalled_panel` keys both callers expose), so it needs no reads,
no tools and no network, and is unit-tested against synthetic caller outputs.

The four honesty rules this module enforces (docs/DIAGNOSTIC_VERDICT_CONTRACT.md)
------------------------------------------------------------------------------------
1. NO PROBABILITY FIELD. The verdict is one of four categorical values; #132 found no data
   to calibrate a probability, so emitting one would be an unbacked number.
2. `NO_KNOWN_MARKER` IS NOT A SUSCEPTIBILITY CALL. We detect known ERG11/FKS1 markers; we do
   not cover efflux (TAC1/MRR1/CDR1), ERG3 or non-target mechanisms, so absence of a marker
   cannot assert susceptibility. The RUO scope string forbids reading it as a clinical "S".
3. `UNCHARACTERIZED_VARIANT` IS A FIRST-CLASS VERDICT, NOT SILENCE. A non-synonymous change
   that is not a known marker gets an explicit verdict, never a dropped/blank result (#135).
4. AN UNCALLED PANEL RESIDUE IS UNRESOLVED, NEVER `NO_KNOWN_MARKER`. If a known-marker
   position could not be read, we cannot rule resistance out there, so the drug-class
   verdict is UNRESOLVED -- the callers' "an uncalled codon is uncalled, never wild-type"
   rule carried up to the verdict level. UNRESOLVED is excluded downstream, never a negative.
"""

RUO_SCOPE = "RUO -- research use only; not a clinical determination"

# The four-value categorical enum. No fifth value, and never a probability.
RESISTANCE_MARKER_DETECTED = "RESISTANCE_MARKER_DETECTED"
UNCHARACTERIZED_VARIANT = "UNCHARACTERIZED_VARIANT"
NO_KNOWN_MARKER = "NO_KNOWN_MARKER"
UNRESOLVED = "UNRESOLVED"

VERDICTS = (
    RESISTANCE_MARKER_DETECTED,
    UNCHARACTERIZED_VARIANT,
    NO_KNOWN_MARKER,
    UNRESOLVED,
)

# gene -> (drug_class, carries a structural pocket verdict?). FKS1 is detection-only:
# echinocandins coordinate no metal, so the CYP51/heme-iron geometry moat does not transfer
# (openafr/fks1_caller.py), and its structural field is ALWAYS null.
_GENE = {
    "ERG11": ("azole", True),
    "FKS1": ("echinocandin", False),
}


def _classify(panel_hits, tokens, uncalled_panel):
    """The four-value decision, in precedence order.

    A positive marker wins even if other positions were uncalled (a detected S639F is a
    detected S639F). Otherwise an uncalled *panel* residue forces UNRESOLVED -- we cannot
    complete the marker check, so we must not emit NO_KNOWN_MARKER (rule 4). Only when every
    panel residue was confidently read do we distinguish an uncharacterized non-synonymous
    change (rule 3) from a clean no-marker result (rule 2).
    """
    if panel_hits:
        return RESISTANCE_MARKER_DETECTED
    if uncalled_panel:
        return UNRESOLVED
    if tokens:
        return UNCHARACTERIZED_VARIANT
    return NO_KNOWN_MARKER


def build_verdict(gene, call=None, unresolved_reason=None, structural=None,
                  provenance=None):
    """Build one verdict object for a single isolate and a single drug-class gene.

    gene               'ERG11' or 'FKS1'.
    call               a caller output dict with keys `tokens`, `panel_hits` and
                       `uncalled_panel` (recaller.call_substitutions /
                       fks1_caller.call_windows). Omit when the isolate is unresolved.
    unresolved_reason  set this (e.g. the str of a ConsensusError, or 'coverage_gap')
                       when the caller could not place the consensus at all; produces a
                       hard UNRESOLVED verdict. Takes precedence over `call`.
    structural         optional ERG11 pocket verdict (#135); forced to None for FKS1.
    provenance         optional dict (reference accession/sha, caller/panel version).

    Never raises on a normal caller output; the caller is expected to translate its own
    hard failures (frameshift/length ConsensusError, orchestration failure) into
    `unresolved_reason` rather than letting them escape here.
    """
    if gene not in _GENE:
        raise ValueError(f"unknown gene {gene!r}; expected one of {sorted(_GENE)}")
    drug_class, carries_structural = _GENE[gene]

    def obj(verdict, resolution, note, called_tokens, structural_out):
        return {
            "drug_class": drug_class,
            "gene": gene,
            "resolution": resolution,
            "resolution_note": note,
            "called_tokens": called_tokens,
            "verdict": verdict,
            "structural": structural_out,
            "provenance": dict(provenance or {}),
            "scope": RUO_SCOPE,
        }

    # Hard unresolved: the caller could not produce a call at all.
    if unresolved_reason is not None:
        return obj(UNRESOLVED, "unresolved", str(unresolved_reason), [], None)

    if call is None:
        raise ValueError("build_verdict needs either `call` or `unresolved_reason`")

    panel_hits = list(call.get("panel_hits") or [])
    tokens = list(call.get("tokens") or [])
    uncalled_panel = list(call.get("uncalled_panel") or [])
    panel_set = set(panel_hits)

    verdict = _classify(panel_hits, tokens, uncalled_panel)

    # called_tokens surfaces every non-synonymous change verbatim, each tagged known vs
    # uncharacterized -- populated even when the overall verdict is UNRESOLVED, so a reader
    # sees what WAS observed alongside the honest "a panel position could not be read".
    called_tokens = [
        {"token": t,
         "class": "known_resistance" if t in panel_set else "uncharacterized"}
        for t in tokens
    ]

    if verdict == UNRESOLVED:
        pos = "+".join(str(p) for p in uncalled_panel)
        note = f"panel_residue_uncalled:{pos}"
        resolution = "unresolved"
    else:
        note = None
        resolution = "resolved"

    structural_out = structural if carries_structural else None
    return obj(verdict, resolution, note, called_tokens, structural_out)


def azole_verdict(call=None, unresolved_reason=None, structural=None, provenance=None):
    """ERG11 -> azole verdict (carries a structural pocket verdict where mappable, #135)."""
    return build_verdict("ERG11", call=call, unresolved_reason=unresolved_reason,
                         structural=structural, provenance=provenance)


def echinocandin_verdict(call=None, unresolved_reason=None, provenance=None):
    """FKS1 -> echinocandin verdict (detection-only; structural is always null)."""
    return build_verdict("FKS1", call=call, unresolved_reason=unresolved_reason,
                         provenance=provenance)
