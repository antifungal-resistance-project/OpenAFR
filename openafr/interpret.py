"""One genotype -> verdict entrypoint (issue #134).

The single interpretation path for the resistance-interpretation engine
(.context/DIAGNOSTICS_CONCEPT.md). It wraps the existing re-callers behind ONE call that
turns raw isolate input -- a consensus sequence OR a variant list, from any WGS / targeted
panel -- into the concordance-only verdict object defined in
docs/DIAGNOSTIC_VERDICT_CONTRACT.md and built by openafr/verdict.py (issue #133).

It REUSES, never forks, the resistance early-warning callers:
  * ERG11 / azole        -> openafr/recaller.py  (call_substitutions)
  * FKS1  / echinocandin -> openafr/fks1_caller.py (call_windows / slice_windows_from_cds)
and, for an ERG11 azole verdict, populates the structural pocket best-guess (issue #135):
an uncharacterized (non-panel) variant is mapped into the modeled azole pocket via
openafr/structural.py and gets a calibrated-LOW, uncharacterized-flagged mechanism guess --
so this engine does not fall silent on a novel variant the way a lookup table does -- while
a caller-supplied structural verdict still passes straight through. Nothing here
re-implements a caller; it only (a) picks the right caller for the gene, (b) accepts the
handful of raw input shapes a real WGS/panel pipeline hands us, (c) translates the
callers' outputs -- and their honest hard-failures -- into the #133 schema, and (d) fills the
structural best-guess for an uncharacterized ERG11 variant.

What "wiring" means here, precisely
-----------------------------------
1. RAW INPUT, NOT PRE-COMPUTED CALLS. verdict.py already turns a caller's output dict into a
   verdict; this module is the step before that -- it runs the caller on a sequence/variant
   list first. That is the whole point of #134: expose a path a clinician's pipeline can
   actually feed, not one that assumes the call is already done.
2. HARD FAILURES BECOME UNRESOLVED, NOT EXCEPTIONS. A ConsensusError (ERG11 frameshift/length)
   or WindowError (FKS1 un-anchorable consensus) is an honest "this isolate could not be
   called", so it is turned into an UNRESOLVED verdict carrying the reason -- never leaked as
   a crash. A ReferenceError is NOT caught: a missing/mismatched pinned reference is an
   environment fault, not a property of the isolate, and must surface loudly.
3. FKS1 MISSING/REFUSED WINDOWS FORCE UNRESOLVED AT THEIR PANEL RESIDUES. call_windows drops
   the residues of a 'missing' or 'refused' window from uncalled_panel (a refused window
   reports no uncalled positions). Left as-is that would let an unread S639 slip out as
   NO_KNOWN_MARKER -- exactly the "an uncalled panel residue is never wild-type" rule the
   contract forbids breaking. So we re-derive uncalled_panel from the per-window statuses
   before building the verdict (normalization, not a fork of the caller).

Honesty caveat on the variant-list input
-----------------------------------------
A bare variant list is a POSITIVE-call format: it says what was found, not what was covered.
We take it as having covered the resistance-panel positions (a variant caller that reports
substitutions sequenced the region), so an absent panel token reads as confidently wild-type
there. When that assumption does not hold -- a targeted panel with a known coverage gap over
a marker -- pass `uncalled=[...]` with those residue positions and they force UNRESOLVED at
the panel, exactly as an uncalled codon does on the sequence path.

Stdlib only, offline. Delegates every base/codon decision to the pinned-reference callers and
every pocket/structural decision to openafr/structural.py (checked-in modeled receptor).
"""
import re

from openafr import fks1_caller, recaller
from openafr import structural as _structural
from openafr.verdict import azole_verdict, echinocandin_verdict

# A canonical protein point-substitution token, e.g. 'Y132F': wild-type, position, mutant.
# Same shape both callers emit and tag (recaller._TOKEN_RE / fks1_caller._TOKEN_RE).
_TOKEN_RE = re.compile(r"^([A-Za-z])(\d+)([A-Za-z*])$")

# Genes this engine can actually call today. CYP51A (the Aspergillus azole target) is named
# in the roadmap but has NO caller yet -- flipping to it pays a documented build cost
# (work/PREREGISTRATION_diagnostics_panel.md), so we refuse it explicitly rather than
# pretend to interpret it.
SUPPORTED_GENES = ("ERG11", "FKS1")


def _pos(token):
    return int(re.search(r"\d+", token).group())


def _normalize_variants(variants, is_panel_token):
    """Turn a raw variant-list input into a caller-shaped (tokens, panel_hits) pair.

    Validates each token's shape (rejecting unparseable ones loudly rather than silently
    dropping a mis-typed marker), discards no-change tokens (wt == mut), de-duplicates and
    sorts by residue, and tags panel hits with the SAME test the caller uses -- so a variant
    list and a sequence produce identical tagging for the same substitutions.
    """
    tokens, bad = [], []
    for raw in variants:
        m = _TOKEN_RE.match(str(raw).strip())
        if not m:
            bad.append(raw)
            continue
        wt, pos, mut = m.group(1).upper(), m.group(2), m.group(3).upper()
        if wt == mut:
            continue                       # not a substitution -- no change at this residue
        tokens.append(f"{wt}{pos}{mut}")
    if bad:
        raise ValueError(
            f"unparseable variant token(s): {bad!r}; expected protein tokens like 'Y132F'"
        )
    tokens = sorted(set(tokens), key=_pos)
    panel_hits = [t for t in tokens if is_panel_token(t)]
    return tokens, panel_hits


def _fks1_uncalled_panel(call):
    """Re-derive FKS1 uncalled_panel so a missing/refused hot-spot window can't leak out as
    NO_KNOWN_MARKER. call_windows already reports the panel residues uncalled *within a
    partial window*; here we additionally treat every panel residue of a window that was
    'missing' (no data) or 'refused' (in-window indel) as uncalled -- we could not read it,
    so we cannot rule its marker out (contract rule 4)."""
    uncalled = set(call["uncalled_panel"])
    for name, win_call in call["windows"].items():
        if win_call["status"] in ("missing", "refused"):
            uncalled.update(fks1_caller.FKS1_WINDOWS[name].panel)
    return sorted(uncalled)


def _provenance(module, extra):
    """Auto-stamp the provenance fields we can assert from code (the pinned reference hash
    and the exact panel positions in play), then layer any caller-supplied fields on top.
    Everything auto-filled is computed from the module, so it can never drift from what
    actually made the call."""
    prov = {
        "reference_sha256": module.REFERENCE_CDS_SHA256,
        "known_panel": sorted(module.RESISTANCE_PANEL),
    }
    if extra:
        prov.update(extra)
    return prov


def _one_input(cds, windows, variants):
    supplied = [name for name, val in
                (("cds", cds), ("windows", windows), ("variants", variants))
                if val is not None]
    if len(supplied) != 1:
        raise ValueError(
            "interpret() needs exactly one of cds=, windows=, variants=; "
            f"got {supplied or 'none'}"
        )
    return supplied[0]


def _erg11_structural(call, structural):
    """Populate the ERG11 structural field with a mechanism-based best-guess for any
    uncharacterized (non-panel) variant the isolate carries (#135) -- the moat: a lookup-table
    tool falls silent on a novel variant, we map it into the pocket and emit a calibrated-low
    structural best-guess instead. Left to `structural` when the caller supplied one
    explicitly; None (honest null) when there is nothing uncharacterized to place."""
    if structural is not None:
        return structural
    panel = set(call.get("panel_hits") or [])
    uncharacterized = [t for t in (call.get("tokens") or []) if t not in panel]
    if not uncharacterized:
        return None
    return _structural.diagnostic_best_guess(uncharacterized)


def _interpret_erg11(kind, cds, variants, uncalled, structural, prov, reference):
    if kind == "windows":
        raise ValueError("windows= is FKS1-only; ERG11 takes a full-length cds= or variants=")

    if kind == "cds":
        try:
            call = recaller.call_substitutions(cds, reference=reference)
        except recaller.ConsensusError as e:
            # Honest hard-failure (frameshift/length): UNRESOLVED, reason carried, not a crash.
            return azole_verdict(unresolved_reason=str(e), provenance=prov)
        return azole_verdict(call, structural=_erg11_structural(call, structural),
                             provenance=prov)

    # variants: build the caller-shaped dict, honouring an explicit panel coverage gap.
    tokens, panel_hits = _normalize_variants(variants, recaller.is_panel_token)
    uncalled_panel = sorted(p for p in (uncalled or []) if p in recaller.RESISTANCE_PANEL)
    call = {"tokens": tokens, "panel_hits": panel_hits, "uncalled_panel": uncalled_panel}
    return azole_verdict(call, structural=_erg11_structural(call, structural),
                         provenance=prov)


def _interpret_fks1(kind, cds, windows, variants, uncalled, prov, reference):
    if kind == "cds":
        # The fragile full-CDS bridge: an un-anchorable consensus is an honest UNRESOLVED.
        try:
            windows = fks1_caller.slice_windows_from_cds(cds)
        except fks1_caller.WindowError as e:
            return echinocandin_verdict(unresolved_reason=str(e), provenance=prov)
        kind = "windows"

    if kind == "windows":
        call = fks1_caller.call_windows(windows, reference=reference)
        norm = {"tokens": call["tokens"], "panel_hits": call["panel_hits"],
                "uncalled_panel": _fks1_uncalled_panel(call)}
        return echinocandin_verdict(norm, provenance=prov)

    # variants
    tokens, panel_hits = _normalize_variants(variants, fks1_caller.is_panel_token)
    uncalled_panel = sorted(p for p in (uncalled or []) if p in fks1_caller.RESISTANCE_PANEL)
    call = {"tokens": tokens, "panel_hits": panel_hits, "uncalled_panel": uncalled_panel}
    return echinocandin_verdict(call, provenance=prov)


def interpret(gene, *, cds=None, windows=None, variants=None, uncalled=None,
              structural=None, provenance=None, reference=None):
    """Interpret one gene's genotype for one isolate -> a #133 verdict object.

    gene        'ERG11' (azole) or 'FKS1' (echinocandin), case-insensitive. 'CYP51A' is
                recognised but refused -- no caller exists for it yet.
    Exactly one raw input:
      cds       consensus nucleotide CDS. ERG11: full-length, in-frame. FKS1: full-length
                CDS (the fragile length-exact bridge; prefer windows= when available).
      windows   FKS1 only -- {window_name: consensus_nt} for the hot-spots (HS1/HS2).
      variants  an iterable of protein substitution tokens (e.g. ['Y132F', 'T123I']) from a
                variant caller / targeted panel.
    uncalled    optional residue positions NOT covered; only meaningful with variants=,
                where panel residues among them force UNRESOLVED (declares a coverage gap a
                bare variant list cannot).
    structural  optional ERG11 pocket verdict override; when omitted, an uncharacterized
                ERG11 variant auto-gets a calibrated-low structural best-guess (#135 /
                openafr/structural.py). Dropped for FKS1 by the verdict contract.
    provenance  optional extra provenance fields, layered over the auto-stamped reference
                hash + panel positions.
    reference   optional pre-loaded (cds, protein) tuple handed to the caller (reuse/tests).

    A full isolate typed for both genes is TWO calls (one per drug-class), matching the
    contract's one-verdict-per-drug-class shape. Returns the verdict dict; raises ValueError
    only on a caller-agnostic misuse (unknown gene, wrong/duplicate input, unparseable
    token). A ConsensusError/WindowError becomes an UNRESOLVED verdict; a ReferenceError
    (broken pinned reference) is left to propagate.
    """
    g = (gene or "").strip().upper()
    if g == "CYP51A":
        raise ValueError(
            "CYP51A (Aspergillus azole target) has no caller yet -- flipping to it pays a "
            "build cost (work/PREREGISTRATION_diagnostics_panel.md); only ERG11 and FKS1 "
            "are wired."
        )
    if g not in SUPPORTED_GENES:
        raise ValueError(f"unknown gene {gene!r}; expected one of {list(SUPPORTED_GENES)}")

    kind = _one_input(cds, windows, variants)

    if g == "ERG11":
        prov = _provenance(recaller, provenance)
        return _interpret_erg11(kind, cds, variants, uncalled, structural, prov, reference)

    prov = _provenance(fks1_caller, provenance)
    return _interpret_fks1(kind, cds, windows, variants, uncalled, prov, reference)
