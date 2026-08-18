"""ERG11 re-caller -- reads/consensus -> azole-resistance substitution call (issue: TODOS top-priority).

Pipeline stage: the single missing piece that turns the early-warning scaffold into
a working early warning. The snapshot store (openafr.earlywarning) leaves `erg11_call`
empty and records `resistance_source = pending:sra-recaller`, because the spike (#20)
proved NCBI runs no AMR pipeline on *C. auris* -- there is no resistance call to read
off the feed. It has to be *manufactured* from the isolate's own ERG11 sequence. This
module is the deterministic, testable heart of that: given an isolate's ERG11 coding
consensus, it emits the substitution tokens (`Y132F,K143R`) that every downstream stage
(emergence -> mapping -> structural so-what -> alert) already consumes.

The boundary this module draws
------------------------------
Reads -> a per-isolate ERG11 *consensus CDS* is the messy, tool-and-network half
(fasterq-dump / minimap2 / bcftools; see scripts/recall_erg11.py, gated on those
binaries exactly like the docking scripts gate on vina). This module starts *after*
that: a consensus CDS in, substitution tokens out. That boundary is deliberate --
everything on this side is stdlib-only and deterministic, so it is unit-tested against
the pinned reference; the tool-orchestration side, which cannot be reproduced without
network + large downloads, is kept out of the tested core and honestly labelled.

The honesty constraints (mirrors emergence.py / mapping.py)
-----------------------------------------------------------
1. WE ASSERT THE REFERENCE, WE DO NOT ASSUME IT. Every call is relative to the pinned
   B8441 ERG11 reference (data/earlywarning/erg11_reference/, accession XM_085597798.1,
   sha256 verified at load). load_reference() refuses a reference whose hotspot residues
   (V125/F126/Y132/K143) do not match the panel's wild-type letters -- a numbering or
   sequence swap is surfaced, never silently mis-called.

2. AN UNCALLED CODON IS UNCALLED, NEVER "WILD TYPE". A codon containing an ambiguous or
   gap base (N, -, IUPAC ambiguity) cannot be translated to a confident residue, so that
   position is reported as *uncalled* -- it never contributes a wild-type "no change".
   An isolate that could not be sequenced across a hotspot says so; it is not quietly
   assumed sensitive.

3. FRAME IS ASSERTED, NOT GUESSED. The consensus must be an in-frame CDS the same length
   as the reference (a length that is not a multiple of 3, or that disagrees with the
   reference length, is a frameshift/indel we cannot represent as point substitutions).
   Such input is refused with the reason, never forced into a frame.

4. WE TAG THE PANEL, WE DO NOT INVENT NOVELTY. Every non-synonymous difference from the
   reference is emitted as a token. Known azole-resistance panel positions are *tagged*
   so the alert can say "this is a named resistance mutation", but non-panel substitutions
   (clade-defining polymorphisms, VUS) are emitted verbatim, never suppressed and never
   asserted to be resistance. Deciding which emitted token is novel/rising is the
   emergence detector's job, not the re-caller's.

Stdlib only. Operates on the reference FASTA and on a consensus CDS string.
"""
import hashlib
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFERENCE_DIR = os.path.join(ROOT, "data", "earlywarning", "erg11_reference")
DEFAULT_REFERENCE = os.path.join(REFERENCE_DIR, "erg11_cds.fasta")

# The pinned reference CDS content hash (bare ACGT string, no header/newlines).
# Must match data/earlywarning/erg11_reference/PROVENANCE.md. A reference whose bytes
# drift from this fails load unless explicitly allowed -- the same freeze contract as
# protocol.yaml and the snapshot INDEX hashes.
REFERENCE_CDS_SHA256 = (
    "764e6d1a429dcfda48e8904d8b438425aee42833d5fe1a354743e07bf629ce11"
)

# The standard genetic code. C. auris uses the standard nuclear code for ERG11 (it is
# NOT a CTG-clade alternative-code gene in the way some Candida serine reassignments are
# discussed -- the reference translation below is checked against the known wild-type
# hotspot residues at load, which would fail if the code were wrong).
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: _AAS[i]
    for i, (a, b, c) in enumerate(
        (a, b, c) for a in _BASES for b in _BASES for c in _BASES
    )
}

# Known C. auris azole-resistance ERG11 substitutions, by residue position. The wild-type
# letter is what the pinned reference MUST show at that position (asserted in
# load_reference). Values are the resistance-associated mutant letters seen in the
# literature; used only to TAG emitted tokens, never to invent them.
RESISTANCE_PANEL = {
    125: ("V", {"A"}),   # V125A -- part of the VF125AL clade-III haplotype
    126: ("F", {"L"}),   # F126L (a.k.a. F126T in some reports; L is the common azole call)
    132: ("Y", {"F"}),   # Y132F -- the dominant global resistance mutation (clades I, IV)
    143: ("K", {"R"}),   # K143R
}

# Bases we treat as confidently called. Anything else in a codon (N, gap, IUPAC
# ambiguity code) makes that codon uncalled -- see honesty constraint 2.
_CONFIDENT = set("ACGT")

# A canonical point-substitution token, e.g. 'Y132F': wild-type, position, mutant.
_TOKEN_RE = re.compile(r"^([A-Za-z])(\d+)([A-Za-z*])$")


def is_panel_token(token):
    """True if `token` (e.g. 'Y132F') is a KNOWN azole-resistance panel substitution.

    Mirrors call_substitutions' panel_hits test exactly -- same wild-type-letter and
    mutant-letter check against RESISTANCE_PANEL -- so a prevalence count rebuilt from
    stored `erg11_call` strings tags precisely the tokens the caller tagged. A substitution
    at a panel position but with a NON-panel mutant letter (e.g. Y132H) is not a panel hit.
    """
    m = _TOKEN_RE.match((token or "").strip())
    if not m:
        return False
    wt, pos, mut = m.group(1).upper(), int(m.group(2)), m.group(3).upper()
    panel = RESISTANCE_PANEL.get(pos)
    return bool(panel and wt == panel[0] and mut in panel[1])


class ReferenceError(ValueError):
    """The reference sequence is missing, malformed, or fails its numbering check."""


class ConsensusError(ValueError):
    """The isolate consensus cannot be placed in the reference frame (length/frame)."""


def _read_fasta_seq(path):
    """Return the single upper-cased sequence in a FASTA file (concatenating wrapped
    lines). Raises ReferenceError if the file is absent or holds no sequence."""
    if not os.path.exists(path):
        raise ReferenceError(
            f"reference not found: {path}\n"
            "  fetch it with the recipe in "
            "data/earlywarning/erg11_reference/PROVENANCE.md"
        )
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq.append(line.strip().upper())
    s = "".join(seq)
    if not s:
        raise ReferenceError(f"reference file has no sequence: {path}")
    return s


def translate(cds):
    """Translate an in-frame CDS to a protein string.

    A codon that is not three confidently-called bases (A/C/G/T) becomes 'X' -- an
    uncalled residue -- rather than being guessed. The trailing stop codon (if present)
    translates to '*'. Length need not be a multiple of 3 here (a trailing partial codon
    is dropped); frame validation is the caller's job (see call_substitutions).
    """
    out = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i + 3]
        if set(codon) <= _CONFIDENT:
            out.append(CODON_TABLE[codon])
        else:
            out.append("X")   # uncalled -- ambiguous/gap base in this codon
    return "".join(out)


def load_reference(path=DEFAULT_REFERENCE, allow_hash_mismatch=False):
    """Load + verify the pinned ERG11 reference. Returns (cds, protein).

    Verifies (a) the CDS content hash against REFERENCE_CDS_SHA256 and (b) that the
    translated hotspot residues match RESISTANCE_PANEL's wild-type letters -- the
    numbering check that makes an emitted 'Y132F' trustworthy. Either failure raises
    ReferenceError (the hash check can be waived with allow_hash_mismatch for a
    deliberate reference update, but the numbering check never can).
    """
    cds = _read_fasta_seq(path)
    if set(cds) - _CONFIDENT:
        raise ReferenceError(
            f"reference CDS has non-ACGT bases: {sorted(set(cds) - _CONFIDENT)}"
        )
    digest = hashlib.sha256(cds.encode()).hexdigest()
    if digest != REFERENCE_CDS_SHA256 and not allow_hash_mismatch:
        raise ReferenceError(
            "reference CDS hash mismatch -- refusing to call against an unpinned "
            f"reference.\n  expected {REFERENCE_CDS_SHA256}\n  got      {digest}\n"
            "  if this is a deliberate reference update, update REFERENCE_CDS_SHA256 "
            "and PROVENANCE.md (and re-verify the numbering)."
        )
    if len(cds) % 3 != 0:
        raise ReferenceError(f"reference CDS length {len(cds)} is not a multiple of 3")
    protein = translate(cds).rstrip("*")
    for pos, (wt, _muts) in sorted(RESISTANCE_PANEL.items()):
        if pos > len(protein):
            raise ReferenceError(
                f"reference protein is too short ({len(protein)} aa) for panel position {pos}"
            )
        observed = protein[pos - 1]
        if observed != wt:
            raise ReferenceError(
                f"numbering check FAILED: panel expects {wt} at residue {pos}, "
                f"reference has {observed}. The reference frame does not match the "
                "resistance-panel numbering -- refusing to call."
            )
    return cds, protein


def _token(wt, pos, mut):
    return f"{wt}{pos}{mut}"


def call_substitutions(consensus_cds, reference=None):
    """Call ERG11 substitutions from an isolate's in-frame consensus CDS.

    Returns a dict:
      tokens            sorted list of substitution tokens (e.g. ['K143R', 'Y132F'])
      panel_hits        subset of tokens at known resistance-panel positions, tagged
      uncalled_positions sorted residue positions that were uncalled (ambiguous codon)
      uncalled_panel    subset of uncalled_positions that fall on panel residues
      n_covered         number of residues confidently called
      n_residues        reference protein length

    `reference` is an optional pre-loaded (cds, protein) tuple; if omitted the pinned
    reference is loaded. Raises ConsensusError if the consensus cannot be placed in the
    reference frame (honesty constraint 3).
    """
    if reference is None:
        reference = load_reference()
    _ref_cds, ref_prot = reference

    cds = (consensus_cds or "").strip().upper().replace("U", "T")
    if not cds:
        raise ConsensusError("empty consensus")
    if len(cds) % 3 != 0:
        raise ConsensusError(
            f"consensus length {len(cds)} is not a multiple of 3 -- a frameshift/indel "
            "cannot be represented as point substitutions (refused, not frame-guessed)."
        )
    iso_prot = translate(cds).rstrip("*")
    if len(iso_prot) != len(ref_prot):
        raise ConsensusError(
            f"consensus protein length {len(iso_prot)} != reference {len(ref_prot)}; "
            "an in-frame full-length ERG11 CDS is required (indels are out of scope)."
        )

    tokens = []
    panel_hits = []
    uncalled = []
    n_covered = 0
    for i, (r, o) in enumerate(zip(ref_prot, iso_prot)):
        pos = i + 1
        if o == "X":                 # uncalled codon -- never counted as wild type
            uncalled.append(pos)
            continue
        n_covered += 1
        if o == r:                   # confidently wild type
            continue
        tok = _token(r, pos, o)
        tokens.append(tok)
        panel = RESISTANCE_PANEL.get(pos)
        if panel and r == panel[0] and o in panel[1]:
            panel_hits.append(tok)

    uncalled_panel = [p for p in uncalled if p in RESISTANCE_PANEL]
    return {
        "tokens": sorted(tokens, key=lambda t: int(re.search(r"\d+", t).group())),
        "panel_hits": sorted(panel_hits, key=lambda t: int(re.search(r"\d+", t).group())),
        "uncalled_positions": uncalled,
        "uncalled_panel": uncalled_panel,
        "n_covered": n_covered,
        "n_residues": len(ref_prot),
    }


def format_call(result):
    """Render call_substitutions() output into an (erg11_call, resistance_source) pair
    ready to drop into a snapshot row.

    - erg11_call is the comma-joined token string emergence.parse_call consumes
      ('Y132F,K143R'), or '' when no substitution was called.
    - resistance_source records HOW and how completely it was called, so an empty call
      is never ambiguous between "sequenced, wild type" and "could not sequence":
        sra-recaller:wild-type            -- full panel covered, no substitution
        sra-recaller:called               -- substitution(s) called
        sra-recaller:partial(<positions>) -- one or more panel residues uncalled
    """
    call = ",".join(result["tokens"])
    if result["uncalled_panel"]:
        pos = "+".join(str(p) for p in result["uncalled_panel"])
        source = f"sra-recaller:partial(panel-uncalled:{pos})"
    elif result["tokens"]:
        source = "sra-recaller:called"
    else:
        source = "sra-recaller:wild-type"
    return call, source
