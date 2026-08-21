"""FKS1 re-caller -- reads/consensus -> echinocandin-resistance substitution call (v2).

Pipeline stage: the detection half of the FKS1/echinocandin early-warning track (v2).
The azole/ERG11 track (openafr/recaller.py) measured an ~80% azole event frequency, which
saturates the azole *emergence* signal (work/RESULTS_prevalence.md); echinocandin (FKS1)
resistance is the dynamic axis with headroom. This module is the FKS1 analog of the ERG11
re-caller: given an isolate's FKS1 consensus over the resistance hot-spots, it emits the
substitution tokens (e.g. 'S639F') that a per-gene emergence/alert pipeline consumes.

It is DETECTION-ONLY by design. FKS1 is a 1,3-beta-glucan synthase; echinocandins do not
coordinate a metal, so the CYP51/heme-iron geometry moat (VISION.md, openafr/structural.py)
does NOT transfer, and no structural "so-what" is claimed here. A downstream FKS1 alert must
say "detection only, no structural verdict", never borrow the azole pocket verdict.

The one design change from the ERG11 re-caller: WINDOWED, not whole-CDS
-----------------------------------------------------------------------
ERG11 is a compact 1575 nt CDS, so recaller.py builds a full-length in-frame consensus and
refuses anything whose length != reference (frameshift = refused, not guessed). FKS1 is a
~5.6 kb gene (1888 aa); all echinocandin resistance lives in two short, conserved hot-spot
windows (HS1 around S639, HS2 around R1354). Calling the whole gene buys nothing clinically
and makes the length-exact frame contract fragile across 5.6 kb -- one low-coverage indel
anywhere would discard an isolate whose hot-spots were perfectly covered. So this caller
operates PER WINDOW: each hot-spot is validated and called independently, and a frameshift
inside one window refuses *that window* with a reason without discarding the other.

The honesty constraints (same bar as recaller.py, applied per window)
---------------------------------------------------------------------
1. WE ASSERT THE REFERENCE, WE DO NOT ASSUME IT. Every call is relative to the pinned
   B8441 GSC1/FKS1 reference (data/earlywarning/fks1_reference/, accession XM_085597048.1,
   sha256 verified at load). load_reference() refuses a reference whose window anchor
   residues (S639, R1354) do not translate to the expected letters -- a numbering or
   sequence swap is surfaced, never silently mis-called. (This numbering was verified by
   translation, not taken on faith -- see PROVENANCE.md; an early wrong-gene guess, locus
   B9J08_004112, was caught exactly this way.)
2. AN UNCALLED CODON IS UNCALLED, NEVER "WILD TYPE". A codon with an ambiguous/gap base is
   reported uncalled for that residue; it never contributes a wild-type "no change".
3. FRAME IS ASSERTED PER WINDOW, NOT GUESSED. A window consensus must be the exact window
   length and in frame; a length mismatch (an indel inside the window) is refused for that
   window with the reason -- never forced into a frame, never fatal to the other window.
4. WE TAG THE PANEL, WE DO NOT INVENT NOVELTY. Every non-synonymous difference in a window
   is emitted verbatim. Known-resistance panel substitutions (S639F/P/Y) are TAGGED; the
   HS2 window is reported but not yet panel-tagged (its mutant letters await a literature
   pin -- a documented gap, not a silent drop), so HS2 substitutions are emitted untagged
   and it is emergence's job to decide which are novel/rising.

Stdlib only. Operates on the reference FASTA and on per-window consensus nucleotide strings.
"""
import hashlib
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFERENCE_DIR = os.path.join(ROOT, "data", "earlywarning", "fks1_reference")
DEFAULT_REFERENCE = os.path.join(REFERENCE_DIR, "fks1_cds.fasta")

# The pinned reference CDS content hash (bare ACGT string, no header/newlines).
# Must match data/earlywarning/fks1_reference/PROVENANCE.md. Same freeze contract as the
# ERG11 reference, protocol.yaml, and the snapshot INDEX hashes.
REFERENCE_CDS_SHA256 = (
    "d50ececa0ff8977b61fce9422d2fc255a5e1d97144dfe6fcb24bb6501fbf15ee"
)

# Standard nuclear genetic code (FKS1/GSC1 uses the standard code; the reference
# translation is checked against the known window anchors at load, which would fail if
# the code were wrong).
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: _AAS[i]
    for i, (a, b, c) in enumerate(
        (a, b, c) for a in _BASES for b in _BASES for c in _BASES
    )
}

_CONFIDENT = set("ACGT")
_TOKEN_RE = re.compile(r"^([A-Za-z])(\d+)([A-Za-z*])$")


class Window:
    """A FKS1 resistance hot-spot: an inclusive 1-based residue range + an anchor residue
    whose wild-type letter the reference MUST translate to (the per-window numbering check).
    `panel` maps a residue position -> set of known-resistance mutant letters to TAG; an
    empty panel means the window is reported but not yet panel-tagged (HS2)."""

    def __init__(self, name, start_res, end_res, anchor_pos, anchor_wt, panel):
        self.name = name
        self.start_res = start_res
        self.end_res = end_res
        self.anchor_pos = anchor_pos
        self.anchor_wt = anchor_wt
        self.panel = panel   # {residue_pos: (wild_type_letter, {mutant letters to tag})}

    @property
    def n_residues(self):
        return self.end_res - self.start_res + 1

    @property
    def nt_len(self):
        return self.n_residues * 3

    def nt_slice(self, cds):
        return cds[(self.start_res - 1) * 3: self.end_res * 3]

    def ref_protein(self, ref_prot):
        return ref_prot[self.start_res - 1: self.end_res]


# The two echinocandin-resistance hot-spots, numbered against XM_085597048.1 and VERIFIED
# by translation (PROVENANCE.md):
#   HS1  residues 635-643  FLTLSLRDP  -- panel S639F/S639P/S639Y (the dominant C. auris call)
#   HS2  residues 1350-1358 DWIRRYTLS -- anchor R1354; mutant letters NOT yet panel-encoded
#     (reported verbatim/untagged until a literature pin fixes the HS2 mutant set -- see
#      constraint 4). Widening the panel is a data edit here, no code change downstream.
FKS1_WINDOWS = {
    "HS1": Window("HS1", 635, 643, 639, "S", {639: ("S", {"F", "P", "Y"})}),
    "HS2": Window("HS2", 1350, 1358, 1354, "R", {}),
}

# Flat position -> (wt, mutant-letters) panel across all windows, for is_panel_token and
# tagging. Rebuilt from the windows so the two can never disagree.
RESISTANCE_PANEL = {
    pos: entry
    for w in FKS1_WINDOWS.values()
    for pos, entry in w.panel.items()
}


class ReferenceError(ValueError):
    """The reference sequence is missing, malformed, or fails its numbering check."""


class WindowError(ValueError):
    """A window consensus cannot be placed in the reference frame (length/frame)."""


def _read_fasta_seq(path):
    if not os.path.exists(path):
        raise ReferenceError(
            f"reference not found: {path}\n"
            "  fetch it with the recipe in "
            "data/earlywarning/fks1_reference/PROVENANCE.md"
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
    """Translate an in-frame CDS to protein. A codon that is not three confidently-called
    bases becomes 'X' (uncalled), never guessed. A trailing partial codon is dropped."""
    out = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i + 3]
        out.append(CODON_TABLE[codon] if set(codon) <= _CONFIDENT else "X")
    return "".join(out)


def load_reference(path=DEFAULT_REFERENCE, allow_hash_mismatch=False):
    """Load + verify the pinned FKS1 reference. Returns (cds, protein).

    Verifies (a) the CDS content hash against REFERENCE_CDS_SHA256 and (b) that every
    window's anchor residue translates to its expected wild-type letter (S639, R1354) --
    the numbering check that makes an emitted 'S639F' trustworthy. Either failure raises
    ReferenceError (the hash check can be waived with allow_hash_mismatch for a deliberate
    reference update, but the numbering check never can)."""
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
            "  if this is a deliberate reference update, update REFERENCE_CDS_SHA256 and "
            "PROVENANCE.md (and re-verify the numbering)."
        )
    if len(cds) % 3 != 0:
        raise ReferenceError(f"reference CDS length {len(cds)} is not a multiple of 3")
    protein = translate(cds).rstrip("*")
    for w in FKS1_WINDOWS.values():
        if w.end_res > len(protein):
            raise ReferenceError(
                f"reference protein ({len(protein)} aa) is too short for window "
                f"{w.name} ({w.start_res}-{w.end_res})"
            )
        observed = protein[w.anchor_pos - 1]
        if observed != w.anchor_wt:
            raise ReferenceError(
                f"numbering check FAILED: window {w.name} expects {w.anchor_wt} at "
                f"residue {w.anchor_pos}, reference has {observed}. The reference frame "
                "does not match the resistance numbering -- refusing to call."
            )
    return cds, protein


def is_panel_token(token):
    """True if `token` (e.g. 'S639F') is a KNOWN echinocandin-resistance panel substitution.
    Mirrors call_windows' tagging exactly, so a prevalence count rebuilt from stored
    fks1_call strings tags precisely the tokens the caller tagged."""
    m = _TOKEN_RE.match((token or "").strip())
    if not m:
        return False
    wt, pos, mut = m.group(1).upper(), int(m.group(2)), m.group(3).upper()
    entry = RESISTANCE_PANEL.get(pos)
    return bool(entry and wt == entry[0] and mut in entry[1])


def _token(wt, pos, mut):
    return f"{wt}{pos}{mut}"


def call_window(window, consensus_window_nt, reference=None):
    """Call substitutions within a single hot-spot window.

    `consensus_window_nt` is the isolate's consensus nucleotides spanning exactly this
    window (window.nt_len bases, in frame, reference-anchored -- the orchestration's job to
    extract). Returns a per-window dict:
      status      'wild-type' | 'called' | 'partial' | 'refused'
      tokens      substitution tokens at ABSOLUTE residue positions (e.g. ['S639F'])
      panel_hits  subset of tokens tagged as known resistance
      uncalled_positions  absolute residue positions that were uncalled
      refused     reason string if status == 'refused', else None
      n_covered / n_residues  coverage within the window

    A frameshift/length mismatch (indel inside the window) is REFUSED with a reason, never
    frame-guessed (constraint 3)."""
    if reference is None:
        reference = load_reference()
    _ref_cds, ref_prot = reference
    ref_win = window.ref_protein(ref_prot)

    def refused(reason):
        return {
            "status": "refused", "tokens": [], "panel_hits": [],
            "uncalled_positions": [], "refused": reason,
            "n_covered": 0, "n_residues": window.n_residues,
        }

    seq = (consensus_window_nt or "").strip().upper().replace("U", "T")
    if not seq:
        return refused("empty window consensus")
    if len(seq) % 3 != 0:
        return refused(
            f"length {len(seq)} not a multiple of 3 -- frameshift/indel, not frame-guessed"
        )
    if len(seq) != window.nt_len:
        return refused(
            f"length {len(seq)} != expected {window.nt_len} for {window.name} "
            f"({window.start_res}-{window.end_res}) -- indel in window, out of scope"
        )
    iso_win = translate(seq)
    tokens, panel_hits, uncalled = [], [], []
    n_covered = 0
    for i, (r, o) in enumerate(zip(ref_win, iso_win)):
        pos = window.start_res + i
        if o == "X":                       # uncalled codon -- never counted as wild type
            uncalled.append(pos)
            continue
        n_covered += 1
        if o == r:
            continue
        tok = _token(r, pos, o)
        tokens.append(tok)
        entry = window.panel.get(pos)
        if entry and r == entry[0] and o in entry[1]:
            panel_hits.append(tok)

    if uncalled:
        status = "partial"
    elif tokens:
        status = "called"
    else:
        status = "wild-type"
    return {
        "status": status, "tokens": tokens, "panel_hits": panel_hits,
        "uncalled_positions": uncalled, "refused": None,
        "n_covered": n_covered, "n_residues": window.n_residues,
    }


def call_windows(window_seqs, reference=None):
    """Call every FKS1 hot-spot from a {window_name: consensus_nt} mapping.

    A window absent from `window_seqs` is reported with status 'missing' (we have no data
    for it) -- never silently assumed wild-type. Returns a dict:
      tokens        sorted absolute tokens across all called windows
      panel_hits    subset tagged as known resistance
      uncalled_panel  panel residue positions that were uncalled
      windows       {name: <call_window result, or {'status':'missing',...}>}
    """
    if reference is None:
        reference = load_reference()
    windows = {}
    for name, w in FKS1_WINDOWS.items():
        if name not in window_seqs or window_seqs.get(name) is None:
            windows[name] = {
                "status": "missing", "tokens": [], "panel_hits": [],
                "uncalled_positions": [], "refused": None,
                "n_covered": 0, "n_residues": w.n_residues,
            }
        else:
            windows[name] = call_window(w, window_seqs[name], reference=reference)

    tokens, panel_hits, uncalled_panel = [], [], []
    for name, res in windows.items():
        tokens.extend(res["tokens"])
        panel_hits.extend(res["panel_hits"])
        for p in res["uncalled_positions"]:
            if p in RESISTANCE_PANEL:
                uncalled_panel.append(p)
    key = lambda t: int(re.search(r"\d+", t).group())
    return {
        "tokens": sorted(tokens, key=key),
        "panel_hits": sorted(panel_hits, key=key),
        "uncalled_panel": sorted(uncalled_panel),
        "windows": windows,
    }


def slice_windows_from_cds(full_cds):
    """Convenience bridge: extract the per-window consensus nucleotides from a full-length,
    in-frame consensus CDS (reference length). This is the FRAGILE path -- it requires the
    whole 5.6 kb consensus to be length-exact -- and exists only for tests / a caller that
    already has a full consensus. The robust orchestration path extracts each window by
    reference coordinate instead. Raises WindowError if the consensus cannot be anchored."""
    reference = load_reference()
    ref_cds = reference[0]
    cds = (full_cds or "").strip().upper().replace("U", "T")
    if len(cds) != len(ref_cds):
        raise WindowError(
            f"full consensus length {len(cds)} != reference {len(ref_cds)}; cannot anchor "
            "windows (use the per-window path instead)."
        )
    return {name: w.nt_slice(cds) for name, w in FKS1_WINDOWS.items()}


def format_call(result):
    """Render call_windows() output into an (fks1_call, resistance_source) snapshot pair.

    - fks1_call is the comma-joined token string ('S639F' or 'S639F,R1354S'), '' if none.
    - resistance_source records HOW and how completely EACH window was called, so an empty
      call is never ambiguous between 'sequenced, wild type' and 'could not sequence':
        sra-fks1-recaller:HS1=wild-type,HS2=wild-type   -- both windows covered, no change
        sra-fks1-recaller:HS1=called,HS2=wild-type      -- a substitution in HS1
        sra-fks1-recaller:HS1=partial(uncalled:639),... -- a hot-spot residue uncalled
        sra-fks1-recaller:HS1=refused(...),HS2=...       -- an indel inside a window
        sra-fks1-recaller:HS1=missing,HS2=...            -- no consensus for that window
    Per-window status (not ERG11's single status) because FKS1 has two independent windows.
    """
    call = ",".join(result["tokens"])
    parts = []
    for name in FKS1_WINDOWS:
        res = result["windows"][name]
        st = res["status"]
        if st == "partial":
            st = "partial(uncalled:%s)" % "+".join(
                str(p) for p in res["uncalled_positions"]
            )
        elif st == "refused":
            st = "refused(%s)" % res["refused"]
        parts.append(f"{name}={st}")
    return call, "sra-fks1-recaller:" + ",".join(parts)
