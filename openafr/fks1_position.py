"""FKS1 novel-variant POSITION CONTEXT (moonshot #4, Rungs 1 + 2).

The problem this closes
-----------------------
A never-seen substitution inside an FKS1 resistance hot-spot flows through the caller and
emergence detector as a bare `UNCHARACTERIZED_VARIANT` (openafr/verdict.py) with
`structural = None` and *no context at all*. A surveillance reader cannot tell whether the
novel change lands squarely on the S639 anchor -- where F/P/Y are canonical echinocandin
resistance markers -- or three residues away in benign hot-spot background. The caller is
blind on first sighting of exactly the variants surveillance exists to catch.

What this module claims -- and the four things it does NOT
--------------------------------------------------------
`position_context(token)` is a DETERMINISTIC DERIVATION over the already-pinned panel and
reference numbering (openafr/fks1_caller.py: FKS1_WINDOWS / RESISTANCE_PANEL / LOW_PPV_PANEL).
It answers "where does this novel residue sit relative to the known resistance geography?"
and NOTHING more. It deliberately does NOT:

  1. claim a drug-FIT change -- the proximity tier (Rung 2) says where the residue sits
     in 3D space relative to the drug, NOT whether the drug still fits;
  2. predict a phenotype (R/S) for the novel letter;
  3. emit a probability;
  4. populate the verdict's `structural` field -- that stays null for FKS1, by design.

So this is strictly MORE than the current `None`, while staying inside the project's
detection-only boundary (openafr/alert.py, data/earlywarning/fks1_reference/PROVENANCE.md).
See work/PREREGISTRATION_fks1_position_context.md for the frozen scope + non-claims.

Stdlib only, pure, offline. Reuses the caller's panel so the two can never disagree.
"""
from openafr.fks1_caller import (
    FKS1_WINDOWS,
    LOW_PPV_PANEL,
    RESISTANCE_PANEL,
    _TOKEN_RE,
    panel_tier,
)

# Position-context tiers, strongest signal first. Derived from the pinned panel, never
# invented -- each says only WHERE the residue sits / WHICH known class it falls in, not
# what it does to the drug.
KNOWN_LOW_PPV_MARKER = "known_low_ppv_marker"            # the token IS a low-PPV panel letter
NOVEL_AT_PANEL_RESIDUE = "novel_substitution_at_panel_residue"  # novel letter, panel residue
ADJACENT_TO_PANEL = "in_hotspot_adjacent_to_panel"
HOTSPOT_BACKGROUND = "in_hotspot_background"

# A novel change within this many residues of a core panel marker is "adjacent"; beyond it
# is hot-spot background. Hot-spot windows are ~9-11 residues wide, so a 2-residue radius
# distinguishes "right next to a known marker" from "elsewhere in the window" without
# implying any mechanistic contact (which we explicitly do not claim).
ADJACENT_RADIUS = 2


def _parse(token):
    """(wt, pos, mut) for a well-formed substitution token, else None. Uppercased."""
    m = _TOKEN_RE.match((token or "").strip())
    if not m:
        return None
    return m.group(1).upper(), int(m.group(2)), m.group(3).upper()


def _window_of(pos):
    """The FKS1_WINDOWS entry whose residue range contains `pos`, or None."""
    for w in FKS1_WINDOWS.values():
        if w.start_res <= pos <= w.end_res:
            return w
    return None


def _core_panel_positions():
    """Residue positions carrying at least one CORE (high-PPV) panel letter.

    A position is core if it has a panel letter that is NOT in LOW_PPV_PANEL -- these are the
    canonical resistance markers (F635, S639, R1354, W691); D642/M690 are low-PPV only and
    are deliberately excluded so "distance to nearest marker" measures distance to a marker
    we would actually call, not to a position we abstain on."""
    out = set()
    for pos, (_wt, letters) in RESISTANCE_PANEL.items():
        if letters - LOW_PPV_PANEL.get(pos, set()):
            out.add(pos)
    return out


_CORE_POSITIONS = _core_panel_positions()


def position_context(token):
    """Position context for a single FKS1 substitution token, or None.

    Returns None when the token is not a well-formed substitution, falls outside every
    hot-spot window, OR is itself a CORE panel marker (those are already
    RESISTANCE_MARKER_DETECTED -- position context is for the UNcharacterized tail). A
    low-PPV panel letter (D642Y, M690I) is NOT core, so it DOES get context: the verdict
    layer abstains on it, and a reader deserves to see it sits on a known-but-low-PPV
    residue rather than in blank background.

    The returned dict (all fields derived, never invented):
      token                  the input token, normalised
      window                 'HS1' | 'HS2' | 'HS3'
      anchor                 the window anchor as '<wt><pos>' (e.g. 'S639')
      anchor_distance        signed residue offset from the anchor (pos - anchor_pos)
      at_panel_residue       True if this residue carries ANY panel letter
      panel_letters_at_position  sorted known panel mutant letters at this residue ([] if none)
      nearest_core_marker    nearest core-marker residue as '<pos>' (e.g. '639'), or None
      distance_to_nearest    |pos - nearest core position|, or None if no core marker in reach
      tier                   NOVEL_AT_PANEL_RESIDUE | ADJACENT_TO_PANEL | HOTSPOT_BACKGROUND
    """
    parsed = _parse(token)
    if parsed is None:
        return None
    wt, pos, mut = parsed

    tier_of_token = panel_tier(token)  # 'core' | 'low_ppv' | None
    # Core panel markers are already DETECTED; context is for the uncharacterised tail.
    if tier_of_token == "core":
        return None

    window = _window_of(pos)
    if window is None:
        return None

    entry = RESISTANCE_PANEL.get(pos)
    panel_letters = sorted(entry[1]) if entry else []
    at_panel_residue = bool(panel_letters)

    nearest = min(_CORE_POSITIONS, key=lambda p: abs(p - pos), default=None)
    dist = abs(nearest - pos) if nearest is not None else None

    if tier_of_token == "low_ppv":
        # The token itself is a known-but-low-PPV panel letter (D642Y, M690I): the verdict
        # layer abstains on it (carriers split R/S). Not novel -- a distinct, honest tier.
        tier = KNOWN_LOW_PPV_MARKER
    elif at_panel_residue:
        tier = NOVEL_AT_PANEL_RESIDUE
    elif dist is not None and dist <= ADJACENT_RADIUS:
        tier = ADJACENT_TO_PANEL
    else:
        tier = HOTSPOT_BACKGROUND

    ctx = {
        "token": f"{wt}{pos}{mut}",
        "window": window.name,
        "anchor": f"{window.anchor_wt}{window.anchor_pos}",
        "anchor_distance": pos - window.anchor_pos,
        "at_panel_residue": at_panel_residue,
        "panel_letters_at_position": panel_letters,
        "nearest_core_marker": str(nearest) if nearest is not None else None,
        "distance_to_nearest": dist,
        "tier": tier,
    }

    # Rung 2: add structural proximity from the 9WZU cryo-EM contact map when available.
    # Import lazily so Rung 1 tests that do not install the contact map still pass.
    try:
        from openafr import fks1_structural
        prox = fks1_structural.structural_proximity(token)
        if prox is not None:
            ctx["proximity"] = prox
    except Exception:
        pass  # contact map absent or malformed -- degrade gracefully to Rung-1-only

    return ctx


def describe(ctx):
    """One plain, defensible sentence for a position-context dict (for alerts / the page).

    States only geography: which hot-spot, offset from the anchor, and the nearest known
    marker. Never a fit or phenotype claim -- the detection-only caveat still carries that.
    """
    if not ctx:
        return ""
    where = f"hot-spot {ctx['window']}"
    off = ctx["anchor_distance"]
    rel = (f"at the {ctx['anchor']} anchor" if off == 0
           else f"{abs(off)} residue(s) {'after' if off > 0 else 'before'} the "
                f"{ctx['anchor']} anchor")
    if ctx["tier"] == KNOWN_LOW_PPV_MARKER:
        so = ("a known-but-low-PPV FKS1 change (carriers split resistant/susceptible), "
              "abstained on rather than called")
    elif ctx["tier"] == NOVEL_AT_PANEL_RESIDUE:
        letters = "/".join(ctx["panel_letters_at_position"])
        one = len(ctx["panel_letters_at_position"]) == 1
        verb = "is a known echinocandin-resistance change" if one else \
            "are known echinocandin-resistance changes"
        so = (f"a novel substitution at residue {ctx['token'][1:-1]}, where {letters} {verb}")
    elif ctx["tier"] == ADJACENT_TO_PANEL:
        so = (f"adjacent (≤{ADJACENT_RADIUS} aa) to known marker position "
              f"{ctx['nearest_core_marker']}")
    else:
        so = "in hot-spot background, not on or next to a known marker position"
    return f"sits in {where}, {rel} — {so}."
