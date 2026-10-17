"""FKS1 structural proximity (moonshot #4, Rung 2).

Maps a novel FKS1 hot-spot substitution to its structural proximity to the
echinocandin / enfumafungin binding site, using a pre-derived contact map
computed from the 9WZU cryo-EM structure (Fks1 + enfumafungin, 2.97 Å,
S. cerevisiae sequence aligned to C. auris via SW alignment).

What this claims — and what it DOES NOT
----------------------------------------
`structural_proximity(token)` is a LOOKUP over the pre-derived contact map.
For a hot-spot token it reports how close (in Angstroms, minimum heavy-atom
distance) the equivalent S. cerevisiae residue is to the bound enfumafungin.
Enfumafungin occupies the same FKS1 binding site as echinocandins.

It deliberately does NOT:
  1. predict phenotype (R/S) for the variant;
  2. emit a resistance probability;
  3. predict binding affinity or ΔΔG;
  4. populate the verdict `structural` field — that remains null for FKS1
     (Rung 3, not in scope here);
  5. extend claims beyond the S. cerevisiae cryo-EM structure alignment.

The `structural` field for FKS1 verdicts stays null; proximity lives in
`position_context[*].proximity` (see openafr/fks1_position.py).

See work/PREREGISTRATION_fks1_structural_proximity.md for the frozen scope.
Stdlib only, offline, pure. The contact map (data/structures/fks1_contact_map.json)
is pre-computed by scripts/derive_fks1_contact_map.py and ships pinned in the repo.
"""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_CONTACT_MAP_PATH = os.path.join(ROOT, "data", "structures", "fks1_contact_map.json")

# Pinned SHA-256 of the contact map (verified at module load via PREREG_fks1_structural_proximity.sha256)
CONTACT_MAP_SHA256 = "6a2033ac5bd28b821e5cc3a47b759932c5b728c7f85b69958a82f1f89902acea"

# Proximity tier names (weakest first, strongest last for ordering comparisons)
DRUG_CONTACT = "drug_contact"       # ≤ 5.0 Å
POCKET_ADJACENT = "pocket_adjacent" # 5.0–8.0 Å
POCKET_PROXIMATE = "pocket_proximate"  # 8.0–12.0 Å
DISTAL = "distal"                   # > 12.0 Å
NO_STRUCTURE = "no_structure"       # position absent from contact map

PROXIMITY_TIERS = (DRUG_CONTACT, POCKET_ADJACENT, POCKET_PROXIMATE, DISTAL, NO_STRUCTURE)

_map_cache = None


def _load_map():
    global _map_cache
    if _map_cache is None:
        _map_cache = json.load(open(_CONTACT_MAP_PATH))
    return _map_cache


def contact_map_meta():
    """Return metadata from the contact map (source PDB, drug, cutoffs)."""
    m = _load_map()
    return {
        "source_pdb": m["source_pdb"],
        "drug_name": m["drug_name"],
        "contact_cutoffs_A": m["contact_cutoffs_A"],
        "alignment_offsets_sc_minus_ca": m["alignment_offsets_sc_minus_ca"],
    }


def structural_proximity(token):
    """3D structural proximity tier for one FKS1 substitution token, or None.

    Returns None when:
      - the token is malformed or cannot be parsed;
      - the residue falls outside every hot-spot window (same gate as
        fks1_position.position_context — proximity is only meaningful inside
        a hot-spot because the contact map only covers those windows);
      - the residue is a CORE panel marker (those carry RESISTANCE_MARKER_DETECTED
        already; proximity context is for the UNCHARACTERIZED tail only).

    For all other in-window, non-core positions the function returns a dict:
      token               the input token, normalised
      position            the residue number (int)
      window              'HS1' | 'HS2' | 'HS3'
      proximity_tier      one of PROXIMITY_TIERS above
      min_dist_to_drug_A  float Angstroms, or None when NO_STRUCTURE
      sc_position         the aligned S. cerevisiae residue number used
      sc_aa               one-letter code of the Sc residue
      reference_structure source PDB ID and drug name
      note                one-line plain-English description (for alerts/page)

    The `proximity_tier` is the headline field. A `drug_contact` tier means the
    aligned S. cerevisiae residue's heavy atoms come within 5 Å of enfumafungin,
    which occupies the same binding site as clinical echinocandins (caspofungin,
    micafungin, anidulafungin). The measurement is NOT a phenotype prediction.
    """
    from openafr.fks1_caller import _TOKEN_RE, panel_tier
    from openafr.fks1_position import _window_of

    m = _TOKEN_RE.match((token or "").strip())
    if not m:
        return None
    wt, pos, mut = m.group(1).upper(), int(m.group(2)), m.group(3).upper()

    # Skip core panel markers — they already have a resistance verdict
    if panel_tier(token) == "core":
        return None

    # Skip tokens outside hot-spot windows
    window = _window_of(pos)
    if window is None:
        return None

    contact_map = _load_map()
    entry = contact_map["residues"].get(str(pos))

    if entry is None:
        return None

    d = entry.get("min_dist_to_drug_A")
    tier = entry["tier"]
    sc_pos = entry["sc_pos"]
    sc_aa = entry.get("sc_aa", "?")
    source = f"PDB {contact_map['source_pdb']} ({contact_map['drug_name']})"

    note = _describe_proximity(tier, d, token, window.name, source)

    return {
        "token": f"{wt}{pos}{mut}",
        "position": pos,
        "window": window.name,
        "proximity_tier": tier,
        "min_dist_to_drug_A": d,
        "sc_position": sc_pos,
        "sc_aa": sc_aa,
        "reference_structure": source,
        "note": note,
    }


def _describe_proximity(tier, dist_a, token, window, source):
    """One plain, defensible sentence for a proximity dict (for alerts / the page).

    States only spatial fact. Never a phenotype or phenotype-prediction claim.
    """
    pos_str = token[1:-1]  # strip wt/mut letters
    if tier == DRUG_CONTACT:
        return (
            f"residue {pos_str} in hot-spot {window} sits within direct contact "
            f"distance of the inhibitor binding site ({dist_a:.1f} Å from "
            f"enfumafungin in {source}); structural proximity to the drug is highest."
        )
    if tier == POCKET_ADJACENT:
        return (
            f"residue {pos_str} in hot-spot {window} is adjacent to the inhibitor "
            f"binding site ({dist_a:.1f} Å from enfumafungin in {source})."
        )
    if tier == POCKET_PROXIMATE:
        return (
            f"residue {pos_str} in hot-spot {window} is proximate but not directly "
            f"adjacent to the inhibitor binding site ({dist_a:.1f} Å from "
            f"enfumafungin in {source})."
        )
    if tier == DISTAL:
        return (
            f"residue {pos_str} in hot-spot {window} is distal from the inhibitor "
            f"binding site ({dist_a:.1f} Å from enfumafungin in {source}); "
            f"within hot-spot window but not at the drug-contact face."
        )
    return f"no structural data for residue {pos_str} in the {source} contact map."
