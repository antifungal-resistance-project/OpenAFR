"""FKS1 structural proximity and contact label (moonshot #4, Rungs 2 + 3).

Rung 2 — `structural_proximity(token)`: minimum heavy-atom distance from each
novel FKS1 hot-spot position to the enfumafungin binding site in the 9WZU
cryo-EM structure (*S. cerevisiae* aligned to *C. auris*).

Rung 3 — `structural_verdict(token)`: a certifiable **structural contact label**
for tokens in `drug_contact` or `pocket_adjacent` proximity tier.  The label
claims spatial co-location only — NOT a resistance prediction, NOT ΔΔG, NOT a
probability. `evidence_class: "proximity_label"` explicitly must NOT trigger
ACT-NOW in alert.py (that tier requires `"docking"` or `"measured"`).

Both functions return None for core panel markers (those carry
RESISTANCE_MARKER_DETECTED already), tokens outside hot-spot windows, and
malformed tokens.

Pre-registrations:
  Rung 2 — work/PREREGISTRATION_fks1_structural_proximity.md (frozen 2026-10-17)
  Rung 3 — work/PREREGISTRATION_fks1_structural_verdict.md (frozen 2026-10-18)
  Cert result — work/RESULTS_fks1_structural_verdict.md (PASS: 4/4 core panel
                positions in drug_contact/pocket_adjacent, 100% ≥ 80% bar)

Stdlib only, offline, pure. The contact map (data/structures/fks1_contact_map.json)
is pre-computed by scripts/derive_fks1_contact_map.py and ships pinned in the repo.
"""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_CONTACT_MAP_PATH = os.path.join(ROOT, "data", "structures", "fks1_contact_map.json")

# Pinned SHA-256 of the contact map, matching PREREG_fks1_structural_proximity.sha256. This is
# a test-pinned integrity constant (asserted by tests/test_fks1_structural.py), not verified at
# module load -- _load_map just reads the committed JSON.
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
    from openafr.fks1_caller import panel_tier, parse_token
    from openafr.fks1_position import _window_of

    parsed = parse_token(token)
    if not parsed:
        return None
    wt, pos, mut = parsed

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


def structural_verdict(token):
    """Structural contact label for a FKS1 VUS token (moonshot #4, Rung 3), or None.

    Returns a verdict dict for tokens in `drug_contact` or `pocket_adjacent` proximity
    tier; None for all others (including core panel markers, `pocket_proximate`, `distal`,
    `no_structure`, out-of-window, malformed).

    The returned dict:
      evidence_class      "proximity_label" — NOT "docking" or "measured"; must NOT
                          trigger ACT-NOW in alert.py
      proximity_tier      "drug_contact" | "pocket_adjacent"
      min_dist_to_drug_A  float Angstroms
      label               human-readable co-location claim (spatial only)
      note                one sentence; states non-claim explicitly

    Pre-registered cert PASS (4/4 core panel positions in contact zone, 100% ≥ 80% bar).
    See work/PREREGISTRATION_fks1_structural_verdict.md and
    work/RESULTS_fks1_structural_verdict.md.
    """
    prox = structural_proximity(token)
    if prox is None:
        return None
    tier = prox["proximity_tier"]
    dist = prox["min_dist_to_drug_A"]
    pos = prox["position"]
    window = prox["window"]
    dist_str = f"{dist:.1f} Å" if dist is not None else "N/A"
    src = prox["reference_structure"]
    caveat = f"S. cerevisiae proxy, ±1–2 Å; spatial co-location only — not a resistance prediction"
    if tier == DRUG_CONTACT:
        label = "at drug-contact residue"
        note = (
            f"Residue {pos} in {window} sits within direct contact distance of the "
            f"enfumafungin binding site ({dist_str} in {src}; {caveat})."
        )
    elif tier == POCKET_ADJACENT:
        label = "adjacent to drug-contact residues"
        note = (
            f"Residue {pos} in {window} is adjacent to the enfumafungin binding site "
            f"({dist_str} in {src}; {caveat})."
        )
    else:
        return None
    return {
        "evidence_class": "proximity_label",
        "proximity_tier": tier,
        "min_dist_to_drug_A": dist,
        "label": label,
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
