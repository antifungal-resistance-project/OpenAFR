"""Tests for openafr/fks1_structural.py (moonshot #4, Rung 2).

Correctness and sanity checks for the pre-derived FKS1 contact map and the
structural_proximity() function. No phenotype or resistance claims are tested
here -- the module makes none. The tests verify:

  1. The contact map is present and has the expected SHA-256 (integrity gate).
  2. Known CORE panel residues (F635, S639, W691, R1354) land in drug_contact or
     pocket_adjacent when queried via a novel letter at that position -- the
     structural validation of the alignment (anti-circularity: the contact zone
     comes from the cryo-EM drug coordinates, not from the panel).
  3. A position outside all hot-spot windows returns None.
  4. A core panel marker token (e.g. S639F) returns None (already DETECTED).
  5. A low-PPV panel token (D642Y) returns a proximity dict with the measured tier.
  6. The proximity dict has the required keys and the tier is a valid PROXIMITY_TIER.
  7. The contact map metadata keys are present.
  8. fks1_position.position_context enriches non-core tokens with a `proximity` key.
"""
import hashlib
import os

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONTACT_MAP_PATH = os.path.join(ROOT, "data", "structures", "fks1_contact_map.json")

from openafr.fks1_structural import (
    CONTACT_MAP_SHA256,
    DISTAL,
    DRUG_CONTACT,
    NO_STRUCTURE,
    POCKET_ADJACENT,
    POCKET_PROXIMATE,
    PROXIMITY_TIERS,
    contact_map_meta,
    structural_proximity,
)


# ── integrity gate ─────────────────────────────────────────────────────────────

def test_contact_map_sha256():
    """The shipped contact map must match the pinned SHA-256."""
    actual = hashlib.sha256(open(CONTACT_MAP_PATH, "rb").read()).hexdigest()
    assert actual == CONTACT_MAP_SHA256, (
        f"fks1_contact_map.json SHA-256 mismatch — "
        f"run scripts/derive_fks1_contact_map.py to regenerate\n"
        f"  expected: {CONTACT_MAP_SHA256}\n"
        f"  got:      {actual}"
    )


# ── metadata ───────────────────────────────────────────────────────────────────

def test_contact_map_meta_keys():
    meta = contact_map_meta()
    for key in ("source_pdb", "drug_name", "contact_cutoffs_A", "alignment_offsets_sc_minus_ca"):
        assert key in meta, f"Missing key {key!r} in contact_map_meta()"
    assert meta["source_pdb"] == "9WZU"
    assert "enfumafungin" in meta["drug_name"].lower() or "Enfumafungin" in meta["drug_name"]


# ── None-return gates ──────────────────────────────────────────────────────────

@pytest.mark.parametrize("token", [
    "S639F",   # core panel marker
    "S639P",   # core panel marker
    "F635C",   # core panel marker
    "F635Y",   # core panel marker
    "W691L",   # core panel marker
    "R1354S",  # core panel marker
])
def test_core_panel_returns_none(token):
    """Core panel markers already have RESISTANCE_MARKER_DETECTED — no proximity."""
    assert structural_proximity(token) is None, f"{token} should return None (core panel)"


@pytest.mark.parametrize("token", [
    "A100G",   # nowhere near a hot-spot
    "V200I",
    "G1L",
    "M1888K",  # last residue, no window
    "notareal",
    "",
    None,
])
def test_out_of_window_returns_none(token):
    assert structural_proximity(token) is None


# ── sanity: known panel residue positions land in high-proximity tiers ────────
# Anti-circularity: the contact zone is derived from cryo-EM drug coordinates
# (independent of the panel). Panel residues landing here is independent evidence.

@pytest.mark.parametrize("token,expected_tiers", [
    ("F635L", {DRUG_CONTACT}),           # F635 site: 2.94 Å, drug_contact
    ("F635A", {DRUG_CONTACT}),
    ("L638A", {DRUG_CONTACT}),           # L638: 2.82 Å, drug_contact (non-panel, in-window)
    ("S639A", {POCKET_ADJACENT}),        # S639 site: 7.16 Å, pocket_adjacent
    ("S639L", {POCKET_ADJACENT}),
    ("W691A", {DRUG_CONTACT}),           # W691 site: 3.66 Å, drug_contact
    ("R1354A", {POCKET_ADJACENT}),       # R1354 site: 7.53 Å, pocket_adjacent
])
def test_panel_residue_sites_in_expected_tiers(token, expected_tiers):
    prox = structural_proximity(token)
    assert prox is not None, f"{token} should return a proximity dict"
    assert prox["proximity_tier"] in expected_tiers, (
        f"{token}: expected tier in {expected_tiers}, got {prox['proximity_tier']!r} "
        f"({prox['min_dist_to_drug_A']} Å)"
    )


# ── low-PPV panel tokens get proximity (verdict layer abstains, not proximity) ─

def test_d642y_gets_proximity():
    prox = structural_proximity("D642Y")
    assert prox is not None
    assert prox["proximity_tier"] in (POCKET_PROXIMATE, POCKET_ADJACENT, DRUG_CONTACT)


def test_m690i_gets_proximity():
    prox = structural_proximity("M690I")
    assert prox is not None
    assert prox["proximity_tier"] in (POCKET_PROXIMATE, POCKET_ADJACENT)


# ── return-dict shape ─────────────────────────────────────────────────────────

def test_proximity_dict_shape():
    prox = structural_proximity("S639L")
    required = {
        "token", "position", "window", "proximity_tier",
        "min_dist_to_drug_A", "sc_position", "sc_aa",
        "reference_structure", "note",
    }
    assert required <= set(prox.keys()), f"Missing keys: {required - set(prox.keys())}"
    assert prox["token"] == "S639L"
    assert prox["position"] == 639
    assert prox["window"] == "HS1"
    assert prox["proximity_tier"] in PROXIMITY_TIERS
    assert isinstance(prox["min_dist_to_drug_A"], float)
    assert prox["min_dist_to_drug_A"] > 0
    assert "9WZU" in prox["reference_structure"]
    assert isinstance(prox["note"], str) and len(prox["note"]) > 0


# ── window membership ─────────────────────────────────────────────────────────

def test_hs3_window_coverage():
    prox = structural_proximity("W691A")
    assert prox["window"] == "HS3"
    prox2 = structural_proximity("M690L")
    assert prox2["window"] == "HS3"


def test_hs2_window_coverage():
    prox = structural_proximity("R1354A")
    assert prox["window"] == "HS2"


def test_hs1_window_coverage():
    prox = structural_proximity("F635L")
    assert prox["window"] == "HS1"


# ── distal positions (panel residues far from drug) ───────────────────────────

@pytest.mark.parametrize("token", [
    "P643A",    # HS1 end: 13.14 Å → distal
    "Y689A",    # HS3 start: 13.07 Å → distal
])
def test_distal_positions(token):
    prox = structural_proximity(token)
    assert prox is not None
    assert prox["proximity_tier"] == DISTAL


# ── fks1_position integration: proximity key added for non-core tokens ─────────

def test_position_context_includes_proximity():
    from openafr.fks1_position import position_context
    ctx = position_context("S639L")
    assert ctx is not None
    assert "proximity" in ctx, "Rung 2 proximity missing from position_context output"
    assert ctx["proximity"]["proximity_tier"] == POCKET_ADJACENT


def test_position_context_no_proximity_for_core():
    from openafr.fks1_position import position_context
    # Core tokens return None from position_context itself
    assert position_context("S639F") is None


def test_position_context_proximity_consistent():
    from openafr.fks1_position import position_context
    ctx = position_context("W691A")
    assert ctx is not None
    assert "proximity" in ctx
    assert ctx["proximity"]["proximity_tier"] == DRUG_CONTACT
