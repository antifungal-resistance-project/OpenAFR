"""Tests for the structural so-what estimator (openafr/structural.py, issue #24).

These lock the actionable tip of the early-warning loop -- the fit verdict the alert
layer (#25) renders. The bar they enforce is the honesty bar: every token resolves to
exactly one of measured | estimate | none, and the most dangerous output ("the drug
still fits") is never emitted without its caveats.

  * the deletion-atom whitelist stays in sync with mapping AND mutate_receptor,
  * a mutation with a real #13 dock is reported as MEASURED, carrying the pre-registered
    numbers and every caveat -- and never a bare "still works",
  * a deletion-buildable pocket contact whose removed atom touches the drug is a WEAKER-FIT
    estimate; a deletion that removes only second-shell atoms is minimal-direct-effect,
  * an unbuildable (add-atom) mutation is a NO-CALL -- plausible if pocket-contact, low
    prior if second-shell -- never an invented number,
  * a non-residue token (TR34) gets no structural fit call,
  * the geometric direction for Y132F agrees with its measured dock direction.
No network: everything runs against the checked-in work/receptor_A.pdb + ref_VT1.pdb.
"""
import importlib.util
import os

import pytest

from openafr import mapping, structural

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RECEPTOR = os.path.join(ROOT, "work", "receptor_A.pdb")
LIGAND = os.path.join(ROOT, "work", "ref_VT1.pdb")


@pytest.fixture(scope="module")
def struct():
    from openafr.pdbqt import iron_position
    residues = mapping.read_residues(RECEPTOR)
    ligand = mapping.read_ligand_atoms(LIGAND)
    fe = iron_position(RECEPTOR)
    return residues, ligand, fe


def _estimate(token, struct):
    residues, ligand, fe = struct
    mv = mapping.map_mutation(token, residues, ligand, fe)
    return structural.estimate_fit_consequence(mv, residues, ligand)


# ---- whitelist stays in sync (no silent drift across the three modules) ------

def test_deletion_atoms_match_mapping_and_mutate_receptor():
    """structural.DELETION_ATOMS keys == mapping.DELETION_MUTATIONS == mutate_receptor
    SUPPORTED, AND the atom sets equal mutate_receptor's. If any of the three changes
    without the others, this fails -- the estimate cannot measure a contact the builder
    would not actually remove."""
    assert set(structural.DELETION_ATOMS) == mapping.DELETION_MUTATIONS
    path = os.path.join(ROOT, "scripts", "mutate_receptor.py")
    spec = importlib.util.spec_from_file_location("mutate_receptor", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert set(mod.SUPPORTED) == set(structural.DELETION_ATOMS)
    assert mod.SUPPORTED == structural.DELETION_ATOMS


# ---- every branch resolves to exactly one confidence -------------------------

@pytest.mark.parametrize("token,evidence,confidence,direction", [
    ("Y132F", "empirical-dock",     "measured", "weaker-fit"),
    ("V125A", "geometric-estimate", "estimate", "minimal-direct-effect"),
    ("F126L", "no-call",            "none",     "uncertain"),
    ("K143R", "no-call",            "none",     "uncertain"),
    ("TR34",  "no-call",            "none",     "uncertain"),
])
def test_canonical_branches(token, evidence, confidence, direction, struct):
    v = _estimate(token, struct)
    assert v["evidence_class"] == evidence
    assert v["confidence"] == confidence
    assert v["direction"] == direction
    assert v["confidence"] in ("measured", "estimate", "none")


# ---- Tier 1: the measured Y132F verdict is complete and honest ---------------

def test_y132f_reports_measured_dock_with_caveats(struct):
    v = _estimate("Y132F", struct)
    assert v["empirical"] is not None
    assert v["empirical"]["gate"] == "FAIL"
    # the pre-registered numbers are surfaced, not summarised away
    assert v["empirical"]["ef1_wildtype"] == 12.68 and v["empirical"]["ef1_mutant"] == 0.00
    # the most dangerous claim never ships bare: the class caveat is always present
    assert any("NOT fluconazole specifically" in c for c in v["caveats"])
    assert len(v["caveats"]) >= 3
    # and it says azoles still reach the iron rather than "the drug is dead"
    assert "not structurally *excluded*" in v["fluconazole_fit_verdict"]


def test_y132f_geometry_agrees_with_dock_direction(struct):
    """The deterministic geometric estimate (removed para-OH contacts the drug) points
    the SAME way as the measured dock (weaker fit) -- an independent corroboration, so a
    later dock-registry edit that flipped the direction would strand this cross-check."""
    residues, ligand, fe = struct
    mv = mapping.map_mutation("Y132F", residues, ligand, fe)
    lost = structural.lost_pocket_contacts(mv, residues, ligand)
    assert [c["atom"] for c in lost] == ["OH"]
    assert lost[0]["dist_to_drug"] <= mapping.CONTACT_CUTOFF
    assert structural.EMPIRICAL_DOCKS["Y132F"]["direction"] == "weaker-fit"


# ---- Tier 2: buildable estimates are directional, never numeric --------------

def test_v125a_lost_contacts_are_second_shell(struct):
    """V125A is deletion-buildable but its removed methyls are all beyond the contact
    cutoff -> no lost drug contact -> minimal-direct-effect, not weaker-fit."""
    residues, ligand, fe = struct
    mv = mapping.map_mutation("V125A", residues, ligand, fe)
    assert structural.lost_pocket_contacts(mv, residues, ligand) == []
    v = structural.estimate_fit_consequence(mv, residues, ligand)
    assert v["direction"] == "minimal-direct-effect"
    assert "second-shell" in v["fluconazole_fit_verdict"]


def test_estimate_never_emits_a_number(struct):
    """A non-empirical verdict must not carry a fabricated magnitude: no empirical block,
    and the reproduce command (dock to measure) is offered instead."""
    v = _estimate("V125A", struct)
    assert v["empirical"] is None
    assert "mutate_receptor.py" in v["fluconazole_fit_verdict"]


# ---- Tier 3: unbuildable and unmappable are honest no-calls ------------------

def test_f126l_pocket_contact_unbuildable_is_plausible_nocall(struct):
    v = _estimate("F126L", struct)
    assert v["pocket_contact"] is True and v["buildable"] is False
    assert "plausible" in v["fluconazole_fit_verdict"].lower()
    assert "repacker" in v["fluconazole_fit_verdict"].lower()


def test_k143r_second_shell_unbuildable_is_low_prior_nocall(struct):
    v = _estimate("K143R", struct)
    assert v["pocket_contact"] is False and v["buildable"] is False
    assert "low prior" in v["fluconazole_fit_verdict"].lower()


def test_tr34_gets_no_structural_call(struct):
    v = _estimate("TR34", struct)
    assert v["mappable"] is False
    assert v["evidence_class"] == "no-call"
    assert "no structural fit call" in v["fluconazole_fit_verdict"]


# ---- the loop: emergence signals map + estimate in one call ------------------

def test_estimate_signals_consumes_emergence_shape():
    """estimate_signals accepts emergence signal dicts (a 'mutation' key) and bare tokens
    alike, returning one verdict per input in order -- the detect -> map -> estimate loop."""
    signals = [{"mutation": "Y132F"}, "F126L", {"mutation": "TR34"}]
    out = structural.estimate_signals(signals, receptor=RECEPTOR, ligand=LIGAND)
    assert [v["token"] for v in out] == ["Y132F", "F126L", "TR34"]
    assert out[0]["confidence"] == "measured"
    assert out[2]["evidence_class"] == "no-call"
