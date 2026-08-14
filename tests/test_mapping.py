"""Tests for mutation -> pocket mapping (openafr/mapping.py, issue #23).

These lock the structural half of the early-warning loop -- the claims the alert
layer (#24/#25) will render:
  * a token is parsed into a residue substitution, or honestly refused,
  * a mapped residue's pocket-contact / iron distance is measured from the real
    modeled structure, not asserted,
  * the wild type at a position is ASSERTED (a numbering/structure mismatch is
    surfaced, never mapped through),
  * non-residue tokens (TR34) and absent residues are refused, not forced,
  * buildability tracks the SAME deterministic-deletion whitelist mutate_receptor
    enforces (so the two cannot silently drift apart),
  * the three canonical C. auris ERG11 mutations get the right structural verdict:
    Y132F/F126L contact the drug; K143R is second-shell.
No network: everything runs against the checked-in work/receptor_A.pdb + ref_VT1.pdb.
"""
import importlib.util
import os

import pytest

from openafr import mapping

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


def _map(token, struct):
    residues, ligand, fe = struct
    return mapping.map_mutation(token, residues, ligand, fe)


# ---- parse_substitution ------------------------------------------------------

def test_parse_canonical():
    assert mapping.parse_substitution("Y132F") == ("TYR", 132, "PHE")
    assert mapping.parse_substitution("k143r") == ("LYS", 143, "ARG")


def test_parse_stop_codon():
    assert mapping.parse_substitution("Y132*") == ("TYR", 132, "*")


def test_parse_rejects_non_substitution():
    for bad in ("TR34", "", "132", "YF", "Z10A", "Y132B", "TR34/L98H"):
        assert mapping.parse_substitution(bad) is None


# ---- honest refusals ---------------------------------------------------------

def test_non_residue_token_unmappable(struct):
    m = _map("TR34", struct)
    assert m["canonical"] is False and m["mappable"] is False
    assert "cannot map" in m["reason"]


def test_absent_residue_unmappable(struct):
    m = _map("Y9999F", struct)
    assert m["canonical"] is True and m["mappable"] is False
    assert "absent" in m["reason"]


def test_wildtype_mismatch_refused(struct):
    # position 132 is TYR in the model; a token claiming ALA there must be refused.
    m = _map("A132F", struct)
    assert m["mappable"] is False and m["wt_matches"] is False
    assert m["observed_wt"] == "TYR" and "mismatch" in m["reason"]


# ---- the three canonical C. auris ERG11 mutations ----------------------------

def test_Y132F_pocket_contact_and_buildable(struct):
    m = _map("Y132F", struct)
    assert m["mappable"] and m["wt_matches"]
    assert m["pocket_contact"] is True and m["dist_to_ligand"] < 4.5
    assert m["build_class"] == "deletion" and m["buildable"] is True


def test_F126L_pocket_contact_but_needs_repacker(struct):
    m = _map("F126L", struct)
    assert m["pocket_contact"] is True and m["dist_to_ligand"] < 4.5
    # Phe->Leu is not a pure deletion; it must be mapped but not buildable.
    assert m["build_class"] == "needs-repacker" and m["buildable"] is False


def test_K143R_second_shell(struct):
    m = _map("K143R", struct)
    assert m["mappable"] is True
    # K143 does not contact the drug -- the honest structural distinction.
    assert m["pocket_contact"] is False and m["dist_to_ligand"] > 4.5
    assert m["build_class"] == "needs-repacker"


# ---- buildability whitelist stays in sync with mutate_receptor ---------------

def test_deletion_whitelist_matches_mutate_receptor():
    """mapping.DELETION_MUTATIONS must equal mutate_receptor.SUPPORTED's key set.

    Loaded by file path (scripts/ is not a package) so a change to one whitelist
    that is not mirrored in the other fails here rather than drifting silently.
    """
    path = os.path.join(ROOT, "scripts", "mutate_receptor.py")
    spec = importlib.util.spec_from_file_location("mutate_receptor", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert set(mod.SUPPORTED.keys()) == mapping.DELETION_MUTATIONS


# ---- map_signals wrapper -----------------------------------------------------

def test_map_signals_accepts_signal_dicts_and_tokens():
    signals = [{"mutation": "Y132F"}, "K143R", {"mutation": "TR34"}]
    verdicts = mapping.map_signals(signals, receptor=RECEPTOR, ligand=LIGAND)
    assert [v["token"] for v in verdicts] == ["Y132F", "K143R", "TR34"]
    assert verdicts[0]["pocket_contact"] is True
    assert verdicts[2]["mappable"] is False
