"""Tests for cross-species CYP51 pocket conservation (openafr/conservation.py, issue #76).

These lock the claims the pocket-conservation map makes:
  * translation of the pinned C. auris ERG11 CDS matches its committed protein length
    and carries V125/F126/Y132/K143 at the C. albicans-frame positions,
  * the Needleman-Wunsch aligner is deterministic and puts the canonical resistance
    residues in register across orthologs (a self-identity and a landmark check),
  * the conservation classifier applies the frozen rule (IDENTICAL / CONSERVATIVE by
    BLOSUM62 >= 1 / DIVERGENT on gap or non-positive substitution),
  * the numbering is ASSERTED -- a structure/reference frame mismatch raises, never
    aligns the wrong columns,
  * the built map reproduces the pre-registered headline: the pocket is 22 residues,
    invariant C. albicans<->C. auris, with divergence only toward Aspergillus and one
    such divergence (G307) in the coordination shell.
No network: everything runs against checked-in structures + pinned sequences.
"""
import os

import pytest

from openafr import conservation as cz

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RECEPTOR = os.path.join(ROOT, "work", "receptor_A.pdb")
LIGAND = os.path.join(ROOT, "work", "ref_VT1.pdb")
DATA = os.path.join(ROOT, "data", "cyp51_orthologs")
AURIS_CDS = os.path.join(ROOT, "data", "earlywarning", "erg11_reference", "erg11_cds.fasta")


@pytest.fixture(scope="module")
def seqs():
    alb = cz.read_fasta(os.path.join(DATA, "P10613.fasta"))[1]
    aspA = cz.read_fasta(os.path.join(DATA, "Q4WNT5.fasta"))[1]
    aspB = cz.read_fasta(os.path.join(DATA, "Q96W81.fasta"))[1]
    auris = cz.translate(cz.read_fasta(AURIS_CDS)[1])
    return dict(alb=alb, aspA=aspA, aspB=aspB, auris=auris)


@pytest.fixture(scope="module")
def result(seqs):
    return cz.build_map(
        RECEPTOR, LIGAND, seqs["alb"],
        verdict_orthologs={"auris": seqs["auris"], "aspergillus_A": seqs["aspA"]},
        supplementary={"aspergillus_B": seqs["aspB"]},
    )


# --- translation ---------------------------------------------------------------

def test_auris_translation_length_and_landmarks(seqs):
    auris = seqs["auris"]
    assert len(auris) == 524                     # 1575-nt CDS -> 524 aa + dropped stop
    assert not auris.endswith("*")
    # canonical resistance residues, C. albicans-frame positions
    assert (auris[124], auris[125], auris[131], auris[142]) == ("V", "F", "Y", "K")


def test_translate_drops_single_trailing_stop():
    assert cz.translate("ATGAAATAA") == "MK"     # M K *  -> "MK"
    assert cz.translate("ATGAAA") == "MK"         # no stop -> unchanged


# --- aligner -------------------------------------------------------------------

def test_alignment_self_identity(seqs):
    assert cz.percent_identity(seqs["alb"], seqs["alb"]) == 100.0


def test_alignment_is_deterministic(seqs):
    a = cz.needleman_wunsch(seqs["alb"], seqs["aspA"])
    b = cz.needleman_wunsch(seqs["alb"], seqs["aspA"])
    assert a == b


def test_alignment_registers_canonical_residues(seqs):
    # In the C. albicans frame, position 132 is Y in every ortholog (a hard landmark).
    for other in ("auris", "aspA"):
        col = cz.reference_column_map(seqs["alb"], seqs[other])
        assert col[132] == "Y", other


# --- classifier ----------------------------------------------------------------

def test_classify_identical():
    assert cz.classify("Y", {"a": "Y", "b": "Y"}) == "IDENTICAL"


def test_classify_conservative_positive_blosum():
    # I<->V scores +3 in BLOSUM62 -> conservative.
    assert cz.classify("I", {"a": "I", "b": "V"}) == "CONSERVATIVE"


def test_classify_divergent_on_negative_blosum():
    # G<->A scores 0 -> not positive -> divergent.
    assert cz.classify("G", {"a": "G", "b": "A"}) == "DIVERGENT"


def test_classify_divergent_on_gap():
    assert cz.classify("Y", {"a": "Y", "b": "-"}) == "DIVERGENT"


# --- numbering assertion (honesty constraint 1) --------------------------------

def test_assert_reference_frame_raises_on_mismatch():
    fake_pocket = [{"pos": 132, "resname": "TYR", "aa1": "Y"}]
    with pytest.raises(ValueError):
        cz.assert_reference_frame(fake_pocket, "A" * 200)   # pos 132 is 'A', not 'Y'


def test_real_structure_passes_reference_frame(seqs, result):
    # build_map already asserts internally; reaching here means the 22 pocket
    # residues all matched P10613. Re-assert explicitly for clarity.
    cz.assert_reference_frame(result["pocket"], seqs["alb"])


# --- the built map reproduces the pre-registered headline ----------------------

def test_pocket_size_and_tally(result):
    s = result["summary"]
    assert s["n_pocket"] == 22
    assert (s["n_identical"], s["n_conservative"], s["n_divergent"]) == (14, 5, 3)


def test_auris_pocket_is_near_invariant(result):
    # The C. albicans -> C. auris pocket differs at exactly one lining residue,
    # I304V, and that lone difference is conservative (never DIVERGENT).
    diffs = [r for r in result["pocket"] if r["orthologs"]["auris"] != r["ref_aa"]]
    assert [r["pos"] for r in diffs] == [304]
    (only,) = diffs
    assert only["ref_aa"] == "I" and only["orthologs"]["auris"] == "V"
    assert only["verdict"] == "CONSERVATIVE"


def test_divergence_is_only_toward_aspergillus(result):
    assert result["summary"]["divergent_positions"] == [303, 307, 380]
    for r in result["pocket"]:
        if r["verdict"] == "DIVERGENT":
            # the divergence never comes from C. auris
            assert r["orthologs"]["auris"] == r["ref_aa"]


def test_coordination_shell_divergence(result):
    assert result["summary"]["coord_shell_divergent"] == [307]
