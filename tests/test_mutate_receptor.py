"""Known-answer tests for the receptor mutation step (scripts/mutate_receptor.py).

The *C. auris* resistance-pocket port grafts ERG11 resistance mutations onto the
validated 5TZ1 chain-A scaffold. The whole point is that it changes *only* the mutated
side chains and nothing else — same heme, same iron, same frame — so the frozen docking
box and gate still apply. These tests lock that:

  - a supported mutation (Y132F) is a pure atom deletion: the phenyl ring keeps its exact
    crystallographic coordinates, only the para-hydroxyl OH is dropped, and every other
    atom in the receptor is byte-identical to the wild type;
  - the heme iron never moves;
  - build-class mutations that add atoms (K143R, F126L) are REFUSED, not approximated —
    guessing a rotamer would break the same reproducibility contract as the frozen protocol;
  - a residue that is not the expected intact wild type is refused (assert, don't assume).
"""
import pathlib

import pytest

from scripts.mutate_receptor import (
    apply_mutation,
    atname,
    iron_xyz,
    parse_spec,
    read_atoms,
    resname,
    resseq,
)

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_ANCHOR = _ROOT / "work" / "receptor_A.pdb"

_needs_anchor = pytest.mark.skipif(
    not _ANCHOR.exists(), reason="needs work/receptor_A.pdb"
)


# --- pure spec parsing (no files) -------------------------------------------------

def test_parse_spec_reads_wildtype_position_and_mutant():
    assert parse_spec("Y132F") == ("TYR", 132, "PHE")
    assert parse_spec("v125a") == ("VAL", 125, "ALA")


@pytest.mark.parametrize("bad", ["Y132", "132F", "YF", "YXXF", "Z132F"])
def test_parse_spec_rejects_malformed_specs(bad):
    with pytest.raises(SystemExit):
        parse_spec(bad)


# --- deterministic deletion on the real anchor ------------------------------------

@_needs_anchor
def test_Y132F_deletes_only_the_hydroxyl_and_renames_to_phe():
    wt = read_atoms(str(_ANCHOR))
    mut = apply_mutation(wt, "TYR", 132, "PHE")

    wt132 = [l for l in wt if resseq(l) == "132"]
    mut132 = [l for l in mut if resseq(l) == "132"]
    assert {resname(l) for l in mut132} == {"PHE"}
    # exactly the OH atom is gone; nothing else about residue 132 changed
    assert {atname(l) for l in wt132} - {atname(l) for l in mut132} == {"OH"}
    # every retained atom keeps its exact coordinates (columns 30:54)
    kept_wt = {atname(l): l[30:54] for l in wt132 if atname(l) != "OH"}
    kept_mut = {atname(l): l[30:54] for l in mut132}
    assert kept_wt == kept_mut


@_needs_anchor
def test_mutation_leaves_every_other_atom_and_the_iron_untouched():
    wt = read_atoms(str(_ANCHOR))
    mut = apply_mutation(wt, "TYR", 132, "PHE")
    # all atoms outside residue 132 are byte-for-byte identical
    other_wt = [l for l in wt if resseq(l) != "132"]
    other_mut = [l for l in mut if resseq(l) != "132"]
    assert other_wt == other_mut
    assert iron_xyz(mut) == iron_xyz(wt) is not None


@_needs_anchor
def test_V125A_deletes_both_gamma_methyls():
    wt = read_atoms(str(_ANCHOR))
    mut = apply_mutation(wt, "VAL", 125, "ALA")
    mut125 = {atname(l) for l in mut if resseq(l) == "125"}
    assert mut125 == {"N", "CA", "C", "O", "CB"}


# --- refusals: honesty over approximation -----------------------------------------

@_needs_anchor
@pytest.mark.parametrize("spec", [("LYS", 143, "ARG"), ("PHE", 126, "LEU")])
def test_build_class_mutations_are_refused_not_approximated(spec):
    wt = read_atoms(str(_ANCHOR))
    with pytest.raises(SystemExit) as e:
        apply_mutation(wt, *spec)
    assert "repacker" in str(e.value)


@_needs_anchor
def test_wrong_wildtype_identity_is_refused():
    # residue 132 is TYR, not VAL; claiming V132A must fail loudly
    wt = read_atoms(str(_ANCHOR))
    with pytest.raises(SystemExit) as e:
        apply_mutation(wt, "VAL", 132, "ALA")
    assert "132" in str(e.value)
