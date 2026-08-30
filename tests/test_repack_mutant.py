"""Known-answer tests for the add-atom side-chain repacker (scripts/repack_mutant.py, #77).

The repacker deliberately crosses the "no modeled geometry" line that mutate_receptor.py holds,
so its guard rails matter: it must own EXACTLY the two add-atom mutants Y132F-style deletion
cannot build (K143R, F126L), refuse everything else — including the deletion-class specs that
belong to mutate_receptor.py — and never conflate a modeled rotamer with an exact one. These
tests pin the pure spec/whitelist logic; they do not require PyMOL (the wizard call is imported
locally inside repack()). The determinism + isolation guarantees are exercised by the build
step itself (repack_mutant.verify_isolated refuses a non-isolated or iron-moving build).
"""
import pytest

from scripts import repack_mutant as rm


def test_parses_the_two_owned_addatom_mutants():
    assert rm.parse_spec("K143R") == ("LYS", 143, "ARG")
    assert rm.parse_spec("F126L") == ("PHE", 126, "LEU")
    assert rm.parse_spec("k143r") == ("LYS", 143, "ARG")   # case-insensitive


def test_refuses_deletion_class_specs_that_belong_to_mutate_receptor():
    # Y132F / V125A are pure deletions — mutate_receptor.py builds them with no rotamer at all;
    # the repacker must NOT quietly take them (that would model geometry that is exactly known).
    for spec in ("Y132F", "V125A"):
        with pytest.raises(SystemExit):
            rm.parse_spec(spec)


def test_refuses_unknown_or_malformed_specs():
    for spec in ("K143", "143R", "Z143R", "K143Z", "", "KR"):
        with pytest.raises(SystemExit):
            rm.parse_spec(spec)


def test_supported_set_is_exactly_the_two_addatom_pairs():
    assert rm.SUPPORTED == {("LYS", "ARG"), ("PHE", "LEU")}


def test_heavy_atom_sets_have_the_right_cardinality():
    # These counts are what the post-repack verification asserts the residue came out with.
    # Lys->Arg: drop CE,NZ, add NE,CZ,NH1,NH2 => net +2 heavy atoms.
    assert len(rm.HEAVY_ATOMS["ARG"]) == len(rm.HEAVY_ATOMS["LYS"]) + 2
    # Phe->Leu: keep CG,CD1,CD2, drop the CE1,CE2,CZ ring atoms => net -3 heavy atoms.
    assert len(rm.HEAVY_ATOMS["LEU"]) == len(rm.HEAVY_ATOMS["PHE"]) - 3
    # backbone atoms present in every residue
    for r in ("LYS", "ARG", "PHE", "LEU"):
        assert {"N", "CA", "C", "O", "CB"} <= rm.HEAVY_ATOMS[r]
