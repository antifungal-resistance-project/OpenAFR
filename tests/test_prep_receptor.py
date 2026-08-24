"""Known-answer test for the receptor-prep extraction (scripts/prep_receptor.py).

The docking receptor work/receptor.pdbqt is built from the chain-A + heme atoms
pulled out of data/structures/5TZ1.pdb. What makes that receptor trustworthy is
that those atoms are the *same* ones the gate grades against, tracked in
work/receptor_A.pdb. This locks that invariant: the extraction must reproduce the
anchor exactly, and it must drop everything that is not chain-A protein or heme
(chain B, waters, the bound VT1 ligand). obabel is not needed — this tests the
selection logic, not the format conversion (that is verified by a VT-1161 redock).
"""
import pathlib

import pytest

from scripts.prep_receptor import atom_key, extract_chain_a, normalize_name_remark

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_RAW = _ROOT / "data" / "structures" / "5TZ1.pdb"
_ANCHOR = _ROOT / "work" / "receptor_A.pdb"

pytestmark = pytest.mark.skipif(
    not (_RAW.exists() and _ANCHOR.exists()),
    reason="needs data/structures/5TZ1.pdb and work/receptor_A.pdb",
)


def test_extraction_reproduces_the_anchor_exactly():
    kept = extract_chain_a(str(_RAW))
    ref = [l for l in open(_ANCHOR) if l[:6].strip() in ("ATOM", "HETATM")]
    assert len(kept) == len(ref)
    # Serial numbers get renumbered from 1 in the anchor; identity is everything else.
    assert [atom_key(l) for l in kept] == [atom_key(l) for l in ref]


def test_extraction_keeps_only_chain_a_protein_and_heme():
    kept = extract_chain_a(str(_RAW))
    assert all(l[21] == "A" for l in kept)                       # chain A only
    hetero = {l[17:20].strip() for l in kept if l[:6].strip() == "HETATM"}
    assert hetero == {"HEM"}                                     # heme kept, water/VT1 dropped
    assert sum(1 for l in kept if l[17:20].strip() == "HEM") == 43


def test_name_remark_is_normalized_so_the_receptor_can_be_hash_pinned(tmp_path):
    """obabel stamps its input filename into the first REMARK, and that input is a random
    temp file — which made two byte-identical receptors hash differently on every rebuild.
    The atom records must be untouched; only that one line is rewritten."""
    p = tmp_path / "receptor.pdbqt"
    body = "ATOM      1  N   LYS A   1      85.031  77.855  30.232  0.00  0.00    +0.000 N \n"
    p.write_text("REMARK  Name = /var/folders/xx/tmpABC123.pdb\nREMARK  keep me\n" + body)
    normalize_name_remark(str(p), "/some/where/data/structures/5TZ1.pdb")
    lines = p.read_text().splitlines(keepends=True)
    assert lines[0] == "REMARK  Name = 5TZ1.pdb\n"
    assert lines[1] == "REMARK  keep me\n"
    assert lines[2] == body


def test_name_remark_normalization_tolerates_a_file_without_one(tmp_path):
    p = tmp_path / "receptor.pdbqt"
    p.write_text("ATOM      1  N   LYS A   1       1.000   2.000   3.000  0.00  0.00 N \n")
    before = p.read_text()
    normalize_name_remark(str(p), "5TZ1.pdb")
    assert p.read_text() == before
