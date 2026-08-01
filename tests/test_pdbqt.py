"""Tests for the shared Vina/PDBQT parsers (openafr/pdbqt.py, task T1).

These four functions were copied by hand into four scripts before T1 consolidated
them. They read docked geometry that the whole result depends on, so lock the
behaviours that matter: hydrogens excluded, best-first score order preserved,
missing iron refused loudly, and no-nitrogen poses reported as unrankable (None).
"""
import pytest

from openafr import pdbqt

_RECEPTOR = (
    "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
    "HETATM 9999 FE   HEM A 500      10.000  10.000  10.000  1.00  0.00          FE\n"
)

_TWO_POSE_LIGAND = (
    "MODEL 1\n"
    "REMARK VINA RESULT:    -9.9      0.000      0.000\n"
    "ATOM      1  N   LIG A   1      11.000  10.000  10.000  1.00  0.00           N\n"
    "HETATM    2  H   LIG A   1      11.500  10.000  10.000  1.00  0.00           H\n"
    "ATOM      3  C   LIG A   1      20.000  20.000  20.000  1.00  0.00           C\n"
    "ENDMDL\n"
    "MODEL 2\n"
    "REMARK VINA RESULT:    -8.1      1.500      2.500\n"
    "ATOM      4  C   LIG A   1      30.000  30.000  30.000  1.00  0.00           C\n"
    "ENDMDL\n"
)


@pytest.fixture
def receptor(tmp_path):
    p = tmp_path / "rec.pdb"
    p.write_text(_RECEPTOR)
    return str(p)


@pytest.fixture
def ligand(tmp_path):
    p = tmp_path / "lig.pdbqt"
    p.write_text(_TWO_POSE_LIGAND)
    return str(p)


def test_iron_position(receptor):
    assert pdbqt.iron_position(receptor) == (10.0, 10.0, 10.0)


def test_iron_position_missing_raises(tmp_path):
    empty = tmp_path / "no_iron.pdb"
    empty.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n")
    with pytest.raises(SystemExit):
        pdbqt.iron_position(str(empty))


def test_read_scores_in_file_order(ligand):
    assert pdbqt.read_scores(ligand) == [-9.9, -8.1]


def test_read_poses_excludes_hydrogens(ligand):
    poses = list(pdbqt.read_poses(ligand))
    assert [num for num, _ in poses] == [1, 2]
    names_pose1 = [a[0] for a in poses[0][1]]
    assert names_pose1 == ["N", "C"]           # the H atom is dropped


def test_min_nitrogen_iron_distance(ligand, receptor):
    fe = pdbqt.iron_position(receptor)
    poses = list(pdbqt.read_poses(ligand))
    # pose 1 has an N at (11,10,10), 1.0 A from the iron
    assert pdbqt.min_nitrogen_iron_distance(poses[0][1], fe) == pytest.approx(1.0)
    # pose 2 has no nitrogen -> unrankable
    assert pdbqt.min_nitrogen_iron_distance(poses[1][1], fe) is None
