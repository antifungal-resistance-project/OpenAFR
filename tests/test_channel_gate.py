"""Known-answer tests for the channel-engagement criterion (scripts/validate_gate_channel.py).

mode F ranks within the azole class on a signal orthogonal to N-Fe distance: how many
channel-lining residues the tail contacts, measured in the coordinating pose. These lock the
pieces the PASS/FAIL result would rest on, so a bug cannot silently corrupt the one graded look:
  * the channel-shell is LABEL-BLIND -- defined from the co-crystal ligand and the receptor
    alone, never from an active/inactive label -- and excludes the heme (a HETATM group);
  * a residue is counted at most ONCE (distinct-residue count, less size-sensitive than atoms),
    and only when a ligand heavy atom is within the fixed 4.5 A cutoff;
  * engagement is read in the COORDINATING pose -- the same min-N-Fe pose mode C selects -- so
    there is no min-vs-max-over-poses degree of freedom;
  * higher engagement ranks FIRST, and a molecule with no coordinating pose is unrankable (last).
No network, no docking; every case is synthetic.
"""
import math

import pytest

from scripts.validate_gate_channel import (
    _pearson, _rank_key, _residualize, channel_shell, coordinating_pose, engagement,
)


def _rec(record, chars):
    """One fixed-column PDB record; assign only the fields the parsers read."""
    s = list(" " * 80)
    s[0:6] = list(f"{record:<6}")
    for a, b, v in chars:
        s[a:b] = list(v)
    return "".join(s)


def _atom(record, name, resname, pos, x, y, z, elem):
    return _rec(record, [(12, 16, f"{name:<4}"), (17, 20, f"{resname:<3}"),
                         (22, 26, f"{pos:>4}"), (30, 38, f"{x:8.3f}"),
                         (38, 46, f"{y:8.3f}"), (46, 54, f"{z:8.3f}"), (76, 78, f"{elem:>2}")])


def _pose_file(tmp_path, models):
    lines = []
    for i, atoms in enumerate(models, 1):
        lines.append(f"MODEL {i}")
        lines += atoms
        lines.append("ENDMDL")
    p = tmp_path / "pose.pdbqt"
    p.write_text("\n".join(lines) + "\n")
    return str(p)


# --- channel_shell: label-blind, heme-excluded ---------------------------------

def test_channel_shell_is_label_blind_and_excludes_heme(tmp_path):
    # A protein residue 3 A from the ref ligand -> in the shell. A far residue -> out.
    # A heme HETATM sitting right on the ligand -> never a shell residue (not a protein ATOM).
    rec = tmp_path / "rec.pdb"
    rec.write_text("\n".join([
        _atom("ATOM", "CA", "PHE", 126, 0.0, 0.0, 3.0, "C"),   # 3.0 A from ligand at origin
        _atom("ATOM", "CA", "ALA", 999, 0.0, 0.0, 40.0, "C"),  # far away
        _atom("HETATM", "FE", "HEM", 601, 0.0, 0.0, 0.0, "FE"),  # on the ligand, but HETATM
    ]) + "\n")
    lig = tmp_path / "lig.pdb"
    lig.write_text(_atom("HETATM", "N1", "VT1", 602, 0.0, 0.0, 0.0, "N") + "\n")

    shell = channel_shell(str(rec), str(lig), cutoff=4.5)
    assert set(shell) == {126}          # only the near protein residue; not 999, not the heme


def test_channel_shell_respects_the_cutoff(tmp_path):
    rec = tmp_path / "rec.pdb"
    rec.write_text("\n".join([
        _atom("ATOM", "CA", "LEU", 5, 0.0, 0.0, 4.4, "C"),   # inside 4.5
        _atom("ATOM", "CA", "LEU", 6, 0.0, 0.0, 4.6, "C"),   # outside 4.5
    ]) + "\n")
    lig = tmp_path / "lig.pdb"
    lig.write_text(_atom("HETATM", "N1", "VT1", 602, 0.0, 0.0, 0.0, "N") + "\n")
    assert set(channel_shell(str(rec), str(lig), cutoff=4.5)) == {5}


# --- engagement: distinct residues, once each, within cutoff -------------------

def _shell(**residues):
    """residues: pos=name->coord list, as {pos: [(atom_name, (x,y,z)), ...]}."""
    return {pos: [(nm, xyz) for nm, xyz in atoms] for pos, atoms in residues.items()}


def test_engagement_counts_distinct_residues_once():
    shell = _shell(
        **{"10": [("CA", (0.0, 0.0, 0.0))],
           "11": [("CA", (10.0, 0.0, 0.0))],
           "12": [("CA", (20.0, 0.0, 0.0))]})
    # Two ligand atoms both near residue 10, one near residue 11, none near 12.
    pose = [("C1", 1.0, 0.0, 0.0), ("C2", 1.5, 0.0, 0.0), ("C3", 10.5, 0.0, 0.0)]
    assert engagement(pose, shell, cutoff=4.5) == 2      # residues 10 and 11, not double-counted


def test_engagement_respects_cutoff():
    shell = _shell(**{"10": [("CA", (0.0, 0.0, 0.0))]})
    assert engagement([("C1", 4.4, 0.0, 0.0)], shell, cutoff=4.5) == 1
    assert engagement([("C1", 4.6, 0.0, 0.0)], shell, cutoff=4.5) == 0


# --- coordinating_pose: the same pose mode C selects ---------------------------

FE = (0.0, 0.0, 0.0)


def test_coordinating_pose_is_the_min_nfe_pose(tmp_path):
    far = [_atom("HETATM", "N1", "LIG", 1, 0.0, 0.0, 5.0, "N"),
           _atom("HETATM", "C1", "LIG", 1, 1.0, 0.0, 5.0, "C")]
    near = [_atom("HETATM", "N1", "LIG", 1, 0.0, 0.0, 2.0, "N"),   # closer N -> this pose wins
            _atom("HETATM", "C1", "LIG", 1, 9.0, 9.0, 9.0, "C")]
    p = _pose_file(tmp_path, [far, near])
    atoms, nfe = coordinating_pose(p, FE)
    assert math.isclose(nfe, 2.0, abs_tol=1e-6)
    assert any(math.isclose(z, 9.0, abs_tol=1e-6) for (_n, _x, _y, z) in atoms)  # the near model


def test_coordinating_pose_none_when_no_nitrogen(tmp_path):
    p = _pose_file(tmp_path, [[_atom("HETATM", "C1", "LIG", 1, 0.0, 0.0, 2.0, "C")]])
    atoms, nfe = coordinating_pose(p, FE)
    assert atoms is None and nfe is None


# --- ranking + diagnostics -----------------------------------------------------

def test_rank_key_ranks_higher_engagement_first_and_unrankable_last():
    rows = [("hi", 8, True, 30, 2.5), ("lo", 2, False, 20, 2.6), ("none", None, False, None, None)]
    keyed = _rank_key(rows)
    assert keyed[0] == ("hi", -8, True) and keyed[1] == ("lo", -2, False)
    assert keyed[2] == ("none", None, False)          # None key -> report places it last
    ranked = sorted([k for k in keyed if k[1] is not None], key=lambda k: k[1])
    assert [k[0] for k in ranked] == ["hi", "lo"]     # higher engagement first


def test_pearson_and_residualize_known_answers():
    assert math.isclose(_pearson([1, 2, 3], [2, 4, 6]), 1.0, abs_tol=1e-9)
    assert math.isclose(_pearson([1, 2, 3], [6, 4, 2]), -1.0, abs_tol=1e-9)
    # feature perfectly linear in size -> residuals all ~0 (no signal beyond size)
    resid = _residualize([2.0, 4.0, 6.0], [1.0, 2.0, 3.0])
    assert all(abs(r) < 1e-9 for r in resid)
