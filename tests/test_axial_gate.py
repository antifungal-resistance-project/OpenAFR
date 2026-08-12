"""Known-answer tests for the axial-geometry criterion (scripts/validate_gate_axial.py).

The axial gate adds a signal orthogonal to N-Fe distance: the angle theta between the
coordinating nitrogen's approach vector and the heme normal, combined as S = |v|*(1+sin theta).
These lock the geometry that the FAIL result rests on — that the normal is the unit porphyrin
perpendicular, that theta folds to the acute value (so the normal's sign cannot be tuned), that
S = |v| for a perfectly axial approach and 2|v| for a fully equatorial one, that theta is
measured on the nitrogen nearest the iron, and that the per-seed minima aggregate by a
None-dropping median.
"""
import math

import pytest

from scripts.validate_gate_axial import _median, heme_normal, per_seed_minima, pose_axial


def _rec(chars):
    """Build one fixed-column PDB record; assign only the fields the parsers read."""
    s = list(" " * 80)
    s[0:6] = list("HETATM")
    for a, b, v in chars:
        s[a:b] = list(v)
    return "".join(s)


def _hem_atom(name, x, y, z):
    return _rec([(12, 16, f"{name:<4}"), (17, 20, "HEM"), (21, 22, "A"),
                 (30, 38, f"{x:8.3f}"), (38, 46, f"{y:8.3f}"), (46, 54, f"{z:8.3f}"),
                 (76, 78, f"{'N':>2}")])


def _n_atom(x, y, z, name="N"):
    return _rec([(12, 16, f"{name:<4}"), (17, 20, "LIG"), (21, 22, "A"),
                 (30, 38, f"{x:8.3f}"), (38, 46, f"{y:8.3f}"), (46, 54, f"{z:8.3f}"),
                 (76, 78, f"{'N':>2}")])


def _pose_file(tmp_path, models):
    """models: list of lists of atom records -> one multi-MODEL pdbqt file."""
    lines = []
    for i, atoms in enumerate(models, 1):
        lines.append(f"MODEL {i}")
        lines += atoms
        lines.append("ENDMDL")
    p = tmp_path / "pose.pdbqt"
    p.write_text("\n".join(lines) + "\n")
    return str(p)


# --- heme normal ---------------------------------------------------------------

def test_heme_normal_is_unit_perpendicular_to_porphyrin_plane(tmp_path):
    # pyrrole Ns on the xy-plane -> normal is +/- z, and always a unit vector
    rec = tmp_path / "rec.pdb"
    rec.write_text("\n".join([
        _hem_atom("NA", 1.0, 0.0, 0.0), _hem_atom("NC", -1.0, 0.0, 0.0),
        _hem_atom("NB", 0.0, 1.0, 0.0), _hem_atom("ND", 0.0, -1.0, 0.0),
    ]) + "\n")
    n = heme_normal(str(rec))
    assert math.isclose(math.sqrt(sum(c * c for c in n)), 1.0, abs_tol=1e-9)
    assert math.isclose(abs(n[2]), 1.0, abs_tol=1e-9)   # perpendicular to xy-plane
    assert math.isclose(n[0], 0.0, abs_tol=1e-9) and math.isclose(n[1], 0.0, abs_tol=1e-9)


def test_heme_normal_missing_atom_raises(tmp_path):
    rec = tmp_path / "rec.pdb"
    rec.write_text(_hem_atom("NA", 1.0, 0.0, 0.0) + "\n")  # only one pyrrole N
    with pytest.raises(SystemExit):
        heme_normal(str(rec))


# --- axial angle + S -----------------------------------------------------------

FE = (0.0, 0.0, 0.0)
NHAT = (0.0, 0.0, 1.0)


def test_axial_nitrogen_gives_theta_zero_and_S_equals_distance(tmp_path):
    p = _pose_file(tmp_path, [[_n_atom(0.0, 0.0, 2.2)]])       # straight up the normal
    (d, theta, S), = pose_axial(p, FE, NHAT)
    assert math.isclose(d, 2.2, abs_tol=1e-6)
    assert math.isclose(theta, 0.0, abs_tol=1e-6)
    assert math.isclose(S, 2.2, abs_tol=1e-6)                  # |v|*(1+sin 0)


def test_equatorial_nitrogen_gives_theta_ninety_and_S_doubles(tmp_path):
    p = _pose_file(tmp_path, [[_n_atom(2.2, 0.0, 0.0)]])       # in-plane, off-axis
    (d, theta, S), = pose_axial(p, FE, NHAT)
    assert math.isclose(theta, 90.0, abs_tol=1e-6)
    assert math.isclose(S, 4.4, abs_tol=1e-6)                  # |v|*(1+sin 90)


def test_normal_sign_is_irrelevant(tmp_path):
    # folding to the acute angle means flipping n leaves S unchanged (nothing to tune)
    p = _pose_file(tmp_path, [[_n_atom(0.5, 0.0, 2.0)]])
    up = pose_axial(p, FE, NHAT)[0]
    down = pose_axial(p, FE, (0.0, 0.0, -1.0))[0]
    assert math.isclose(up[1], down[1], abs_tol=1e-9)
    assert math.isclose(up[2], down[2], abs_tol=1e-9)


def test_theta_measured_on_nitrogen_nearest_iron(tmp_path):
    # a far axial N and a near equatorial N -> the NEAR one is scored (theta ~ 90)
    p = _pose_file(tmp_path, [[_n_atom(0.0, 0.0, 9.0), _n_atom(2.0, 0.0, 0.0)]])
    (d, theta, S), = pose_axial(p, FE, NHAT)
    assert math.isclose(d, 2.0, abs_tol=1e-6)                  # nearest N, not the far one
    assert math.isclose(theta, 90.0, abs_tol=1e-6)


# --- per-seed minima + aggregation --------------------------------------------

def test_per_seed_minima_take_each_column_independently(tmp_path):
    # pose A: close but equatorial; pose B: far but axial. min S, min theta, min dist
    # can each come from a different pose, per the pre-registration.
    p = _pose_file(tmp_path, [
        [_n_atom(2.0, 0.0, 0.0)],     # d=2.0  theta=90  S=4.0
        [_n_atom(0.0, 0.0, 3.0)],     # d=3.0  theta=0   S=3.0
    ])
    minS, minTheta, minDist = per_seed_minima(p, FE, NHAT)
    assert math.isclose(minS, 3.0, abs_tol=1e-6)              # from pose B
    assert math.isclose(minTheta, 0.0, abs_tol=1e-6)         # from pose B
    assert math.isclose(minDist, 2.0, abs_tol=1e-6)          # from pose A


def test_pose_with_no_nitrogen_is_unrankable(tmp_path):
    p = tmp_path / "empty.pdbqt"
    p.write_text("MODEL 1\nENDMDL\n")
    assert per_seed_minima(str(p), FE, NHAT) == (None, None, None)


def test_median_drops_none_and_takes_the_middle():
    assert _median({"42": 3.1, "1": 3.0, "2": None, "3": 3.2}) == 3.1  # median of [3.0,3.1,3.2]
    assert _median({"42": None, "1": None}) is None                    # -> ranked LAST
