"""Tests for the GNINA --score_only parser (openafr/gnina.py).

The metal-aware baseline comparison (work/PREREGISTRATION_baselines.md) ranks held-out
azoles by GNINA CNNaffinity on the same run-2 poses. Everything except the GNINA binary
itself runs on any host, so the parser is the load-bearing, testable piece: it must map
GNINA stdout to one dict per pose IN FILE ORDER (so it aligns 1:1 with pdbqt.read_poses),
pick the best CNNaffinity as higher-is-better, and honour the iron-bound `keep` mask used
for mode GC. Lock those behaviours.
"""
from openafr.gnina import best_cnnaffinity, parse_score_only

# Two poses as GNINA prints them under --score_only, with banner noise between blocks.
_TWO_POSE = """\
                     _             _
Using random seed: 42

## Ligand LIG
Affinity: -7.12 (kcal/mol)
CNNscore: 0.9123
CNNaffinity: 6.234
CNNvariance: 0.451

Affinity: -5.40 (kcal/mol)
CNNscore: 0.2010
CNNaffinity: 4.100
CNNvariance: 0.300
"""


def test_parse_preserves_pose_order_and_fields():
    poses = parse_score_only(_TWO_POSE)
    assert len(poses) == 2
    assert poses[0] == {"affinity": -7.12, "cnnscore": 0.9123, "cnnaffinity": 6.234}
    assert poses[1]["cnnaffinity"] == 4.100  # second block, in file order


def test_cnnaffinity_line_closes_a_pose():
    # A block missing Affinity/CNNscore still counts once CNNaffinity is seen.
    text = "CNNaffinity: 3.5\nCNNaffinity: 9.9\n"
    poses = parse_score_only(text)
    assert [p["cnnaffinity"] for p in poses] == [3.5, 9.9]
    assert poses[0]["affinity"] is None and poses[0]["cnnscore"] is None


def test_empty_output_is_no_poses():
    assert parse_score_only("no scores here\n") == []


def test_best_cnnaffinity_is_higher_is_better():
    poses = parse_score_only(_TWO_POSE)
    assert best_cnnaffinity(poses) == 6.234          # max, not min or first
    assert best_cnnaffinity([]) is None


def test_best_cnnaffinity_respects_keep_mask():
    poses = parse_score_only(_TWO_POSE)
    # Only the SECOND pose is iron-bound: the best eligible CNNaffinity is the lower value.
    assert best_cnnaffinity(poses, keep=[False, True]) == 4.100
    # No eligible pose -> None (molecule ranked LAST, never dropped, by the gate).
    assert best_cnnaffinity(poses, keep=[False, False]) is None
