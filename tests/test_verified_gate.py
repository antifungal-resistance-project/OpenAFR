"""Known-answer tests for the verified-inactive gate (scripts/validate_gate_verified.py).

This gate grades look #6 — the first on an independent dataset — so the two things it adds
on top of the run-2 grading math have to be right:

  * the tip statistic. On presumed decoys, six outranked the best true active. A raw count
    is not comparable between pools of different sizes, so the reported statistic is the
    best active's PERCENTILE rank; these lock the percentile, the rank, the count, and the
    all-decoys-at-the-top and no-rankable-active edge cases.
  * the pre-specified azole/non-azole subgroup split, which decides whether a FAIL reads as
    "geometry cannot separate actives from inactives" or the narrower and more likely
    "geometry cannot separate working azoles from failed azole analogues".
"""
import pytest

from scripts.validate_gate_verified import azole_bearing, tip_statistics


def _row(name, key, is_active):
    return (name, key, is_active)


def test_tip_statistics_on_a_clean_top():
    """Best active first: rank 1, nothing above it, top percentile."""
    ordered = [_row("act", 2.5, True)] + [_row(f"d{i}", 3.0 + i, False) for i in range(9)]
    rank, above, pct = tip_statistics(ordered)
    assert (rank, above) == (1, 0)
    assert pct == 10.0


def test_tip_statistics_counts_inactives_above_the_best_active():
    """The published failure mode: decoys crowding the tip."""
    ordered = [_row(f"d{i}", 2.0 + i * 0.1, False) for i in range(6)] + [_row("act", 2.7, True)]
    rank, above, pct = tip_statistics(ordered)
    assert (rank, above) == (7, 6)
    assert pct == 100.0


def test_tip_percentile_is_pool_size_neutral():
    """Six inactives above the best active means something different in a pool of 10 than in
    a pool of 1000 — which is exactly why the percentile is the reported statistic."""
    small = [_row(f"d{i}", i, False) for i in range(6)] + [_row("act", 9, True)] + \
            [_row(f"e{i}", 10 + i, False) for i in range(3)]
    big = [_row(f"d{i}", i, False) for i in range(6)] + [_row("act", 9, True)] + \
          [_row(f"e{i}", 10 + i, False) for i in range(93)]
    assert tip_statistics(small)[1] == tip_statistics(big)[1] == 6
    assert tip_statistics(small)[2] == 70.0
    assert tip_statistics(big)[2] == pytest.approx(7.0)


def test_tip_statistics_with_no_rankable_active():
    rank, above, pct = tip_statistics([_row("d1", 1.0, False), _row("d2", 2.0, False)])
    assert (rank, above, pct) == (None, None, None)


def test_azole_bearing_splits_the_warhead_from_the_rest(tmp_path):
    """Fluconazole-style triazoles and imidazoles are in; a plain benzamide is out."""
    p = tmp_path / "set.smi"
    p.write_text(
        "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F\ttriazole_one\n"       # 1,2,4-triazole
        "Clc1ccc(Cn2ccnc2)cc1\timidazole_one\n"                    # imidazole
        "O=C(Nc1ccccc1)c1ccccc1\tno_azole\n"                       # no azole ring
        "malformed line without a tab\n"
        "not-a-smiles\tunparseable\n")
    got = azole_bearing(str(p))
    assert got == {"triazole_one", "imidazole_one"}


def test_azole_bearing_on_the_real_actives():
    """Every held-out active carries the warhead — the subgroup split's own sanity check."""
    got = azole_bearing("data/ligands/actives_holdout_final.smi")
    assert len(got) == 7
