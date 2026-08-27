"""Known-answer tests for the resistance-retention classifier
(scripts/candidate_resistance.py, issues #69/#77).

These lock the per-candidate verdict rule frozen in
work/PREREGISTRATION_candidate_resistance.md, BEFORE any Y132F distance was read:
  - LOST  iff the molecule coordinates the Y132F iron in < MIN_COORD_SEEDS seeds
          (the informative demotion; checked before the shift, since a shift is only
          meaningful once the molecule reliably reaches the iron);
  - RETAINED  iff it coordinates in >= MIN_COORD_SEEDS AND |shift| <= RETENTION_TOL;
  - SHIFTED-away / SHIFTED-closer  iff it coordinates but moves by more than the tolerance,
          with the direction reported (never silently praising a molecule that moved closer).
"""
from scripts.candidate_resistance import MIN_COORD_SEEDS, RETENTION_TOL, classify


def test_reliable_and_unmoved_is_retained():
    # lersivirine's real shape: 5/5, +0.064 A shift
    verdict, shift = classify(y132f_cons=2.636, y132f_ncoord=5, wt_cons=2.572)
    assert verdict == "RETAINED"
    assert round(shift, 3) == 0.064


def test_too_few_coordinating_seeds_is_lost():
    # vatalanib's real shape: only 2/5 seeds coordinate on Y132F -> LOST regardless of the min
    verdict, shift = classify(y132f_cons=2.510, y132f_ncoord=2, wt_cons=2.479)
    assert verdict == "LOST"
    assert shift is None  # no shift is claimed when the molecule does not reliably coordinate


def test_lost_takes_precedence_over_a_small_shift():
    # even a near-zero shift cannot rescue a molecule that fails the seed count
    verdict, _ = classify(y132f_cons=2.50, y132f_ncoord=MIN_COORD_SEEDS - 1, wt_cons=2.50)
    assert verdict == "LOST"


def test_large_move_away_is_shifted_away():
    verdict, shift = classify(y132f_cons=3.00 - 0.01, y132f_ncoord=5, wt_cons=2.40)
    # 2.99 - 2.40 = 0.59 > 0.5 tol
    assert verdict == "SHIFTED-away"
    assert shift > RETENTION_TOL


def test_large_move_closer_is_shifted_closer_not_silently_praised():
    verdict, shift = classify(y132f_cons=2.00, y132f_ncoord=5, wt_cons=2.60)
    assert verdict == "SHIFTED-closer"
    assert shift < -RETENTION_TOL


def test_boundary_shift_at_tolerance_is_retained():
    # exactly RETENTION_TOL is within tolerance (<=), so still RETAINED
    verdict, _ = classify(y132f_cons=2.50 + RETENTION_TOL, y132f_ncoord=5, wt_cons=2.50)
    assert verdict == "RETAINED"
