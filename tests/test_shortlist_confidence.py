"""Known-answer tests for the shortlist-confidence logic (scripts/shortlist_confidence.py):
seed stability (#68) and human-CYP51 selectivity (#73).

These lock the pure decision rules frozen in the two pre-registrations:
  - consensus = MIN over seeds that COORDINATE (N-Fe < FE_CUTOFF); a seed above the cutoff is
    not a coordinating seed and must not count toward n_coord;
  - a molecule that coordinates in no seed reports n_coord=0 (falls back to the closest
    approach seen, so it stays informative but is never called "reaching the iron");
  - STABLE iff it coordinates in EVERY seed AND spread <= SPREAD_MAX;
  - selectivity margin = human - fungal, POSITIVE meaning the fungal iron is reached closer.
"""
from scripts.shortlist_confidence import FE_CUTOFF, SPREAD_MAX, assess, consensus


def test_consensus_min_over_coordinating_seeds_only():
    # 3.70 is above the 3.0 cutoff -> not a coordinating seed, excluded from min and spread
    cons, spread, n_coord, n_seeds = consensus({"42": 2.72, "1": 2.58, "2": 3.70, "3": 2.69})
    assert cons == 2.58
    assert round(spread, 2) == 0.14          # max-min over {2.72,2.58,2.69}, not incl. 3.70
    assert (n_coord, n_seeds) == (3, 4)


def test_all_seeds_coordinate_is_stable():
    r = assess("x", {"42": 2.71, "1": 2.71, "2": 2.70, "3": 2.72, "4": 2.69}, {})
    assert r["fungal_seeds_coordinating"] == "5/5"
    assert r["stable"] is True


def test_single_dock_artifact_is_unstable():
    # LDN-27219's real shape: tight in 2 seeds, ~5.5 A (non-coordinating) in 3 -> NOT stable
    byseed = {"42": 2.603, "4": 2.66, "1": 5.574, "2": 5.509, "3": 5.566}
    r = assess("LDN-27219", byseed, {})
    assert r["fungal_consensus"] == 2.603     # still the achievable min...
    assert r["fungal_seeds_coordinating"] == "2/5"
    assert r["stable"] is False               # ...but it did not reproduce -> demoted


def test_large_spread_breaks_stability_even_if_all_coordinate():
    # every seed coordinates (all < FE_CUTOFF) but the spread exceeds SPREAD_MAX
    lo, hi = 1.85, 2.95                       # both < 3.0, spread 1.10 > SPREAD_MAX (1.0)
    assert hi < FE_CUTOFF and (hi - lo) > SPREAD_MAX
    r = assess("y", {"42": lo, "1": hi, "2": 2.3, "3": 2.4, "4": 2.5}, {})
    assert r["fungal_seeds_coordinating"] == "5/5"
    assert r["stable"] is False


def test_never_coordinates_reports_zero_and_not_reaching():
    cons, spread, n_coord, n_seeds = consensus({"42": 4.2, "1": 3.9})
    assert cons == 3.9                        # closest approach seen, still informative
    assert n_coord == 0                       # but zero coordinating seeds
    r = assess("z", {"42": 4.2, "1": 3.9}, {})
    assert r["reaches_fungal_iron"] is False


def test_selectivity_margin_sign_and_reach():
    # reaches fungal iron closely, never reaches human iron -> positive margin, human 0/5
    r = assess("sel", {"42": 2.5, "1": 2.5, "2": 2.5, "3": 2.5, "4": 2.5},
               {"42": 4.3, "1": 4.4, "2": 4.5, "3": 4.3, "4": 4.4})
    assert r["reaches_human_iron"] is False
    assert r["selectivity_margin"] > 0        # fungal reached closer than human
    assert r["stable"] is True


def test_non_selective_azole_hits_human_too():
    # a fluconazole-like baseline: coordinates both irons about equally -> margin ~ 0, human 5/5
    r = assess("flu", {"42": 2.74, "1": 2.73, "2": 2.71, "3": 2.74, "4": 2.76},
               {"42": 2.82, "1": 2.82, "2": 2.86, "3": 2.77, "4": 2.78})
    assert r["reaches_human_iron"] is True
    assert abs(r["selectivity_margin"]) < 0.2

    assert FE_CUTOFF == 3.0
