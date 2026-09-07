"""Known-answer tests for the coordinating-nitrogen classifier (openafr/coordinator.py).

These lock the pure topological rule frozen in work/PREREGISTRATION_coordinator.md:
  - a nitrogen with ONE heavy neighbour bonded shorter than 1.25 A is a nitrile (C#N);
  - a nitrogen with ONE longer (single-bond) neighbour is an amine;
  - a nitrogen with TWO+ neighbours at aromatic/amide bond length is a ring nitrogen;
  - 'in-mechanism' (azole-like) == aromatic-ring only;
  - coordinating_nitrogen() labels exactly the iron-nearest nitrogen — the atom the
    validated criterion (openafr.pdbqt.min_nitrogen_iron_distance) measures.
"""
import math

from openafr.coordinator import (
    BAND,
    classify_nitrogen,
    coordinating_nitrogen,
    coordinator_summary,
    min_aromatic_nitrogen_iron_distance,
    nearest_by_kind,
)


def _nitrile():
    """C-C#N laid on the x-axis: N at x=0, nitrile C at 1.16, chain C at 2.63."""
    return [("N", 0.0, 0.0, 0.0), ("C", 1.16, 0.0, 0.0), ("C", 2.63, 0.0, 0.0)]


def _pyridine_n():
    """A ring nitrogen with two aromatic C neighbours at ~1.34 A."""
    return [("N", 0.0, 0.0, 0.0), ("C", 1.34, 0.0, 0.0), ("C", -0.67, 1.16, 0.0)]


def _amine():
    """A terminal amine: one C neighbour at a single-bond 1.47 A."""
    return [("N", 0.0, 0.0, 0.0), ("C", 1.47, 0.0, 0.0)]


def test_nitrile_is_nitrile():
    assert classify_nitrogen(_nitrile(), 0) == "nitrile"


def test_pyridine_is_aromatic_ring():
    assert classify_nitrogen(_pyridine_n(), 0) == "aromatic-ring"


def test_amine_is_amine_not_nitrile():
    # one neighbour, but the bond is a single bond (1.47 A) -> not a triple bond
    assert classify_nitrogen(_amine(), 0) == "amine"


def test_isolated_nitrogen_is_other():
    # a nitrogen with no heavy neighbour inside the bond cutoff
    assert classify_nitrogen([("N", 0.0, 0.0, 0.0), ("C", 5.0, 0.0, 0.0)], 0) == "amide-or-other"


def test_coordinating_picks_iron_nearest_nitrogen():
    # two nitrogens on separated rays from Fe: an aromatic ring N at 2.5 A (+x),
    # a nitrile N at 2.0 A (+y). The criterion measures the NEAREST -> the nitrile,
    # and flags it out-of-mechanism.
    fe = (0.0, 0.0, 0.0)
    atoms = [
        ("N", 2.5, 0.0, 0.0), ("C", 3.84, 0.0, 0.0), ("C", 2.5, 1.34, 0.0),  # ring N (idx 0)
        ("N", 0.0, 2.0, 0.0), ("C", 0.0, 3.16, 0.0),                         # nitrile N (idx 3)
    ]
    c = coordinating_nitrogen(atoms, fe)
    assert c["index"] == 3
    assert c["kind"] == "nitrile"
    assert c["in_mechanism"] is False
    assert c["distance"] == math.dist((0.0, 2.0, 0.0), fe)


def test_coordinating_none_when_no_nitrogen():
    assert coordinating_nitrogen([("C", 0.0, 0.0, 0.0)], (1.0, 0.0, 0.0)) is None


def test_nearest_by_kind_separates_the_two_nitrogens():
    # the in-mechanism (aromatic) distance behind an out-of-mechanism (nitrile) rank
    fe = (0.0, 0.0, 0.0)
    atoms = [
        ("N", 2.79, 0.0, 0.0), ("C", 4.13, 0.0, 0.0), ("C", 2.79, 1.34, 0.0),  # ring N (+x)
        ("N", 0.0, 2.42, 0.0), ("C", 0.0, 3.58, 0.0),                          # nitrile N (+y)
    ]
    best = nearest_by_kind(atoms, fe)
    assert best["nitrile"] == math.dist((0.0, 2.42, 0.0), fe)
    assert best["aromatic-ring"] == math.dist((2.79, 0.0, 0.0), fe)
    # the reported rank rides on the nitrile, but the aromatic N is farther
    assert best["nitrile"] < best["aromatic-ring"]


def test_aromatic_min_ignores_a_nearer_nitrile():
    # a nitrile N sits nearest the iron, an aromatic ring N farther. The any-N criterion
    # would report the nitrile; the aromatic-restricted one reports the ring N.
    fe = (0.0, 0.0, 0.0)
    atoms = [
        ("N", 2.79, 0.0, 0.0), ("C", 4.13, 0.0, 0.0), ("C", 2.79, 1.34, 0.0),  # ring N
        ("N", 0.0, 2.42, 0.0), ("C", 0.0, 3.58, 0.0),                          # nitrile N (nearer)
    ]
    assert min_aromatic_nitrogen_iron_distance(atoms, fe) == math.dist((2.79, 0.0, 0.0), fe)


def test_aromatic_min_none_when_only_nonaromatic_nitrogens():
    # a molecule whose only nitrogen is a nitrile cannot coordinate azole-like -> None
    fe = (0.0, 0.0, 0.0)
    atoms = [("N", 2.0, 0.0, 0.0), ("C", 3.16, 0.0, 0.0)]
    assert min_aromatic_nitrogen_iron_distance(atoms, fe) is None


# --- coordinator_summary: the four pre-registered verdicts (PREREGISTRATION_coordinator.md) ---
#
# Fe sits at the origin. Pose builders below place a nitrogen of a given KIND at a chosen
# distance from the iron: a ring N carries two ~1.34 A C neighbours (aromatic-ring), a
# nitrile N carries one 1.16 A C neighbour (out-of-mechanism). BAND = (2.47, 2.88) is the
# validated actives' iron-bound distance window.

_FE = (0.0, 0.0, 0.0)


def _ring_at(d):
    """An aromatic ring N at distance d from Fe on +x (two C neighbours at ~1.34 A)."""
    return [("N", d, 0.0, 0.0), ("C", d + 1.34, 0.0, 0.0), ("C", d, 1.34, 0.0)]


def _nitrile_at(d):
    """A nitrile N at distance d from Fe on +y (one triple-bonded C at 1.16 A)."""
    return [("N", 0.0, d, 0.0), ("C", 0.0, d + 1.16, 0.0)]


def test_summary_in_mechanism_when_nearest_n_is_aromatic():
    # the reported (iron-nearest) coordinator is an aromatic ring N inside the band
    s = coordinator_summary([_ring_at(2.6)], _FE)
    assert s["verdict"] == "in-mechanism"
    assert s["reported_kind"] == "aromatic-ring"
    assert s["in_mechanism"] is True
    assert s["reported_distance"] == 2.6
    assert s["aromatic_distance"] == 2.6


def test_summary_out_of_mechanism_when_nitrile_reports_and_aromatic_misses_band():
    # nitrile nearer the iron sets the rank; the only aromatic N sits past the band -> artifact
    pose = _nitrile_at(2.0) + _ring_at(3.5)   # aromatic 3.5 A is above BAND (2.88)
    s = coordinator_summary([pose], _FE)
    assert s["verdict"] == "out-of-mechanism"
    assert s["reported_kind"] == "nitrile"
    assert s["in_mechanism"] is False
    assert s["reported_distance"] == 2.0
    assert s["aromatic_distance"] == 3.5


def test_summary_dual_mode_when_nitrile_reports_but_aromatic_reaches_band():
    # nitrile sets the (inflated) rank, but an aromatic ring N still reaches inside the band
    pose = _nitrile_at(2.0) + _ring_at(2.7)   # aromatic 2.7 A is within BAND
    s = coordinator_summary([pose], _FE)
    assert s["verdict"] == "dual-mode"
    assert s["reported_kind"] == "nitrile"
    assert BAND[0] <= s["aromatic_distance"] <= BAND[1]


def test_summary_no_coordination_when_no_nitrogen_reaches():
    # no nitrogen in any pose -> nothing coordinates the iron
    poses = [[("C", 1.0, 0.0, 0.0)], [("C", 0.0, 2.0, 0.0), ("O", 0.0, 3.0, 0.0)]]
    s = coordinator_summary(poses, _FE)
    assert s["verdict"] == "no-coordination"
    assert s["reported_kind"] is None
    assert s["in_mechanism"] is None
    assert s["reported_distance"] is None
    assert s["aromatic_distance"] is None


def test_summary_fuses_across_poses_reported_and_aromatic_from_different_poses():
    # the reported (nearest) coordinator and the band-reaching aromatic N come from
    # DIFFERENT poses: the fusion must take the min of each across the whole ensemble.
    pose_a = _nitrile_at(2.0)          # nearest coordinator overall (a nitrile)
    pose_b = _ring_at(2.7)             # aromatic N reaching the band, in another pose
    s = coordinator_summary([pose_a, pose_b], _FE)
    assert s["reported_kind"] == "nitrile"
    assert s["reported_distance"] == 2.0
    assert s["aromatic_distance"] == 2.7
    assert s["verdict"] == "dual-mode"


# --- the type-aware coordinator criterion (min_coordinating_nitrogen_iron_distance) ---
# Hypothesis of work/PREREGISTRATION_coordinator_type.md: admit aromatic-ring AND nitrile
# (both type-II-capable donors, and a real active — luliconazole — uses its nitrile), but
# NOT amine/azide. These lock that the criterion admits exactly that kind set.

def _amine_at(d):
    """An sp3 amine N at distance d from Fe on +z (one single-bonded C at 1.47 A)."""
    return [("N", 0.0, 0.0, d), ("C", 0.0, 0.0, d + 1.47)]


def test_type_aware_admits_aromatic():
    from openafr.coordinator import min_coordinating_nitrogen_iron_distance
    assert min_coordinating_nitrogen_iron_distance(_ring_at(2.6), _FE) == 2.6


def test_type_aware_admits_nitrile():
    # the key difference from the aromatic-only criterion: a nitrile still sets the distance
    from openafr.coordinator import min_coordinating_nitrogen_iron_distance
    assert min_coordinating_nitrogen_iron_distance(_nitrile_at(2.4), _FE) == 2.4
    assert min_aromatic_nitrogen_iron_distance(_nitrile_at(2.4), _FE) is None


def test_type_aware_excludes_amine():
    # an amine nearer the iron must NOT set the distance (out-of-mechanism kind)
    from openafr.coordinator import min_coordinating_nitrogen_iron_distance
    assert min_coordinating_nitrogen_iron_distance(_amine_at(2.0), _FE) is None


def test_type_aware_picks_nitrile_over_nearer_amine():
    # amine sits nearest (2.0) but is inadmissible; the nitrile (2.5) sets the distance
    from openafr.coordinator import min_coordinating_nitrogen_iron_distance
    pose = _amine_at(2.0) + _nitrile_at(2.5)
    assert min_coordinating_nitrogen_iron_distance(pose, _FE) == 2.5


def test_type_aware_sits_between_any_n_and_aromatic_only():
    from openafr.coordinator import min_coordinating_nitrogen_iron_distance
    from openafr.pdbqt import min_nitrogen_iron_distance
    # amine nearest (1.8), then nitrile (2.4), then aromatic (2.9)
    pose = _amine_at(1.8) + _nitrile_at(2.4) + _ring_at(2.9)
    assert min_nitrogen_iron_distance(pose, _FE) == 1.8            # any-N: the amine
    assert min_coordinating_nitrogen_iron_distance(pose, _FE) == 2.4  # type-aware: the nitrile
    assert min_aromatic_nitrogen_iron_distance(pose, _FE) == 2.9   # aromatic-only: the ring N


def test_type_aware_kinds_argument_is_honoured():
    # restricting kinds to aromatic-only reproduces the aromatic criterion
    from openafr.coordinator import min_coordinating_nitrogen_iron_distance
    pose = _nitrile_at(2.4) + _ring_at(2.9)
    assert min_coordinating_nitrogen_iron_distance(pose, _FE, kinds=("aromatic-ring",)) == 2.9
