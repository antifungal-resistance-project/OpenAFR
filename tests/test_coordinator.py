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
    classify_nitrogen,
    coordinating_nitrogen,
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
