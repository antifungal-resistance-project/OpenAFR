"""Which nitrogen actually reaches the heme iron — the mechanism behind a rank.

The validated criterion (openafr.pdbqt.min_nitrogen_iron_distance) is the closest
approach of *any* nitrogen in a pose to the heme iron. That is deliberately
element-only: it does not ask *which* nitrogen. But the method was validated on the
azole mechanism — an aromatic sp2 ring nitrogen (imidazole/triazole/pyridine)
donating its lone pair to the iron. A molecule that instead points a **nitrile**
(C#N), an **amine**, or an **amide** nitrogen at the iron can post an identical tight
N-Fe distance while coordinating through chemistry the method never validated. Its
rank is then an out-of-mechanism artifact, not azole-mimetic coordination.

This module classifies the coordinating nitrogen of a docked pose **purely from the
pose geometry** — no atom-order bookkeeping through obabel, no external labels — by
its bonded-neighbour topology, which is decisive:

    aromatic ring N  (pyridine/triazole/imidazole)  2 heavy neighbours, C-N ~1.34 A
    nitrile N        (C#N)                           1 heavy neighbour,  C#N ~1.16 A
    amine N          (-NH2 / -NR2)                   1-3 neighbours,     C-N ~1.47 A
    amide N          (-C(=O)N)                       2-3 neighbours,     C-N ~1.34 A

The one load-bearing threshold is BOND_NITRILE (1.25 A), and it sits in the wide,
empty gap between a triple bond (~1.16 A) and every single/aromatic C-N bond
(>=1.34 A) — there is nothing to tune it to. `classify_nitrogen` is deterministic and
pure; `coordinating_nitrogen` applies it to the iron-nearest nitrogen of one pose.
"""
import math

# A covalent bond between two heavy atoms is shorter than this; anything farther is a
# non-bonded contact. C-N (all orders), C-C, C-O, C-S single bonds are all < 1.85 A.
BOND_CUTOFF = 1.85
# A C#N triple bond is ~1.16 A; the next-shortest N bond (aromatic/amide C-N) is ~1.34 A.
# This cutoff sits in the empty gap between them, so it is not a tunable knob.
BOND_NITRILE = 1.25


def _neighbours(atoms, i, cutoff=BOND_CUTOFF):
    """Bonded heavy-atom neighbours of atom i: [(index, name, distance), ...]."""
    xi = atoms[i][1:]
    out = []
    for j, a in enumerate(atoms):
        if j == i:
            continue
        d = math.dist(xi, a[1:])
        if d < cutoff:
            out.append((j, a[0], d))
    return out


def classify_nitrogen(atoms, i):
    """Classify nitrogen atom i by its bonded-neighbour topology in a docked pose.

    atoms: [(atom_name, x, y, z), ...] heavy atoms of one pose (openafr.pdbqt.read_poses).
    Returns one of: 'aromatic-ring', 'nitrile', 'amine', 'amide-or-other'.

    Pure geometry — a nitrile is the only nitrogen with a single sub-1.25 A (triple)
    bond; a ring/amide nitrogen carries two heavy neighbours; a terminal amine carries
    one longer (single-bond) neighbour. 'in-mechanism' (azole-like) == 'aromatic-ring'.
    """
    nb = _neighbours(atoms, i)
    if not nb:
        return "amide-or-other"          # isolated N (no bonded heavy atom resolved)
    shortest = min(d for _, _, d in nb)
    if len(nb) == 1:
        return "nitrile" if shortest < BOND_NITRILE else "amine"
    return "aromatic-ring" if shortest < BOND_NITRILE + 0.20 else "amide-or-other"


IN_MECHANISM = "aromatic-ring"


def coordinating_nitrogen(atoms, fe):
    """The nitrogen that sets this pose's N-Fe distance, with its identity.

    Returns None if the pose has no nitrogen, else
    {'index', 'distance', 'kind', 'in_mechanism'} for the iron-nearest nitrogen —
    exactly the atom openafr.pdbqt.min_nitrogen_iron_distance measures, now labelled.
    """
    ns = [k for k, a in enumerate(atoms) if a[0].startswith("N")]
    if not ns:
        return None
    k = min(ns, key=lambda k: math.dist(atoms[k][1:], fe))
    kind = classify_nitrogen(atoms, k)
    return {
        "index": k,
        "distance": math.dist(atoms[k][1:], fe),
        "kind": kind,
        "in_mechanism": kind == IN_MECHANISM,
    }


def min_aromatic_nitrogen_iron_distance(atoms, fe):
    """Closest approach of an AROMATIC-ring nitrogen to the heme iron, or None.

    The in-mechanism twin of openafr.pdbqt.min_nitrogen_iron_distance: it restricts the
    minimum to the nitrogen kind the azole mechanism was validated on (an aromatic ring N
    donating its lone pair), so a nitrile / amine / azide nitrogen that happens to sit
    nearer the iron cannot set the distance. None means the pose has no aromatic-ring
    nitrogen reaching the pocket — the molecule cannot coordinate azole-like at all.
    """
    best = None
    for k, a in enumerate(atoms):
        if not a[0].startswith("N"):
            continue
        if classify_nitrogen(atoms, k) != IN_MECHANISM:
            continue
        d = math.dist(a[1:], fe)
        if best is None or d < best:
            best = d
    return best


# The nitrogen kinds that can donate a lone pair to heme Fe in a type-II mechanism the
# method has direct evidence for: an aromatic ring N (the azole mechanism, sp2 lone pair)
# and a nitrile N (sp lone pair, roughly colinear for sigma-donation). The nitrile is here
# not by analogy but because a *real active*, luliconazole, coordinates the iron through
# its nitrile in-pose (work/RESULTS_aromatic_criterion.md) — so excluding it demotes a
# known active. An sp3 amine (protonated / bulky at physiological pH) and a terminal azide
# are the kinds with no such evidence; they are the out-of-mechanism artifacts the
# any-N criterion silently rewards. This set is the pre-registered hypothesis of
# work/PREREGISTRATION_coordinator_type.md, NOT the frozen gate criterion.
COORDINATOR_KINDS = ("aromatic-ring", "nitrile")


def min_coordinating_nitrogen_iron_distance(atoms, fe, kinds=COORDINATOR_KINDS):
    """Closest approach to the heme iron of a *type-II-capable coordinator* nitrogen, or None.

    A middle ground between `openafr.pdbqt.min_nitrogen_iron_distance` (any N — silently
    rewards amine/azide artifacts) and `min_aromatic_nitrogen_iron_distance` (aromatic only
    — demotes the nitrile-coordinating active luliconazole). It admits exactly the nitrogen
    kinds in `kinds` (default `COORDINATOR_KINDS` = aromatic-ring + nitrile) and ignores the
    rest, so an amine/azide N that happens to sit nearest cannot set the distance while a
    nitrile — a mechanism a real active uses — still can. None means no admissible
    coordinator reaches the pocket in this pose. Pure geometry; classification is
    `classify_nitrogen`.
    """
    best = None
    for k, a in enumerate(atoms):
        if not a[0].startswith("N"):
            continue
        if classify_nitrogen(atoms, k) not in kinds:
            continue
        d = math.dist(a[1:], fe)
        if best is None or d < best:
            best = d
    return best


BAND = (2.47, 2.88)   # validated actives' iron-bound band (work/RESULTS_all_azoles.md)


def coordinator_summary(poses, fe, band=BAND):
    """Fuse the coordinator identity over all poses of one molecule into a dossier block.

    poses: iterable of pose atom lists ([(name, x, y, z), ...]). Returns
    {reported_kind, in_mechanism, reported_distance, aromatic_distance, verdict}, where the
    verdict follows PREREGISTRATION_coordinator.md:
      - in-mechanism     : the reported (nearest-N) coordinator is an aromatic ring N;
      - dual-mode        : reported coordinator is non-aromatic, but an aromatic ring N still
                           reaches inside `band` in some pose (rank inflated, azole-like mode
                           available);
      - out-of-mechanism : reported coordinator is non-aromatic and no aromatic ring N reaches
                           the band (the rank is an artifact);
      - no-coordination  : no nitrogen reached the pocket in any pose.
    """
    reported = None
    aromatic = None
    for atoms in poses:
        c = coordinating_nitrogen(atoms, fe)
        if c and (reported is None or c["distance"] < reported["distance"]):
            reported = c
        d = min_aromatic_nitrogen_iron_distance(atoms, fe)
        if d is not None and (aromatic is None or d < aromatic):
            aromatic = d

    if reported is None:
        verdict = "no-coordination"
    elif reported["in_mechanism"]:
        verdict = "in-mechanism"
    elif aromatic is not None and band[0] <= aromatic <= band[1]:
        verdict = "dual-mode"
    else:
        verdict = "out-of-mechanism"

    return {
        "reported_kind": reported["kind"] if reported else None,
        "in_mechanism": reported["in_mechanism"] if reported else None,
        "reported_distance": reported["distance"] if reported else None,
        "aromatic_distance": aromatic,
        "verdict": verdict,
    }


def nearest_by_kind(atoms, fe):
    """Best (smallest) Fe distance achieved by each nitrogen KIND in one pose.

    {kind: min_distance}. Lets a caller ask 'how close does this molecule's *aromatic*
    nitrogen get, ignoring a nitrile that happens to sit nearer' — the in-mechanism
    distance behind an out-of-mechanism rank.
    """
    best = {}
    for k, a in enumerate(atoms):
        if not a[0].startswith("N"):
            continue
        d = math.dist(a[1:], fe)
        kind = classify_nitrogen(atoms, k)
        if kind not in best or d < best[kind]:
            best[kind] = d
    return best
