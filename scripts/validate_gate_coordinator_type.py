"""Criterion comparison — does a TYPE-AWARE coordinator criterion (aromatic ring N +
nitrile N, excluding amine/azide) beat both the frozen any-N gate and the aromatic-only fix?

Pipeline stage: validate/ (a criterion-quality comparison, NOT a new gate; protocol.yaml is
read-only here and never modified — the frozen gate stands on the any-N criterion).

Background. The frozen gate ranks by `openafr.pdbqt.min_nitrogen_iron_distance` (nearest N of
ANY type). Two prior findings bracket it:
  - the any-N criterion silently rewards out-of-mechanism amine/azide/nitrile coordinators on
    a novel library (work/RESULTS_coordinator.md);
  - restricting to the AROMATIC ring N only LOWERS the run-2 gate AUC 0.794 -> 0.729, because a
    real active, luliconazole, coordinates via its NITRILE in-pose (work/RESULTS_aromatic_criterion.md).
The pre-registered hypothesis (work/PREREGISTRATION_coordinator_type.md) is the middle ground:
admit aromatic-ring AND nitrile (both type-II-capable donors; the nitrile is a mechanism a real
active uses) but NOT amine/azide. It should rescue luliconazole (beating aromatic-only) while
still pushing amine/azide decoys away (a shot at beating any-N).

Two ways to run:
  1. From docked poses (needs one cloud dock of the run-2 gate set, exhaustiveness 32):
        python scripts/validate_gate_coordinator_type.py work/screen2
     This ALSO writes a per-molecule, per-nitrogen-kind distance CACHE
     (work/candidates/run2_coordinator_distances.json). That cache is a sufficient statistic
     for ANY kind-set criterion, so it is the LAST time a criterion idea needs a fresh dock.
  2. From that cache, fully local, no dock:
        python scripts/validate_gate_coordinator_type.py --cache work/candidates/run2_coordinator_distances.json

A molecule with no admissible coordinator reaching the pocket is `unrankable` and ranked LAST,
never dropped (the run-2 rule verbatim; dropping unrankable rows is the easiest way to fake a pass).
"""
import argparse
import json
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import protocol
from openafr.coordinator import COORDINATOR_KINDS, classify_nitrogen
from openafr.pdbqt import iron_position, read_poses
from openafr.scoring import report

ACTIVES = "data/ligands/actives_holdout_final.smi"
CACHE = "work/candidates/run2_coordinator_distances.json"
KINDS = ("aromatic-ring", "nitrile", "amine", "amide-or-other")


def _per_kind_min_over_poses(pose_file, fe):
    """{kind: smallest Fe distance any N of that kind achieves across all poses}."""
    best = {}
    for _num, atoms in read_poses(pose_file):
        import math
        for k, a in enumerate(atoms):
            if not a[0].startswith("N"):
                continue
            d = math.dist(a[1:], fe)
            kind = classify_nitrogen(atoms, k)
            if kind not in best or d < best[kind]:
                best[kind] = d
    return best


def build_cache(screen, receptor, actives_file, out=CACHE):
    """Dock dir -> per-molecule {kind: min distance} cache. Written once, reused forever."""
    fe = iron_position(receptor)
    actives = {l.split("\t")[1].strip() for l in open(actives_file) if "\t" in l}
    rows = {}
    for f in sorted(os.listdir(screen)):
        if not f.endswith(".pdbqt"):
            continue
        name = f[:-6]
        rows[name] = {"is_active": name in actives,
                      "kind_min": _per_kind_min_over_poses(os.path.join(screen, f), fe)}
    pathlib.Path(out).parent.mkdir(parents=True, exist_ok=True)
    json.dump({"receptor": receptor, "actives_file": actives_file, "molecules": rows},
              open(out, "w"), indent=2)
    print(f"wrote per-kind distance cache: {out}  ({len(rows)} molecules)")
    return rows


def _min_over_kinds(kind_min, kinds):
    """Smallest distance over the admitted kinds, or None if none is present."""
    ds = [kind_min[k] for k in kinds if k in kind_min]
    return min(ds) if ds else None


def grade(molecules):
    """molecules: {name: {'is_active': bool, 'kind_min': {kind: dist}}}. Grade 3 criteria."""
    gate = protocol.load()["gate"]
    min_auc, min_ef1 = gate["min_auc"], gate["min_ef1"]
    ef_fractions, bedroc_alpha = gate["ef_fractions"], gate["bedroc_alpha"]

    def rows(kinds):
        return [(n, _min_over_kinds(m["kind_min"], kinds), m["is_active"])
                for n, m in molecules.items()]

    any_rows = rows(KINDS)                    # any nitrogen — the frozen gate
    arom_rows = rows(("aromatic-ring",))      # aromatic only — the fix that fell
    type_rows = rows(COORDINATOR_KINDS)       # aromatic + nitrile — the hypothesis

    n_act = sum(1 for m in molecules.values() if m["is_active"])
    print("=" * 74)
    print("CRITERION COMPARISON — any-N vs aromatic-N vs type-aware (aromatic+nitrile)")
    print("=" * 74)
    protocol.verify()
    print(f"actives : {n_act}   pool : {len(molecules)}")
    print(f"pass bar: AUC >= {min_auc} AND EF@1% >= {min_ef1}x\n")

    a = report("any-N (FROZEN gate criterion — reference)", any_rows, False,
               ef_fractions, bedroc_alpha)
    b = report("aromatic-N only (the fix that fell to 0.729)", arom_rows, False,
               ef_fractions, bedroc_alpha)
    c = report("type-aware: aromatic + nitrile (THE HYPOTHESIS)", type_rows, True,
               ef_fractions, bedroc_alpha)

    print("\n" + "=" * 74)
    print(f"any-N        : AUC {a[0]:.3f}   EF@1% {a[1][0]:.2f}x")
    print(f"aromatic-N   : AUC {b[0]:.3f}   EF@1% {b[1][0]:.2f}x")
    print(f"type-aware   : AUC {c[0]:.3f}   EF@1% {c[1][0]:.2f}x")
    beats_any = c[0] >= a[0] - 1e-9
    beats_arom = c[0] >= b[0] - 1e-9
    passes = c[0] >= min_auc and c[1][0] >= min_ef1
    print(f"type-aware >= any-N     : {'YES' if beats_any else 'NO'}  ({c[0]-a[0]:+.3f})")
    print(f"type-aware >= aromatic  : {'YES' if beats_arom else 'NO'}  ({c[0]-b[0]:+.3f})")
    print(f"type-aware still PASSES : {'YES' if passes else 'NO'}")
    # Pre-registered adoption rule: adopt for NOVEL-CANDIDATE ranking only (never the frozen
    # gate) iff it does not fall below any-N AND still passes. Beating aromatic-only alone is
    # necessary but not sufficient.
    print(f"VERDICT (adopt for novel-candidate ranking): "
          f"{'YES' if (beats_any and passes) else 'NO'}")
    return 0 if (beats_any and passes) else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("screen", nargs="?", help="docked pose dir (run-2 gate set); builds+uses the cache")
    ap.add_argument("--cache", help="grade from an existing per-kind distance cache (no dock)")
    ap.add_argument("--receptor", default="work/receptor_A.pdb")
    ap.add_argument("--actives", default=ACTIVES)
    args = ap.parse_args()
    if args.cache:
        molecules = json.load(open(args.cache))["molecules"]
    elif args.screen:
        molecules = build_cache(args.screen, args.receptor, args.actives)
    else:
        ap.error("give a docked pose dir (to build the cache) or --cache <json> (to re-grade)")
    return grade(molecules)


if __name__ == "__main__":
    sys.exit(main())
