"""Coordinator-identity axis — is a candidate's iron-reach azole-like, or an artifact?

Pipeline stage: candidate confidence (per-candidate axis, annotation not a gate).

The validated criterion measures the closest approach of ANY nitrogen to the heme iron
(openafr.pdbqt.min_nitrogen_iron_distance). This axis re-reads the same frozen consensus
docks and asks *which* nitrogen wins that minimum — the validated azole mechanism (an
aromatic ring N) or an out-of-mechanism nitrile/amine/amide N that posts an identical
distance through chemistry the method never validated.

For each candidate, over all 20 poses x 5 seeds, it reports:
  - reported coordinator kind (the min-N-Fe nitrogen of the frozen seed-42 dock) + modal
    kind across seeds;
  - the best (smallest) Fe distance an AROMATIC-ring nitrogen reaches (the in-mechanism
    distance behind the rank);
  - a verdict: in-mechanism / dual-mode / out-of-mechanism (per PREREGISTRATION_coordinator.md).

Known azole controls (scorecard.CONTROLS) MUST come back in-mechanism — they coordinate
through a triazole/imidazole ring N by construction; a control labelled out-of-mechanism
means the classifier is wrong, not the molecule. Controls are exempt from demotion.

Pure classification lives in openafr/coordinator.py (locked by tests/test_coordinator.py);
this script only reads pose files and tabulates.
"""
import argparse
import collections
import json
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr.coordinator import IN_MECHANISM, coordinating_nitrogen, nearest_by_kind
from openafr.pdbqt import iron_position, read_poses
from scripts.scorecard import CONTROLS

# The validated actives' iron-bound band (work/RESULTS_all_azoles.md): an aromatic N that
# reaches inside it coordinates like the real drugs. Used to separate dual-mode (aromatic N
# still lands in-band) from out-of-mechanism (it does not).
BAND = (2.47, 2.88)
FE_CUTOFF = 3.0
DEFAULT_SEEDS = ["42", "1", "2", "3", "4"]


def _pose_path(docks, seed, name):
    return os.path.join(docks, f"s{seed}", f"{name}.pdbqt")


def analyze_candidate(docks, name, fe, seeds=DEFAULT_SEEDS):
    """Coordinator-identity readout for one candidate across its per-seed docks.

    Returns a dict, or None if the frozen seed's pose file is missing.
    """
    frozen = _pose_path(docks, seeds[0], name)
    if not os.path.isfile(frozen):
        return None

    # reported coordinator: the min-N-Fe nitrogen of the frozen (first) seed's dock.
    reported = _reported_coordinator(frozen, fe)

    # modal kind across seeds (each seed's own min-N-Fe pose).
    per_seed_kind = {}
    best_by_kind = {}          # smallest Fe distance per nitrogen KIND across all poses/seeds
    for seed in seeds:
        p = _pose_path(docks, seed, name)
        if not os.path.isfile(p):
            continue
        seed_min = None
        for _, atoms in read_poses(p):
            c = coordinating_nitrogen(atoms, fe)
            if c and (seed_min is None or c["distance"] < seed_min["distance"]):
                seed_min = c
            for kind, d in nearest_by_kind(atoms, fe).items():
                if kind not in best_by_kind or d < best_by_kind[kind]:
                    best_by_kind[kind] = d
        if seed_min:
            per_seed_kind[seed] = seed_min["kind"]

    modal = collections.Counter(per_seed_kind.values()).most_common(1)
    modal_kind = modal[0][0] if modal else None
    in_mech_dist = best_by_kind.get(IN_MECHANISM)
    reported_dist = reported["distance"] if reported else None
    gap = (in_mech_dist - reported_dist) if (in_mech_dist is not None
                                             and reported_dist is not None) else None

    return {
        "name": name,
        "is_control": name in CONTROLS,
        "reported_kind": reported["kind"] if reported else None,
        "reported_distance": reported_dist,
        "reported_in_mechanism": reported["in_mechanism"] if reported else None,
        "modal_kind": modal_kind,
        "per_seed_kind": per_seed_kind,
        "best_by_kind": {k: round(v, 3) for k, v in sorted(best_by_kind.items())},
        "in_mechanism_distance": in_mech_dist,
        "gap": gap,
        "verdict": _verdict(reported, in_mech_dist),
    }


def _reported_coordinator(pose_file, fe):
    """The min-N-Fe coordinating nitrogen over all poses of one dock file."""
    best = None
    for _, atoms in read_poses(pose_file):
        c = coordinating_nitrogen(atoms, fe)
        if c and (best is None or c["distance"] < best["distance"]):
            best = c
    return best


def _verdict(reported, in_mech_dist):
    """Frozen verdict rule (PREREGISTRATION_coordinator.md). Pure."""
    if reported is None:
        return "no-coordination"
    if reported["in_mechanism"]:
        return "in-mechanism"
    # reported coordinator is non-aromatic -> is an azole-like mode nonetheless available?
    if in_mech_dist is not None and BAND[0] <= in_mech_dist <= BAND[1]:
        return "dual-mode"
    return "out-of-mechanism"


def read_names(smi_path):
    names = []
    for line in open(smi_path):
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 2:
            names.append(parts[1])
    return names


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--docks", default="work/candidates/dock_coordinator",
                    help="parent dir with per-seed pose subdirs s<seed>/")
    ap.add_argument("--ligands", default="work/candidates/shortlist_focus.smi",
                    help="SMILES<TAB>name library (for candidate names/order)")
    ap.add_argument("--receptor-pdb", default="work/receptor_A.pdb",
                    help="receptor carrying the heme iron")
    ap.add_argument("--seeds", nargs="+", default=DEFAULT_SEEDS)
    ap.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    args = ap.parse_args()

    fe = iron_position(args.receptor_pdb)
    rows = []
    for name in read_names(args.ligands):
        r = analyze_candidate(args.docks, name, fe, args.seeds)
        if r is not None:
            rows.append(r)

    rows.sort(key=lambda r: (r["reported_distance"] is None, r["reported_distance"]))

    if args.json:
        print(json.dumps({"axis": "coordinator_identity", "band": list(BAND),
                          "seeds": args.seeds, "candidates": rows}, indent=2))
        return

    print(f"Heme iron at ({fe[0]:.2f}, {fe[1]:.2f}, {fe[2]:.2f})  "
          f"validated aromatic-N band {BAND[0]:.2f}-{BAND[1]:.2f} A\n")
    hdr = f"{'candidate':16}{'report N-Fe':>12}{'coord kind':>16}{'aromatic N-Fe':>15}{'verdict':>17}"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        rd = "-" if r["reported_distance"] is None else f"{r['reported_distance']:.3f}"
        am = "-" if r["in_mechanism_distance"] is None else f"{r['in_mechanism_distance']:.3f}"
        tag = r["verdict"] + ("*" if r["is_control"] else "")
        print(f"{r['name']:16}{rd:>12}{r['reported_kind'] or '-':>16}{am:>15}{tag:>17}")
    print("\n* = known azole control (validates the classifier; exempt from demotion)")


if __name__ == "__main__":
    main()
