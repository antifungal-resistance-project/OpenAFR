"""Per-candidate stereochemistry axis (issue #85).

Grades whether each shortlist rank is stereochemistry-explicit and, where it is not,
whether the iron-coordinating geometry survives across all enumerated isomers. This is a
per-candidate CONFIDENCE AXIS, not a gate — it annotates, never promotes. Criterion + bar
are frozen in work/PREREGISTRATION_candidate_stereo.md
(sha work/PREREG_candidate_stereo.sha256) BEFORE any per-isomer distance was read.

Inputs (all frozen upstream):
  - work/candidates/stereo_audit.json          the deterministic enumeration
                                               (scripts/enumerate_stereo.py)
  - work/candidates/stereo_ligands/stereo_map.json   isomer tag -> candidate + SMILES
  - work/candidates/dock_stereo/s<seed>/*.pdbqt      per-isomer poses (screen_consensus.sh
                                               into the WT receptor work/receptor.pdbqt)

Per-isomer consensus and per-candidate spread reuse the frozen consensus() logic from
scripts/shortlist_confidence.py unchanged. Verdicts (frozen):
  EXPLICIT           supplied SMILES pins a single isomer; rank is stereo-explicit (no dock).
  ROBUST             every isomer coordinates (>= MIN_COORD_SEEDS) AND iso_spread <= ISO_TOL.
  ISOMER-DEPENDENT   isomers disagree (some coordinate, some don't) OR iso_spread > ISO_TOL.

Usage:
    python scripts/candidate_stereo.py \
        [--audit work/candidates/stereo_audit.json] \
        [--stereo-map work/candidates/stereo_ligands/stereo_map.json] \
        [--dock work/candidates/dock_stereo] \
        [--receptor work/receptor.pdbqt] [--json]
"""
import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr.pdbqt import iron_position
# Reuse the frozen consensus logic unchanged — do not redefine it here.
from scripts.shortlist_confidence import FE_CUTOFF, consensus, per_seed_distances

ISO_TOL = 1.0        # A; per-candidate isomer spread beyond this = ISOMER-DEPENDENT (frozen)
MIN_COORD_SEEDS = 3  # coordinating seeds for an isomer to count as reaching the iron (frozen)


def classify(iso_conses):
    """(verdict, iso_spread) from a list of (iso_tag, cons, ncoord) for one candidate."""
    coord = [c for _t, c, n in iso_conses if n >= MIN_COORD_SEEDS and c is not None]
    all_coord = all(n >= MIN_COORD_SEEDS for _t, _c, n in iso_conses)
    spread = (max(coord) - min(coord)) if len(coord) > 1 else 0.0
    if not all_coord:
        return "ISOMER-DEPENDENT", spread
    if spread > ISO_TOL:
        return "ISOMER-DEPENDENT", spread
    return "ROBUST", spread


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--audit", default="work/candidates/stereo_audit.json")
    ap.add_argument("--stereo-map", default="work/candidates/stereo_ligands/stereo_map.json")
    ap.add_argument("--dock", default="work/candidates/dock_stereo")
    # Iron position is read from the .pdb anchor (same coordinate frame as the docking
    # receptor work/receptor.pdbqt), matching how every other candidate axis locates the iron.
    ap.add_argument("--receptor", default="work/receptor_A.pdb")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    audit = json.load(open(args.audit))
    smap = {}
    if pathlib.Path(args.stereo_map).exists():
        smap = json.load(open(args.stereo_map))

    # Per-isomer consensus from the dock tree (only exists for ambiguous candidates).
    byseed = {}
    if pathlib.Path(args.dock).exists():
        fe = iron_position(args.receptor)
        byseed = per_seed_distances(args.dock, fe)

    # Group docked isomers by their candidate.
    iso_by_cand = {}
    for tag, cons_seed in byseed.items():
        cand = smap.get(tag, {}).get("candidate")
        if cand is None:
            continue
        cons, spread, ncoord, nseed = consensus(cons_seed)
        iso_by_cand.setdefault(cand, []).append({
            "isomer": tag,
            "smiles": smap[tag]["smiles"],
            "consensus": None if cons is None else round(cons, 3),
            "seeds_coordinating": f"{ncoord}/{nseed}",
            "spread": round(spread, 3),
            "per_seed": {k: (None if v is None else round(v, 3))
                         for k, v in sorted(cons_seed.items())},
        })

    recs = []
    for c in audit["candidates"]:
        name = c["name"]
        if c["stereo_explicit"]:
            recs.append({
                "name": name,
                "stereo_explicit": True,
                "n_isomers": c["n_isomers"],
                "n_unassigned": c["n_unassigned"],
                "n_unassigned_bonds": c["n_unassigned_bonds"],
                "verdict": "EXPLICIT",
                "iso_spread": None,
                "isomers": [],
            })
            continue
        isos = sorted(iso_by_cand.get(name, []), key=lambda r: r["isomer"])
        triples = [(r["isomer"], r["consensus"],
                    int(r["seeds_coordinating"].split("/")[0])) for r in isos]
        verdict, spread = classify(triples) if triples else ("UNDOCKED", None)
        recs.append({
            "name": name,
            "stereo_explicit": False,
            "n_isomers": c["n_isomers"],
            "n_unassigned": c["n_unassigned"],
            "n_unassigned_bonds": c["n_unassigned_bonds"],
            "verdict": verdict,
            "iso_spread": None if spread is None else round(spread, 3),
            "isomers": isos,
        })

    out = {"fe_cutoff": FE_CUTOFF, "iso_tol": ISO_TOL,
           "min_coord_seeds": MIN_COORD_SEEDS,
           "n_stereo_explicit": sum(1 for r in recs if r["stereo_explicit"]),
           "candidates": recs}

    if args.json:
        print(json.dumps(out, indent=2))
        return 0

    print("=" * 90)
    print("CANDIDATE STEREOCHEMISTRY — supplied-SMILES isomer audit + isomer dock (#85)")
    print(f"receptor: {args.receptor}   cutoff {FE_CUTOFF} A   "
          f"iso_tol {ISO_TOL} A   min {MIN_COORD_SEEDS}/5 seeds")
    print("=" * 90)
    print("%-16s %8s %10s %10s   %s" %
          ("candidate", "isomers", "unassigned", "iso_spread", "verdict"))
    print("-" * 90)
    for r in recs:
        sp = "-" if r["iso_spread"] is None else "%.3f" % r["iso_spread"]
        print("%-16s %8d %10d %10s   %s" %
              (r["name"], r["n_isomers"], r["n_unassigned"], sp, r["verdict"]))
        for iso in r["isomers"]:
            print("    %-24s %s A  %s" %
                  (iso["isomer"], iso["consensus"], iso["seeds_coordinating"]))
    print("-" * 90)
    print(f"stereo-explicit: {out['n_stereo_explicit']}/{len(recs)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
