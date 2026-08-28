"""Per-candidate protonation/tautomer axis (issue #84).

Grades whether each shortlist rank is protonation/tautomer-explicit and, where it is not,
whether the iron-coordinating geometry survives across all enumerated microstates. This is a
per-candidate CONFIDENCE AXIS, not a gate — it annotates, never promotes. Criterion + bar are
frozen in work/PREREGISTRATION_candidate_protomer.md (sha work/PREREG_candidate_protomer.sha256)
BEFORE any per-microstate distance was read. Twin of the stereochemistry axis (#85),
scripts/candidate_stereo.py — same machinery, protonation/tautomer instead of stereoisomers.

Inputs (all frozen upstream):
  - work/candidates/protomer_audit.json               the deterministic enumeration
                                                      (scripts/enumerate_protomers.py)
  - work/candidates/protomer_ligands/protomer_map.json  microstate tag -> candidate + SMILES
  - work/candidates/dock_protomer/s<seed>/*.pdbqt     per-microstate poses (screen_consensus.sh
                                                      into the WT receptor work/receptor.pdbqt)

Per-microstate consensus and per-candidate spread reuse the frozen consensus() logic from
scripts/shortlist_confidence.py unchanged. Verdicts (frozen):
  EXPLICIT               single relevant microstate at pH 7.4; rank is state-explicit (no dock).
  ROBUST                 every microstate coordinates (>= MIN_COORD_SEEDS) AND spread <= STATE_TOL.
  MICROSTATE-DEPENDENT   microstates disagree (some coordinate, some don't) OR spread > STATE_TOL.

Usage:
    python scripts/candidate_protomer.py \
        [--audit work/candidates/protomer_audit.json] \
        [--protomer-map work/candidates/protomer_ligands/protomer_map.json] \
        [--dock work/candidates/dock_protomer] \
        [--receptor work/receptor_A.pdb] [--json]
"""
import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr.pdbqt import iron_position
# Reuse the frozen consensus logic unchanged — do not redefine it here.
from scripts.shortlist_confidence import FE_CUTOFF, consensus, per_seed_distances

STATE_TOL = 1.0      # A; per-candidate microstate spread beyond this = MICROSTATE-DEPENDENT (frozen)
MIN_COORD_SEEDS = 3  # coordinating seeds for a microstate to count as reaching the iron (frozen)


def classify(state_conses):
    """(verdict, state_spread) from a list of (tag, cons, ncoord) for one candidate."""
    coord = [c for _t, c, n in state_conses if n >= MIN_COORD_SEEDS and c is not None]
    all_coord = all(n >= MIN_COORD_SEEDS for _t, _c, n in state_conses)
    spread = (max(coord) - min(coord)) if len(coord) > 1 else 0.0
    if not all_coord:
        return "MICROSTATE-DEPENDENT", spread
    if spread > STATE_TOL:
        return "MICROSTATE-DEPENDENT", spread
    return "ROBUST", spread


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--audit", default="work/candidates/protomer_audit.json")
    ap.add_argument("--protomer-map",
                    default="work/candidates/protomer_ligands/protomer_map.json")
    ap.add_argument("--dock", default="work/candidates/dock_protomer")
    # Iron position is read from the .pdb anchor (same coordinate frame as the docking
    # receptor work/receptor.pdbqt), matching how every other candidate axis locates the iron.
    ap.add_argument("--receptor", default="work/receptor_A.pdb")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    audit = json.load(open(args.audit))
    pmap = {}
    if pathlib.Path(args.protomer_map).exists():
        pmap = json.load(open(args.protomer_map))

    # Per-microstate consensus from the dock tree (only exists for multi-state candidates).
    byseed = {}
    if pathlib.Path(args.dock).exists():
        fe = iron_position(args.receptor)
        byseed = per_seed_distances(args.dock, fe)

    # Group docked microstates by their candidate.
    state_by_cand = {}
    for tag, cons_seed in byseed.items():
        cand = pmap.get(tag, {}).get("candidate")
        if cand is None:
            continue
        cons, spread, ncoord, nseed = consensus(cons_seed)
        state_by_cand.setdefault(cand, []).append({
            "state": tag,
            "smiles": pmap[tag]["smiles"],
            "is_reference": pmap[tag].get("is_reference", False),
            "consensus": None if cons is None else round(cons, 3),
            "seeds_coordinating": f"{ncoord}/{nseed}",
            "spread": round(spread, 3),
            "per_seed": {k: (None if v is None else round(v, 3))
                         for k, v in sorted(cons_seed.items())},
        })

    recs = []
    for c in audit["candidates"]:
        name = c["name"]
        if c["protomer_explicit"]:
            recs.append({
                "name": name,
                "protomer_explicit": True,
                "n_states": c["n_states"],
                "n_protonation_states": c["n_protonation_states"],
                "verdict": "EXPLICIT",
                "state_spread": None,
                "states": [],
            })
            continue
        states = sorted(state_by_cand.get(name, []), key=lambda r: r["state"])
        triples = [(r["state"], r["consensus"],
                    int(r["seeds_coordinating"].split("/")[0])) for r in states]
        verdict, spread = classify(triples) if triples else ("UNDOCKED", None)
        recs.append({
            "name": name,
            "protomer_explicit": False,
            "n_states": c["n_states"],
            "n_protonation_states": c["n_protonation_states"],
            "verdict": verdict,
            "state_spread": None if spread is None else round(spread, 3),
            "states": states,
        })

    out = {"fe_cutoff": FE_CUTOFF, "state_tol": STATE_TOL,
           "min_coord_seeds": MIN_COORD_SEEDS,
           "n_protomer_explicit": sum(1 for r in recs if r["protomer_explicit"]),
           "candidates": recs}

    if args.json:
        print(json.dumps(out, indent=2))
        return 0

    print("=" * 92)
    print("CANDIDATE PROTONATION / TAUTOMER — pH-7.4 microstate audit + microstate dock (#84)")
    print(f"receptor: {args.receptor}   cutoff {FE_CUTOFF} A   "
          f"state_tol {STATE_TOL} A   min {MIN_COORD_SEEDS}/5 seeds")
    print("=" * 92)
    print("%-16s %8s %10s %10s   %s" %
          ("candidate", "states", "protStates", "spread", "verdict"))
    print("-" * 92)
    for r in recs:
        sp = "-" if r["state_spread"] is None else "%.3f" % r["state_spread"]
        print("%-16s %8d %10d %10s   %s" %
              (r["name"], r["n_states"], r["n_protonation_states"], sp, r["verdict"]))
        for st in r["states"]:
            ref = " (ref)" if st["is_reference"] else ""
            print("    %-24s %s A  %s%s" %
                  (st["state"], st["consensus"], st["seeds_coordinating"], ref))
    print("-" * 92)
    print(f"protomer-explicit: {out['n_protomer_explicit']}/{len(recs)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
