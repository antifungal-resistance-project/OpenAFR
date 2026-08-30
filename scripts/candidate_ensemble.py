"""Per-candidate ensemble-consistency (issue #70).

Grades whether the shortlist candidates keep their iron-coordinating geometry across the
EXHAUSTIVE crystallographic *C. albicans* CYP51 ensemble {5TZ1, 5FSA, 5V5Z}, rather than
only in the single 5TZ1 conformer they were ranked against. This is a per-candidate
CONFIDENCE AXIS, not a gate — it demotes/annotates, never promotes. Criterion + bar are
frozen in work/PREREGISTRATION_candidate_ensemble.md
(sha work/PREREG_candidate_ensemble.sha256) BEFORE any 5FSA/5V5Z distance was read.

Method (all frozen in the pre-registration):
  - Per conformer, per molecule: the SAME consensus() logic as scripts/shortlist_confidence.py
    (MIN nitrogen-to-iron distance over the seeds that coordinated, < FE_CUTOFF 3.0 A, with
    n_coord/n_seeds) over the per-seed dock tree.
  - 5TZ1 = the already-frozen `per_seed_fungal` in work/candidates/shortlist_confidence.json
    (NOT re-docked). 5FSA / 5V5Z = newly docked, work/candidates/dock_ens_{5FSA,5V5Z}/s*.
  - A conformer "coordinates" for a molecule iff n_coord >= MIN_COORD_SEEDS (3/5).
  - ensemble_min  = MIN conformer consensus over conformers where the molecule had ANY
                    coordinating pose (openafr.ensemble.aggregate, mode="min").
  - ensemble_mean = mean of the same (diagnostic; never a verdict).
  - n_conf_coord  = how many conformers coordinate (>= MIN_COORD_SEEDS).

Self-docking guard (a no-op for this ligand set, applied for honesty): a co-crystal azole is
never scored against the conformer it was solved in. None of the 16 shortlist molecules is
posaconazole (5FSA) or itraconazole (5V5Z), so nothing is excluded here.

Pre-committed per-candidate classification (on n_conf_coord over the conformers considered):
  ROBUST             coordinates in ALL considered conformers (the strongest hypothesis)
  CONSISTENT         coordinates in a majority but not all
  CONFORMER-FRAGILE  coordinates in exactly 1 (single-snapshot artifact; the informative demotion)
  NON-COORDINATING   coordinates in 0

Ensemble-consistency != potency: look #9 showed the ensemble does not lift within-azole ranking
to threshold, so ROBUST is a single-snapshot-artifact filter, not a binding/potency claim. The
CONFORMER-FRAGILE / NON-COORDINATING verdict is the informative demotion this axis surfaces.

Usage:
    python scripts/candidate_ensemble.py \
        [--dock-5fsa work/candidates/dock_ens_5FSA] \
        [--dock-5v5z work/candidates/dock_ens_5V5Z] \
        [--receptor-5fsa work/receptor_5FSA_A.pdb] \
        [--receptor-5v5z work/receptor_5V5Z_A.pdb] \
        [--baseline work/candidates/shortlist_confidence.json] [--json]
"""
import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr.ensemble import aggregate
from openafr.pdbqt import iron_position
# Reuse the frozen consensus logic unchanged — do not redefine it here.
from scripts.shortlist_confidence import FE_CUTOFF, consensus, per_seed_distances

MIN_COORD_SEEDS = 3     # coordinating seeds required for a conformer to count (frozen)

# Self-docking guard: co-crystal azole -> the conformer it was solved in (never scored there).
SELF_DOCK_EXCLUDE = {
    "posaconazole": {"5FSA"},
    "itraconazole": {"5V5Z"},
}


def classify(n_conf_coord, n_conf_considered):
    """Pre-committed ensemble-consistency verdict from the conformer coordination count."""
    if n_conf_considered == 0:
        return "UNRANKABLE"
    if n_conf_coord == 0:
        return "NON-COORDINATING"
    if n_conf_coord >= n_conf_considered:
        return "ROBUST"
    if n_conf_coord == 1:
        return "CONFORMER-FRAGILE"
    return "CONSISTENT"


def conformer_consensus(byseed):
    """(coordinating_min_or_None, n_coord, n_seeds) for one molecule in one conformer.

    Value is the consensus MIN over coordinating seeds when the molecule reached the iron in
    at least one seed, else None (so a non-coordinating conformer drops out of ensemble_min).
    """
    cons, _spread, n_coord, n_seeds = consensus(byseed)
    return (cons if n_coord > 0 else None), n_coord, n_seeds


def build_records(baseline_path, per_conf_byseed):
    """per_conf_byseed: {conformer_id: {name: {seed: best_N_Fe}}} for the docked conformers.

    5TZ1 is read from the frozen baseline JSON (per_seed_fungal), never re-docked.
    """
    baseline = json.load(open(baseline_path))
    tz1_byseed = {r["name"]: {k: v for k, v in r["per_seed_fungal"].items()}
                  for r in baseline["candidates"]}
    tz1_wt = {r["name"]: r["fungal_consensus"] for r in baseline["candidates"]}

    conformers = ["5TZ1"] + list(per_conf_byseed.keys())
    names = sorted(set(tz1_byseed))

    recs = []
    for name in names:
        skip = SELF_DOCK_EXCLUDE.get(name, set())
        per_conf = {}
        per_conf["5TZ1"] = conformer_consensus(tz1_byseed.get(name, {}))
        for cid, byname in per_conf_byseed.items():
            per_conf[cid] = conformer_consensus(byname.get(name, {}))

        considered = [c for c in conformers if c not in skip]
        vals = [per_conf[c][0] for c in considered]           # None-safe list for aggregation
        n_conf_coord = sum(1 for c in considered if per_conf[c][1] >= MIN_COORD_SEEDS)
        verdict = classify(n_conf_coord, len(considered))

        recs.append({
            "name": name,
            "wt_consensus": tz1_wt.get(name),
            "ensemble_min": _r(aggregate(vals, "min")),
            "ensemble_mean": _r(aggregate(vals, "mean")),
            "n_conf_coord": f"{n_conf_coord}/{len(considered)}",
            "verdict": verdict,
            "per_conformer": {
                c: {"consensus": _r(per_conf[c][0]),
                    "seeds_coordinating": f"{per_conf[c][1]}/{per_conf[c][2]}"}
                for c in considered
            },
            "excluded_conformers": sorted(skip),
        })
    # Sort by WT consensus (the shortlist order); unrankable last.
    recs.sort(key=lambda r: (r["wt_consensus"] is None, r["wt_consensus"] or 9e9))
    return recs


def _r(v):
    return None if v is None else round(v, 3)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dock-5fsa", default="work/candidates/dock_ens_5FSA")
    ap.add_argument("--dock-5v5z", default="work/candidates/dock_ens_5V5Z")
    ap.add_argument("--receptor-5fsa", default="work/receptor_5FSA_A.pdb")
    ap.add_argument("--receptor-5v5z", default="work/receptor_5V5Z_A.pdb")
    ap.add_argument("--baseline", default="work/candidates/shortlist_confidence.json")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    per_conf_byseed = {
        "5FSA": per_seed_distances(args.dock_5fsa, iron_position(args.receptor_5fsa)),
        "5V5Z": per_seed_distances(args.dock_5v5z, iron_position(args.receptor_5v5z)),
    }
    recs = build_records(args.baseline, per_conf_byseed)

    if args.json:
        print(json.dumps({"fe_cutoff": FE_CUTOFF, "min_coord_seeds": MIN_COORD_SEEDS,
                          "conformers": ["5TZ1", "5FSA", "5V5Z"], "candidates": recs},
                         indent=2))
        return 0

    print("=" * 100)
    print("CANDIDATE ENSEMBLE CONSISTENCY — C. albicans CYP51 {5TZ1, 5FSA, 5V5Z} (#70)")
    print(f"coordination cutoff {FE_CUTOFF} A   conformer counts if >= {MIN_COORD_SEEDS}/5 seeds")
    print("=" * 100)
    print("%-16s %8s %10s %10s %8s   %s" %
          ("candidate", "WT_min", "ens_min", "ens_mean", "coord", "verdict"))
    print("-" * 100)
    for r in recs:
        wm = "-" if r["wt_consensus"] is None else "%.3f" % r["wt_consensus"]
        em = "-" if r["ensemble_min"] is None else "%.3f" % r["ensemble_min"]
        me = "-" if r["ensemble_mean"] is None else "%.3f" % r["ensemble_mean"]
        print("%-16s %8s %10s %10s %8s   %s" %
              (r["name"], wm, em, me, r["n_conf_coord"], r["verdict"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
