"""Consensus validation gate — implements work/PREREGISTRATION_consensus.md.

Tests whether a *variance-reduced* iron-approach criterion ranks the run-2 held-out
actives better (or at least as well) as the single-search criterion, on the SAME blind
wild-type held-out gate that validated the original.

Consensus criterion (fixed in the pre-registration, chosen on principle not tuned):
    for each molecule, dock it under R independent Vina searches (different seeds; box,
    exhaustiveness and num_modes unchanged from the frozen protocol), take the per-seed
    best N-Fe distance, and rank by the MINIMUM across seeds — a less-biased estimate of
    the achievable coordination, which is what the criterion always meant to measure. The
    run-to-run spread (max - min) is reported as a per-molecule confidence signal.

The single-seed-42 slice is graded alongside as the apples-to-apples baseline (seed 42 is
the frozen protocol seed, so that slice reproduces the original single-search gate).

Pass bar and grading math are identical to run 2 (from protocol.yaml): AUC >= 0.70 AND
EF@1% >= 5.0x on mode C. Unrankable actives are ranked LAST, never dropped.

Usage:
    python scripts/validate_gate_consensus.py [CONSENSUS_DIR] [RECEPTOR]
        CONSENSUS_DIR  parent dir with per-seed subdirs s<seed>/  (default work/screen_consensus_wt)
        RECEPTOR       structure to read the heme iron from       (default work/receptor_A.pdb)
"""
import glob
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import protocol
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses
from openafr.scoring import report
from scripts.validate_gate2 import ACTIVES, check_prereg

PREREG = "work/PREREGISTRATION_consensus.md"
PREREG_SHA = "1c6ea042ab85a2ac567b8095003c0a2a3c6340af68adc66c531ad2e33cc40230"
BASELINE_SEED = "42"  # the frozen protocol seed


def consensus_stats(byseed):
    """(consensus_min, spread) from a {seed: best_N_Fe} dict.

    consensus = min over the seeds that produced a coordinating pose (a less-biased
    estimate of the achievable minimum); spread = max - min over those seeds (0 if fewer
    than two). Seeds with no coordinating pose (None) are ignored, never counted as
    coordinating; if no seed coordinated, consensus is None (ranked last by the gate)."""
    vals = [v for v in byseed.values() if v is not None]
    if not vals:
        return None, 0.0
    return min(vals), (max(vals) - min(vals) if len(vals) > 1 else 0.0)


def best_distance(pose_path, fe):
    """Minimum nitrogen-to-iron distance over all poses in one docked file, or None."""
    best = None
    for _num, atoms in read_poses(pose_path):
        d = min_nitrogen_iron_distance(atoms, fe)
        if d is not None and (best is None or d < best):
            best = d
    return best


def main():
    parent = sys.argv[1] if len(sys.argv) > 1 else "work/screen_consensus_wt"
    receptor = sys.argv[2] if len(sys.argv) > 2 else "work/receptor_A.pdb"
    fe = iron_position(receptor)
    actives = {l.split("\t")[1].strip() for l in open(ACTIVES) if "\t" in l}

    # per-seed pose dirs are s<digits>/; exclude screen.sh's s<seed>_pdbqt/ intermediate
    # dirs (which hold pre-docking input poses, not docked results).
    seed_dirs = sorted(d for d in glob.glob(os.path.join(parent, "s*"))
                       if os.path.isdir(d) and os.path.basename(d)[1:].isdigit())
    if not seed_dirs:
        sys.exit(f"no per-seed subdirs (s<seed>/) under {parent} — run scripts/screen_consensus.sh")
    seeds = [os.path.basename(d)[1:] for d in seed_dirs]

    gate = protocol.load()["gate"]
    min_auc, min_ef1 = gate["min_auc"], gate["min_ef1"]
    rep = lambda label, rows, gating=False: report(
        label, rows, gating, gate["ef_fractions"], gate["bedroc_alpha"])

    print("=" * 70)
    print("CONSENSUS GATE — PRE-REGISTERED WILD-TYPE HELD-OUT RE-VALIDATION")
    print("=" * 70)
    print(f"receptor         : {receptor}")
    print(f"seeds            : {seeds}  (baseline slice = seed {BASELINE_SEED})")
    check_prereg(PREREG, PREREG_SHA)
    protocol.verify()

    # molecule -> {seed: best N-Fe}; a molecule missing from a seed dir simply has fewer
    # replicates (never silently treated as coordinating).
    per_mol = {}
    for d in seed_dirs:
        seed = os.path.basename(d)[1:]
        for f in sorted(os.listdir(d)):
            if not f.endswith(".pdbqt"):
                continue
            name = f[:-6]
            per_mol.setdefault(name, {})[seed] = best_distance(os.path.join(d, f), fe)

    cons_rows, base_rows, spreads = [], [], {}
    for name, byseed in per_mol.items():
        cons, spreads[name] = consensus_stats(byseed)
        cons_rows.append((name, cons, name in actives))
        base_rows.append((name, byseed.get(BASELINE_SEED), name in actives))

    seen = {r[0] for r in cons_rows if r[2]}
    if seen != actives:
        print(f"\n!! MISSING ACTIVES: {actives - seen} — results not valid")

    c = rep("MODE C (CONSENSUS) — min N-Fe over seeds  [PRIMARY — gates]", cons_rows, gating=True)
    b = rep(f"MODE C (BASELINE, seed {BASELINE_SEED}) — single-search min N-Fe", base_rows)

    print("\nheld-out actives: consensus vs single-seed baseline")
    print("  %-16s %10s %8s %10s" % ("active", "consensus", "spread", f"seed{BASELINE_SEED}"))
    for a in sorted(actives, key=lambda a: dict((r[0], r[1]) for r in cons_rows).get(a) or 9e9):
        cd = dict((r[0], r[1]) for r in cons_rows).get(a)
        bd = dict((r[0], r[1]) for r in base_rows).get(a)
        print("  %-16s %10s %8s %10s" % (
            a, "-" if cd is None else "%.3f" % cd,
            "%.3f" % spreads.get(a, 0.0),
            "-" if bd is None else "%.3f" % bd))

    print("\n" + "=" * 70)
    passed = c[0] >= min_auc and c[1][0] >= min_ef1
    print(f"CONSENSUS mode C: AUC {c[0]:.3f} (need >= {min_auc})   EF@1% {c[1][0]:.2f}x (need >= {min_ef1})")
    print(f"BASELINE  mode C: AUC {b[0]:.3f}                 EF@1% {b[1][0]:.2f}x")
    if passed:
        print("GATE: PASS — consensus criterion supported on the wild-type held-out set.")
        return 0
    print("GATE: FAIL — consensus criterion not supported on the held-out set.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
