"""External-holdout gate — implements work/PREREGISTRATION_external.md (issue #82).

Grades the FROZEN run-2 protocol, unchanged, on a never-scored active/decoy holdout
(data/ligands/external_actives.smi + external_decoys.smi), to test whether the geometric
criterion generalizes to data it has never touched — closing the overfitting objection.

Primary (gates): mode C — rank by minimum azole-nitrogen-to-heme-iron distance.
Secondary (report only): mode A — best Vina affinity (the run-2 worse-than-random baseline).

Endpoints AUC and EF@5% (issue #82); EF@1% printed for continuity but not an endpoint
(RESULTS_reliability). Pass bar: mode-C AUC >= 0.70 (protocol.yaml gate.min_auc), qualified
by the bootstrap 95% CI per the pre-registered three-outcome interpretation. Permutation test
and bootstrap CI come from the same openafr.scoring functions every prior gate uses.

Unrankable actives are ranked LAST, never dropped — the same anti-inflation rule as run 2.
Re-checks the pre-registration and protocol hashes before reporting.

Usage:
    python scripts/validate_gate_external.py [SCREEN_DIR] [RECEPTOR_PDB]
        SCREEN_DIR   dir of docked *.pdbqt poses   (default work/screen_external)
        RECEPTOR_PDB structure to read the iron from (default work/receptor_A.pdb)
"""
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import protocol
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses, read_scores
from openafr.scoring import bootstrap_auc_ci, permutation_p, report
from scripts.validate_gate2 import check_prereg

PREREG = "work/PREREGISTRATION_external.md"
PREREG_SHA = "66df34cb0fe4861bfffb7c209f26b88e318e53d29dd4037b6163c0aaa056f44c"
ACTIVES_FILE = "data/ligands/external_actives.smi"
DECOYS_FILE = "data/ligands/external_decoys.smi"


def _names(path):
    return {l.split("\t")[1].strip() for l in open(path) if "\t" in l}


def _ordered_flags(rows):
    """rows: (name, key, is_active), lower key better, None last. -> flags best-first."""
    ok = sorted((r for r in rows if r[1] is not None), key=lambda r: r[1])
    bad = [r for r in rows if r[1] is None]
    return [r[2] for r in ok + bad]


def classify_outcome(auc, ci_lo, min_auc):
    """The pre-registered three-outcome rule (work/PREREGISTRATION_external.md).

    FAIL              — point estimate below the bar; overfitting not excluded.
    PASS              — point estimate AND the whole bootstrap CI clear the bar.
    SUPPORTED-BUT-WIDE — point estimate clears it but the CI dips below.
    """
    if auc < min_auc:
        return "FAIL"
    return "PASS" if ci_lo >= min_auc else "SUPPORTED-BUT-WIDE"


def main():
    screen = sys.argv[1] if len(sys.argv) > 1 else "work/screen_external"
    receptor = sys.argv[2] if len(sys.argv) > 2 else "work/receptor_A.pdb"

    fe = iron_position(receptor)
    actives = _names(ACTIVES_FILE)
    expected = actives | _names(DECOYS_FILE)     # every frozen molecule, pose or not
    gate = protocol.load()["gate"]
    min_auc = gate["min_auc"]
    ef_fractions, bedroc_alpha = gate["ef_fractions"], gate["bedroc_alpha"]

    print("=" * 70)
    print("EXTERNAL HOLDOUT — PRE-REGISTERED GENERALIZATION TEST (issue #82)")
    print("=" * 70)
    print(f"receptor         : {receptor}")
    check_prereg(PREREG, PREREG_SHA)
    protocol.verify()
    print(f"held-out actives : {len(actives)}")
    print(f"pass bar         : mode-C AUC >= {min_auc}")

    # Seed rows from the frozen .smi lists, not just the poses on disk: a molecule that failed
    # to dock (no .pdbqt) must be ranked LAST as unrankable, never dropped — the pre-registered
    # anti-inflation rule. Dropping a failed active is the easiest way to fake a passing gate.
    dist_rows, score_rows, unranked = [], [], []
    for name in sorted(expected):
        is_act = name in actives
        p = os.path.join(screen, name + ".pdbqt")
        best_d, best_s = None, None
        if os.path.exists(p):
            sc = read_scores(p)
            best_s = sc[0] if sc else None
            for _num, atoms in read_poses(p):
                d = min_nitrogen_iron_distance(atoms, fe)
                if d is not None and (best_d is None or d < best_d):
                    best_d = d
        else:
            unranked.append(name)
        dist_rows.append((name, best_d, is_act))
        score_rows.append((name, best_s, is_act))

    if unranked:
        na = sum(1 for n in unranked if n in actives)
        print(f"no pose (docking failed): {len(unranked)} ({na} actives) — ranked LAST per pre-registration")

    rep = lambda label, rows, g=False: report(label, rows, g, ef_fractions, bedroc_alpha)
    c = rep("MODE C — rank by iron-approach distance", dist_rows, g=True)
    rep("MODE A — rank by Vina score", score_rows)

    flags = _ordered_flags(dist_rows)
    pval = permutation_p(flags)
    lo, hi = bootstrap_auc_ci(flags)
    positions = [i + 1 for i, fl in enumerate(flags) if fl]

    print("\nStatistical support (mode C):")
    print(f"  permutation test (20,000 shuffles): p = {pval:.4f}")
    print(f"  bootstrap 95% CI over actives     : AUC {lo:.3f} - {hi:.3f}")
    print(f"  active rank positions (of {len(flags)}) : {positions}")

    auc, efs = c[0], c[1]
    ef5 = efs[ef_fractions.index(0.05)]
    outcome = classify_outcome(auc, lo, min_auc)
    print("\n" + "=" * 70)
    print(f"MODE C: AUC {auc:.3f} (need >= {min_auc})   EF@5% {ef5:.2f}x   "
          f"CI [{lo:.3f}, {hi:.3f}]")
    if outcome == "FAIL":
        print("OUTCOME: FAIL — the criterion does NOT generalize to a never-scored set.")
        print("Per pre-registration: do NOT re-tune the criterion to recover a pass.")
        return 1
    if outcome == "PASS":
        print("OUTCOME: PASS — generalizes; the whole 95% CI clears the bar (overfitting excluded).")
    else:
        print("OUTCOME: SUPPORTED-BUT-WIDE — point estimate passes, CI lower bound below the bar.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
