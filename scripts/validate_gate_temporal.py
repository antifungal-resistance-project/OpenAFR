"""Temporal-holdout gate — implements work/PREREGISTRATION_temporal.md (issue #82, stronger
follow-up to the disjoint-population holdout #108).

Grades the FROZEN run-2 protocol, unchanged, on a publication-date holdout
(data/ligands/temporal_actives.smi + temporal_decoys.smi): the newest-literature azoles
(earliest antifungal / CYP51 document year >= 2024) that no gate has ever docked, versus
property-matched decoys. This closes the one honest residual the #108 holdout left on the
record ("a clean never-touched holdout, but not a later time slice") — the criterion was
defined on the classic approved azoles, and this tests whether it ranks recent research
chemistry it never saw. The grading math, the pass bar and the three-outcome interpretation
are reused verbatim from the external gate (#108); only the held-out set changes.

Primary (gates): mode C — rank by minimum azole-nitrogen-to-heme-iron distance.
Secondary (report only): mode A — best Vina affinity (the run-2 worse-than-random baseline).

Usage:
    python scripts/validate_gate_temporal.py [SCREEN_DIR] [RECEPTOR_PDB]
        SCREEN_DIR   dir of docked *.pdbqt poses   (default work/screen_temporal)
        RECEPTOR_PDB structure to read the iron from (default work/receptor_A.pdb)
"""
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import protocol
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses, read_scores
from openafr.scoring import bootstrap_auc_ci, permutation_p, report
from scripts.make_temporal_holdout import ACTIVES_OUT, CUTOFF_YEAR, DECOYS_OUT
from scripts.validate_gate2 import check_prereg
from scripts.validate_gate_external import _names, _ordered_flags, classify_outcome

PREREG = "work/PREREGISTRATION_temporal.md"
PREREG_SHA = "33ff95584aac263bc498eb9f5cda8707ca0dbd42e14b460f7d42288f3363c558"


def main():
    screen = sys.argv[1] if len(sys.argv) > 1 else "work/screen_temporal"
    receptor = sys.argv[2] if len(sys.argv) > 2 else "work/receptor_A.pdb"

    print("=" * 70)
    print("TEMPORAL HOLDOUT — PRE-REGISTERED PUBLICATION-DATE TEST (issue #82)")
    print("=" * 70)
    print(f"receptor         : {receptor}")
    check_prereg(PREREG, PREREG_SHA)

    if not (os.path.exists(ACTIVES_OUT) and os.path.exists(DECOYS_OUT)):
        sys.exit(f"no temporal set — build it first:\n"
                 f"  python scripts/make_temporal_holdout.py actives\n"
                 f"  python scripts/make_temporal_holdout.py decoys")
    if not os.path.isdir(screen) or not any(f.endswith(".pdbqt") for f in os.listdir(screen)):
        sys.exit(f"no docked poses in {screen} — dock the set:\n"
                 f"  RECEPTOR=work/receptor.pdbqt scripts/screen.sh work/temporal_ligands {screen}")

    fe = iron_position(receptor)
    actives = _names(ACTIVES_OUT)
    expected = actives | _names(DECOYS_OUT)
    gate = protocol.load()["gate"]
    min_auc = gate["min_auc"]
    ef_fractions, bedroc_alpha = gate["ef_fractions"], gate["bedroc_alpha"]

    protocol.verify()
    print(f"held-out actives : {len(actives)}   (earliest antifungal year > {CUTOFF_YEAR})")
    print(f"pass bar         : mode-C AUC >= {min_auc}")

    # Seed rows from the frozen .smi lists, not just the poses on disk: a molecule that failed
    # to dock (no .pdbqt) is ranked LAST as unrankable, never dropped — the anti-inflation rule.
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
    print(f"MODE C: AUC {auc:.3f} (need >= {min_auc})   EF@5% {ef5:.2f}x   CI [{lo:.3f}, {hi:.3f}]")
    if outcome == "FAIL":
        print("OUTCOME: FAIL — the criterion does NOT generalize to newer-literature azoles.")
        print("Per pre-registration: do NOT re-tune the criterion to recover a pass.")
        return 1
    if outcome == "PASS":
        print("OUTCOME: PASS — generalizes across time; the whole 95% CI clears the bar.")
    else:
        print("OUTCOME: SUPPORTED-BUT-WIDE — point estimate passes, CI lower bound below the bar.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
