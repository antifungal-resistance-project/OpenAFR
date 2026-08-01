"""The validation gate: can this pipeline re-discover the drugs we already know work?

Pipeline stage: validate/

This is the load-bearing module of the whole project. It answers one question:

    When the 8 known azoles are hidden in a haystack of 399 property-matched decoys,
    does the pipeline rank the real drugs near the top?

If it cannot, nothing downstream is trustworthy and the gate exits NONZERO so no
screen can run on top of it.

Two ranking modes are reported side by side, because the difference is the finding:

    score-only   : rank by Vina affinity alone
    constrained  : discard any pose whose azole nitrogen is not within 3 A of the
                   heme iron, then rank the survivors

Metrics (all from rdkit.ML.Scoring, the canonical implementations):
    AUC     - overall ranking power. 0.5 = coin flip, 1.0 = perfect.
    EF@x%   - enrichment factor: how many times more actives appear in the top x%
              than random chance would give. EF 10 at 1% = 10x better than random.
    BEDROC  - like AUC but weights the top of the list heavily (alpha=20), which is
              what matters when you can only afford to test the top few hundred.
"""
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses, read_scores
from openafr.scoring import EF_FRACTIONS, metrics

FE_CUTOFF = 3.0

# Gate thresholds, declared BEFORE looking at results (pre-registration).
# A pipeline that cannot beat these has not earned the right to screen new molecules.
MIN_AUC = 0.70
MIN_EF1 = 5.0


def summarize(pdbqt, fe):
    """Return (best_score, best_iron_bound_score or None)."""
    scores = read_scores(pdbqt)
    if not scores:
        return None, None
    best_constrained = None
    for (num, atoms), score in zip(read_poses(pdbqt), scores):
        d = min_nitrogen_iron_distance(atoms, fe)
        if d is not None and d < FE_CUTOFF:
            if best_constrained is None or score < best_constrained:
                best_constrained = score
    return scores[0], best_constrained


def report(label, rows, n_actives_total):
    """rows: (name, score, is_active) for molecules that produced a usable score."""
    rows = sorted(rows, key=lambda r: r[1])          # more negative = better
    n_act = sum(r[2] for r in rows)
    auc, efs, bedroc = metrics([r[2] for r in rows])
    print(f"\n{label}")
    print(f"  molecules ranked : {len(rows)}  ({n_act} actives, {len(rows)-n_act} decoys)")
    if n_act < n_actives_total:
        print(f"  NOTE: {n_actives_total - n_act} active(s) produced no usable pose "
              f"and are excluded — this LOWERS the apparent difficulty.")
    print(f"  AUC              : {auc:.3f}   (0.5 = random)")
    for f, ef in zip(EF_FRACTIONS, efs):
        print(f"  EF@{int(f*100):>3}%          : {ef:.2f}x random")
    print(f"  BEDROC (a=20)    : {bedroc:.3f}")
    print(f"  actives in top 10: {sum(r[2] for r in rows[:10])}/{n_act}")
    return auc, efs, bedroc, rows


def main():
    screen_dir = sys.argv[1] if len(sys.argv) > 1 else "work/screen"
    receptor = sys.argv[2] if len(sys.argv) > 2 else "work/receptor_A.pdb"
    fe = iron_position(receptor)

    active_names = {l.split("\t")[1].strip() for l in open("data/ligands/actives.smi")
                    if "\t" in l}

    score_rows, constrained_rows = [], []
    missing_active = []
    for f in sorted(os.listdir(screen_dir)):
        if not f.endswith(".pdbqt"):
            continue
        name = f[:-6]
        is_active = name in active_names
        best, constrained = summarize(os.path.join(screen_dir, f), fe)
        if best is not None:
            score_rows.append((name, best, is_active))
        if constrained is not None:
            constrained_rows.append((name, constrained, is_active))
        elif is_active:
            missing_active.append(name)

    n_act = len(active_names)
    print("=" * 68)
    print("VALIDATION GATE — can the pipeline re-discover the known azoles?")
    print("=" * 68)
    print(f"receptor        : {receptor}")
    print(f"heme iron       : ({fe[0]:.2f}, {fe[1]:.2f}, {fe[2]:.2f})")
    print(f"metal cutoff    : N-Fe < {FE_CUTOFF} A")
    print(f"pre-registered  : AUC >= {MIN_AUC}, EF@1% >= {MIN_EF1}x")

    a1 = report("MODE A — rank by docking score only", score_rows, n_act)
    a2 = report("MODE B — metal constraint, then rank by score", constrained_rows, n_act)

    if missing_active:
        print(f"\n  actives with NO iron-bound pose (excluded from mode B): "
              f"{', '.join(missing_active)}")

    print("\n" + "=" * 68)
    best_auc = max(a1[0], a2[0])
    best_ef1 = max(a1[1][0], a2[1][0])
    passed = best_auc >= MIN_AUC and best_ef1 >= MIN_EF1
    print(f"best AUC {best_auc:.3f} (need >= {MIN_AUC})   "
          f"best EF@1% {best_ef1:.2f}x (need >= {MIN_EF1})")
    if passed:
        print("GATE: PASS — the pipeline recovers known drugs; screening is justified.")
        return 0
    print("GATE: FAIL — the pipeline cannot recover known drugs.")
    print("Screening new molecules on top of this would produce untrustworthy rankings.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
