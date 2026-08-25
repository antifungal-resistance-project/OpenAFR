"""Metal-aware baseline gate — implements work/PREREGISTRATION_baselines.md.

Answers the one objection run 2 leaves open: the iron-approach criterion beats *Vina's*
score, but Vina is known-broken on metal sites — does it also beat a competent, metal-aware
scorer? GNINA (CNN rescoring) is that opponent. Same held-out set, same poses, one new
number to sort by.

Three rankings, all on the SAME run-2 held-out poses:

    C  (ours)     iron filter (N-Fe < 3.0 A), then best Vina affinity among survivors
                  — the run-2 PRIMARY criterion, reproduced here
    G  (opponent) best GNINA CNNaffinity over ALL poses — no mechanistic prior
    GC (probe)    best GNINA CNNaffinity among the iron-bound poses only

Headline comparison G vs C, endpoints AUC and EF@5% (EF@1% is excluded — RESULTS_reliability
established it is not robust). Interpretation is pre-committed in the pre-registration; this
script only reports which outcome occurred. It does NOT re-open run 2.

Grading math, the ef_fractions and BEDROC alpha all come from the frozen protocol.yaml, the
same source every other gate reads. Unrankable actives (no rankable pose) are ranked LAST,
never dropped. This script re-checks the pre-registration SHA-256 and the protocol hash
before reporting.

Usage:
    python scripts/rescore_gnina.py [GNINA_DIR] [POSE_DIR] [RECEPTOR]
        GNINA_DIR  dir of per-ligand GNINA --score_only outputs (<name>.gnina)
                                                        (default work/gnina_scores)
        POSE_DIR   dir of the matching Vina .pdbqt poses  (default work/screen2)
        RECEPTOR   structure to read the heme iron from    (default work/receptor_A.pdb)
"""
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import protocol
from openafr.gnina import best_cnnaffinity, parse_score_only
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses, read_scores
from openafr.scoring import report
from scripts.validate_gate2 import ACTIVES, check_prereg

PREREG = "work/PREREGISTRATION_baselines.md"
PREREG_SHA = "374ba516b83b73b7e2907a9dae7166d2abddc3713132870b5a1313c11a5994d7"
FE_CUTOFF = 3.0  # angstroms; the same metal-coordination cutoff every gate uses


def molecule_keys(pose_path, gnina_path, fe):
    """Return (ours_C, gnina_G, gnina_GC) ranking keys for one molecule (lower = better).

    ours_C  : best (most negative) Vina affinity among iron-bound poses, or None.
    gnina_G : -(best CNNaffinity over ALL poses), or None.
    gnina_GC: -(best CNNaffinity among iron-bound poses), or None.

    Poses and GNINA blocks are aligned 1:1 by input-file order. If they do not align
    (different pose counts) the molecule's GNINA keys are None — never guessed — so a
    parsing mismatch degrades honestly rather than silently mis-scoring.
    """
    poses = list(read_poses(pose_path))
    vina = read_scores(pose_path)
    iron_bound = [min_nitrogen_iron_distance(atoms, fe) for _num, atoms in poses]
    keep = [d is not None and d < FE_CUTOFF for d in iron_bound]

    # Mode C — reproduce the run-2 primary: min Vina score among iron-bound poses.
    ours = None
    for score, k in zip(vina, keep):
        if k and (ours is None or score < ours):
            ours = score

    gnina_G = gnina_GC = None
    if gnina_path and os.path.exists(gnina_path):
        blocks = parse_score_only(open(gnina_path).read())
        if len(blocks) == len(poses) and blocks:
            g = best_cnnaffinity(blocks)
            gc = best_cnnaffinity(blocks, keep=keep)
            gnina_G = None if g is None else -g
            gnina_GC = None if gc is None else -gc
        else:
            print(f"  !! pose/GNINA count mismatch for {os.path.basename(pose_path)}: "
                  f"{len(poses)} poses vs {len(blocks)} GNINA blocks — GNINA excluded for it")
    return ours, gnina_G, gnina_GC


def _verdict(label, our, opp, higher_is=None):
    """Pre-committed G-vs-C readout on one endpoint (both are 'higher = better' metrics)."""
    if our > opp:
        return f"{label}: C {our:.3f} > G {opp:.3f}  -> geometry ahead"
    if opp > our:
        return f"{label}: G {opp:.3f} > C {our:.3f}  -> GNINA ahead"
    return f"{label}: C {our:.3f} = G {opp:.3f}  -> tie"


def main():
    gnina_dir = sys.argv[1] if len(sys.argv) > 1 else "work/gnina_scores"
    pose_dir = sys.argv[2] if len(sys.argv) > 2 else "work/screen2"
    receptor = sys.argv[3] if len(sys.argv) > 3 else "work/receptor_A.pdb"

    if not os.path.isdir(pose_dir):
        sys.exit(f"no pose dir: {pose_dir} — run scripts/run_screen2.sh first")
    fe = iron_position(receptor)
    actives = {l.split("\t")[1].strip() for l in open(ACTIVES) if "\t" in l}

    gate = protocol.load()["gate"]
    min_auc = gate["min_auc"]
    rep = lambda label, rows, gating=False: report(
        label, rows, gating, gate["ef_fractions"], gate["bedroc_alpha"])
    ef5_idx = gate["ef_fractions"].index(0.05)

    print("=" * 70)
    print("METAL-AWARE BASELINE GATE — GNINA vs the iron-approach criterion")
    print("=" * 70)
    print(f"receptor         : {receptor}")
    print(f"heme iron        : ({fe[0]:.2f}, {fe[1]:.2f}, {fe[2]:.2f})")
    print(f"metal cutoff     : N-Fe < {FE_CUTOFF} A")
    print(f"GNINA outputs    : {gnina_dir}")
    check_prereg(PREREG, PREREG_SHA)
    protocol.verify()

    have_gnina = os.path.isdir(gnina_dir)
    if not have_gnina:
        print(f"\nNOTE: {gnina_dir} absent — GNINA not yet run. Reproducing mode C only;")
        print("      run scripts/run_gnina_rescore.sh on a GNINA host to fill G/GC.")

    c_rows, g_rows, gc_rows = [], [], []
    for f in sorted(os.listdir(pose_dir)):
        if not f.endswith(".pdbqt"):
            continue
        name = f[:-6]
        is_act = name in actives
        gpath = os.path.join(gnina_dir, name + ".gnina") if have_gnina else None
        ours, g, gc = molecule_keys(os.path.join(pose_dir, f), gpath, fe)
        # Append every molecule with its (possibly None) key; scoring.report() ranks
        # None-key rows LAST rather than dropping them — the pre-registration's
        # "unrankable ranked last, never dropped" rule, enforced identically to run 2.
        c_rows.append((name, ours, is_act))
        g_rows.append((name, g, is_act))
        gc_rows.append((name, gc, is_act))

    seen = {r[0] for r in c_rows if r[2]}
    if seen != actives:
        print(f"\n!! actives missing from mode C: {actives - seen} — check the pose set")

    c = rep("MODE C (OURS) — iron filter, then Vina affinity  [run-2 primary, reproduced]",
            c_rows, gating=True)
    if not any(r[1] is not None for r in g_rows):
        print("\nGNINA modes not available — comparison incomplete. Nothing headlined.")
        return 2
    g = rep("MODE G (GNINA) — best CNNaffinity, all poses  [PRIMARY OPPONENT]", g_rows)
    gc = rep("MODE GC (GNINA) — best CNNaffinity among iron-bound poses  [probe]", gc_rows)

    print("\n" + "=" * 70)
    print("PRE-COMMITTED COMPARISON  (see work/PREREGISTRATION_baselines.md)")
    print("  headline: G vs C at AUC and EF@5%  (EF@1% excluded — not robust)")
    print("  " + _verdict("AUC  ", c[0], g[0]))
    print("  " + _verdict("EF@5%", c[1][ef5_idx], g[1][ef5_idx]))
    print(f"  (probe) GC vs C: AUC {gc[0]:.3f} vs {c[0]:.3f} | "
          f"EF@5% {gc[1][ef5_idx]:.2f}x vs {c[1][ef5_idx]:.2f}x")

    c_wins = c[0] > g[0] and c[1][ef5_idx] > g[1][ef5_idx]
    g_wins = g[0] > c[0] and g[1][ef5_idx] > c[1][ef5_idx]
    print()
    if c_wins:
        print("OUTCOME: C > G on both endpoints — geometry beats a metal-aware learned scorer,")
        print("not just broken vanilla Vina. Run-2 framing stands and strengthens.")
    elif g_wins:
        print("OUTCOME: G > C on both endpoints — a metal-aware scorer wins. Per the")
        print("pre-registration, retire 'geometry beats the scoring function' as the headline;")
        print("lead with the self-validating gate. The criterion remains a cheap first-pass filter.")
    else:
        print("OUTCOME: mixed / tie on the endpoints — report as a TIE: the free, interpretable")
        print("criterion matches the learned scorer. Confirm against bootstrap CIs before")
        print("calling any endpoint a win (prereg: overlap = tie).")
    print(f"\n(point estimates only; n={len(actives)} actives — CIs are wide. AUC bar for "
          f"reference: {min_auc}.)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
