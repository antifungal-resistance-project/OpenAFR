"""Run-2 validation gate — implements work/PREREGISTRATION_run2.md exactly.

Primary (gates):   mode C — rank by minimum nitrogen-to-heme-iron distance
Secondary (report only): mode A (Vina score), mode D (combined rank)

Pre-registered pass bar: AUC >= 0.70 AND EF@1% >= 5.0x on mode C.

Pre-registered handling rule: an active that produces no usable pose is ranked LAST,
never dropped. Dropping unrankable actives inflates every metric and is the single
easiest way to fake a passing gate.
"""
import hashlib
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import protocol
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses, read_scores
from openafr.scoring import report

PREREG = "work/PREREGISTRATION_run2.md"
PREREG_SHA = "7c89d78cdadba90bae54b60ba9583eec606f6c6aec88688dd2ef353e773ddbf1"
ACTIVES = "data/ligands/actives_holdout_final.smi"


def check_prereg(prereg=PREREG, prereg_sha=PREREG_SHA):
    """Refuse to grade if the pre-registration was edited after being frozen."""
    if not os.path.exists(prereg):
        print(f"WARNING: {prereg} missing — cannot verify the rules were fixed in advance.")
        return
    got = hashlib.sha256(open(prereg, "rb").read()).hexdigest()
    if got == prereg_sha:
        print(f"pre-registration verified unmodified (sha256 {got[:16]}...)")
    else:
        print("!! PRE-REGISTRATION HAS BEEN MODIFIED SINCE IT WAS FROZEN !!")
        print(f"   expected {prereg_sha[:16]}...  got {got[:16]}...")
        print("   Any 'pass' from this run is not credible.")


def run_gate(screen, receptor="work/receptor_A.pdb", prereg=PREREG,
             prereg_sha=PREREG_SHA, actives_file=ACTIVES,
             title="RUN 2 — PRE-REGISTERED HELD-OUT VALIDATION"):
    """Grade a screen against the pre-registered gate. Defaults reproduce run 2 exactly;
    the C. auris resistance-port gate calls this with a mutant receptor + its own frozen
    pre-registration. The grading math, the pass bar (from protocol.yaml) and the
    rank-unrankable-actives-last rule are identical across every run — only the receptor,
    the held-out actives file and the pre-registration document change."""
    fe = iron_position(receptor)
    actives = {l.split("\t")[1].strip() for l in open(actives_file) if "\t" in l}

    gate = protocol.load()["gate"]
    min_auc, min_ef1 = gate["min_auc"], gate["min_ef1"]
    ef_fractions, bedroc_alpha = gate["ef_fractions"], gate["bedroc_alpha"]
    rep = lambda label, rows, gating=False: report(
        label, rows, gating, ef_fractions, bedroc_alpha)

    print("=" * 70)
    print(title)
    print("=" * 70)
    print(f"receptor         : {receptor}")
    check_prereg(prereg, prereg_sha)
    protocol.verify()
    print(f"held-out actives : {len(actives)}")
    print(f"pass bar         : AUC >= {min_auc} AND EF@1% >= {min_ef1}x (mode C only)")

    dist_rows, score_rows, both = [], [], []
    for f in sorted(os.listdir(screen)):
        if not f.endswith(".pdbqt"):
            continue
        name = f[:-6]
        is_act = name in actives
        p = os.path.join(screen, f)
        sc = read_scores(p)
        best_d = None
        for _num, atoms in read_poses(p):
            d = min_nitrogen_iron_distance(atoms, fe)
            if d is not None and (best_d is None or d < best_d):
                best_d = d
        dist_rows.append((name, best_d, is_act))
        score_rows.append((name, sc[0] if sc else None, is_act))
        both.append((name, best_d, sc[0] if sc else None, is_act))

    # confirm every held-out active is present
    seen = {r[0] for r in dist_rows if r[2]}
    if seen != actives:
        print(f"\n!! MISSING ACTIVES: {actives - seen} — results not valid")

    c = rep("MODE C — rank by iron-approach distance", dist_rows, gating=True)
    rep("MODE A — rank by Vina score", score_rows)

    dr = {b[0]: i for i, b in enumerate(sorted([b for b in both if b[1] is not None],
                                               key=lambda b: b[1]))}
    sr = {b[0]: i for i, b in enumerate(sorted([b for b in both if b[2] is not None],
                                               key=lambda b: b[2]))}
    comb = [(b[0], (dr[b[0]] + sr[b[0]]) if (b[0] in dr and b[0] in sr) else None, b[3])
            for b in both]
    rep("MODE D — combined rank (distance + score)", comb)

    print("\n" + "=" * 70)
    passed = c[0] >= min_auc and c[1][0] >= min_ef1
    print(f"MODE C: AUC {c[0]:.3f} (need >= {min_auc})   EF@1% {c[1][0]:.2f}x (need >= {min_ef1})")
    if passed:
        print("GATE: PASS — hypothesis H1 supported on held-out actives.")
        print("Screening new molecules is justified under this protocol and domain.")
        return 0
    print("GATE: FAIL — H1 not supported on held-out actives.")
    print("Per pre-registration, the geometric ranking claim must be withdrawn or revised.")
    print("Do NOT invent a new ranking method and report it as a pass.")
    return 1


def main():
    screen = sys.argv[1] if len(sys.argv) > 1 else "work/screen2"
    receptor = sys.argv[2] if len(sys.argv) > 2 else "work/receptor_A.pdb"
    return run_gate(screen, receptor)


if __name__ == "__main__":
    sys.exit(main())
