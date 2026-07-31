"""Run-2 validation gate — implements work/PREREGISTRATION_run2.md exactly.

Primary (gates):   mode C — rank by minimum nitrogen-to-heme-iron distance
Secondary (report only): mode A (Vina score), mode D (combined rank)

Pre-registered pass bar: AUC >= 0.70 AND EF@1% >= 5.0x on mode C.

Pre-registered handling rule: an active that produces no usable pose is ranked LAST,
never dropped. Dropping unrankable actives inflates every metric and is the single
easiest way to fake a passing gate.
"""
import hashlib
import math
import os
import sys

from rdkit.ML.Scoring import Scoring

MIN_AUC = 0.70
MIN_EF1 = 5.0
EF_FRACTIONS = [0.01, 0.05, 0.10]
BEDROC_ALPHA = 20.0
PREREG = "work/PREREGISTRATION_run2.md"
PREREG_SHA = "7c89d78cdadba90bae54b60ba9583eec606f6c6aec88688dd2ef353e773ddbf1"


def check_prereg():
    """Refuse to grade if the pre-registration was edited after being frozen."""
    if not os.path.exists(PREREG):
        print(f"WARNING: {PREREG} missing — cannot verify the rules were fixed in advance.")
        return
    got = hashlib.sha256(open(PREREG, "rb").read()).hexdigest()
    if got == PREREG_SHA:
        print(f"pre-registration verified unmodified (sha256 {got[:16]}...)")
    else:
        print("!! PRE-REGISTRATION HAS BEEN MODIFIED SINCE IT WAS FROZEN !!")
        print(f"   expected {PREREG_SHA[:16]}...  got {got[:16]}...")
        print("   Any 'pass' from this run is not credible.")


def poses(path):
    cur = None
    for l in open(path):
        if l.startswith("MODEL"):
            cur = []
        elif l.startswith("ENDMDL"):
            yield cur
        elif l.startswith(("ATOM", "HETATM")):
            nm = l[12:16].strip()
            el = (l[76:78].strip() or nm[0]).upper()
            if el != "H":
                cur.append((nm, float(l[30:38]), float(l[38:46]), float(l[46:54])))


def vina_scores(path):
    return [float(l.split()[3]) for l in open(path) if "REMARK VINA RESULT" in l]


def iron(receptor):
    for l in open(receptor):
        if l.startswith("HETATM") and l[76:78].strip().upper() == "FE":
            return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    raise SystemExit("no heme iron in receptor")


def metrics(ordered_flags):
    r = [[1 if a else 0] for a in ordered_flags]
    return (Scoring.CalcAUC(r, 0),
            [Scoring.CalcEnrichment(r, 0, [f])[0] for f in EF_FRACTIONS],
            Scoring.CalcBEDROC(r, 0, BEDROC_ALPHA))


def report(label, ranked, gating=False):
    """ranked: list of (name, key, is_active), lower key = better. None key = ranked last."""
    ok = [r for r in ranked if r[1] is not None]
    bad = [r for r in ranked if r[1] is None]
    ordered = sorted(ok, key=lambda r: r[1]) + bad          # unrankable go last
    flags = [r[2] for r in ordered]
    auc, efs, bedroc = metrics(flags)
    n_act = sum(flags)
    print(f"\n{label}{'   [PRIMARY — gates the result]' if gating else '   (secondary, reported only)'}")
    print(f"  ranked           : {len(ordered)}  ({n_act} actives, {len(ordered)-n_act} decoys)")
    if bad:
        nb = sum(1 for r in bad if r[2])
        print(f"  unrankable       : {len(bad)} molecules ({nb} actives) — placed LAST per pre-registration")
    print(f"  AUC              : {auc:.3f}")
    for f, ef in zip(EF_FRACTIONS, efs):
        print(f"  EF@{int(f*100):>3}%          : {ef:.2f}x")
    print(f"  BEDROC (a=20)    : {bedroc:.3f}")
    top = [r[0] for r in ordered[:15] if r[2]]
    print(f"  actives in top 15: {len(top)}/{n_act}  {top if top else ''}")
    return auc, efs, bedroc


def main():
    screen = sys.argv[1] if len(sys.argv) > 1 else "work/screen2"
    receptor = sys.argv[2] if len(sys.argv) > 2 else "work/receptor_A.pdb"
    fe = iron(receptor)
    actives = {l.split("\t")[1].strip()
               for l in open("data/ligands/actives_holdout_final.smi") if "\t" in l}

    print("=" * 70)
    print("RUN 2 — PRE-REGISTERED HELD-OUT VALIDATION")
    print("=" * 70)
    check_prereg()
    print(f"held-out actives : {len(actives)}")
    print(f"pass bar         : AUC >= {MIN_AUC} AND EF@1% >= {MIN_EF1}x (mode C only)")

    dist_rows, score_rows, both = [], [], []
    for f in sorted(os.listdir(screen)):
        if not f.endswith(".pdbqt"):
            continue
        name = f[:-6]
        is_act = name in actives
        p = os.path.join(screen, f)
        sc = vina_scores(p)
        best_d = None
        for pose in poses(p):
            ns = [a for a in pose if a[0].startswith("N")]
            if not ns:
                continue
            d = min(math.dist(a[1:], fe) for a in ns)
            if best_d is None or d < best_d:
                best_d = d
        dist_rows.append((name, best_d, is_act))
        score_rows.append((name, sc[0] if sc else None, is_act))
        both.append((name, best_d, sc[0] if sc else None, is_act))

    # confirm every held-out active is present
    seen = {r[0] for r in dist_rows if r[2]}
    if seen != actives:
        print(f"\n!! MISSING ACTIVES: {actives - seen} — results not valid")

    c = report("MODE C — rank by iron-approach distance", dist_rows, gating=True)
    a = report("MODE A — rank by Vina score", score_rows)

    dr = {b[0]: i for i, b in enumerate(sorted([b for b in both if b[1] is not None],
                                               key=lambda b: b[1]))}
    sr = {b[0]: i for i, b in enumerate(sorted([b for b in both if b[2] is not None],
                                               key=lambda b: b[2]))}
    comb = [(b[0], (dr[b[0]] + sr[b[0]]) if (b[0] in dr and b[0] in sr) else None, b[3])
            for b in both]
    report("MODE D — combined rank (distance + score)", comb)

    print("\n" + "=" * 70)
    passed = c[0] >= MIN_AUC and c[1][0] >= MIN_EF1
    print(f"MODE C: AUC {c[0]:.3f} (need >= {MIN_AUC})   EF@1% {c[1][0]:.2f}x (need >= {MIN_EF1})")
    if passed:
        print("GATE: PASS — hypothesis H1 supported on held-out actives.")
        print("Screening new molecules is justified under this protocol and domain.")
        return 0
    print("GATE: FAIL — H1 not supported on held-out actives.")
    print("Per pre-registration, the geometric ranking claim must be withdrawn or revised.")
    print("Do NOT invent a new ranking method and report it as a pass.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
