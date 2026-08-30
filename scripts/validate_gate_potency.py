"""Continuous-potency gate — implements work/PREREGISTRATION_potency.md (#81).

Does the N–Fe geometry track *how potent* a molecule is, not merely whether it is labelled
active? This grader joins the frozen potency datasets (median pChEMBL per compound) to the
docked closest-N–Fe distance of the SAME compound, and reports the pre-registered Spearman
correlation with a one-sided permutation p and a bootstrap CI. It refuses to certify anything
if the pre-registration, the protocol, the receptor anchor or either dataset moved since freeze
(the same anti-tuning guarantee as every other gate here).

H1 (the only thing that gates): over the ENZYME in-domain set, Spearman ρ(closest N–Fe,
median pChEMBL) is negative — tighter iron approach ↔ higher potency.
PASS bar: ρ ≤ −0.30 AND one-sided permutation p < 0.05.

Usage (needs the docked screens; see the header of work/RESULTS_potency.md for the commands):
    conda run -n openafr python scripts/validate_gate_potency.py \
        --enzyme-screen work/potency_screen_enzyme \
        --wholecell-screen work/potency_screen_wholecell
"""
import argparse
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from openafr import potency, protocol
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses

PREREG = "work/PREREGISTRATION_potency.md"
PREREG_SHA = "5d3416d9907ffa9c2e3f0a66d43a121340d92e1eccf8958012ba0752aaede34a"
SHAFILE = "work/PREREG_potency.sha256"
RECEPTOR = "work/receptor_A.pdb"
ENZYME_TSV = "data/ligands/potency_enzyme.tsv"
WHOLECELL_TSV = "data/ligands/potency_wholecell.tsv"

FE_CUTOFF = 3.0          # closest N–Fe ≤ this Å = coordinates the iron (frozen protocol)
RHO_BAR = -0.30          # pre-registered H1 pass bar (direction + magnitude)
P_BAR = 0.05
N_PERM = 100000
N_BOOT = 10000
SEED = 42


def _sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def verify_freeze():
    """Re-check every hash pinned in work/PREREG_potency.sha256; refuse to certify on drift."""
    ok = True
    got = _sha(PREREG)
    print(f"pre-registration : {PREREG}")
    if got == PREREG_SHA:
        print(f"  verified unmodified (sha256 {got[:16]}...)")
    else:
        ok = False
        print("  !! PRE-REGISTRATION MODIFIED SINCE FREEZE — any pass is not credible")
        print(f"     expected {PREREG_SHA[:16]}...  got {got[:16]}...")
    for line in open(SHAFILE):
        want, _, path = line.strip().partition("  ")
        if not path or path == PREREG:
            continue
        got = _sha(path)
        mark = "ok" if got == want else "!! MOVED"
        if got != want:
            ok = False
        print(f"  {mark:9} {path}  ({got[:12]}...)")
    protocol.verify()
    return ok


def load_pchembl(tsv):
    """frozen .tsv -> ({cid: median_pchembl}, {cid: verdict})."""
    pch, verdict = {}, {}
    with open(tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(header)}
        for line in fh:
            r = line.rstrip("\n").split("\t")
            if r[ci["in_domain"]] != "True":
                continue
            cid = r[ci["cid"]]
            pch[cid] = float(r[ci["pchembl"]])
            verdict[cid] = r[ci["verdict"]]
    return pch, verdict


def closest_n_fe(screen, fe):
    """screen dir of docked *.pdbqt -> {cid: min N–Fe distance over its poses}."""
    out = {}
    if not os.path.isdir(screen):
        print(f"  (screen dir {screen} not present — skipping this set)")
        return out
    for f in sorted(os.listdir(screen)):
        if not f.endswith(".pdbqt"):
            continue
        cid = f[:-6]
        best = None
        for _num, atoms in read_poses(os.path.join(screen, f)):
            d = min_nitrogen_iron_distance(atoms, fe)
            if d is not None and (best is None or d < best):
                best = d
        if best is not None:
            out[cid] = best
    return out


def correlate(label, cids, dist, pch, *, gate=False):
    """Report Spearman ρ, one-sided permutation p and bootstrap CI for one subset."""
    xs = [dist[c] for c in cids]
    ys = [pch[c] for c in cids]
    n = len(cids)
    print(f"\n{label}  (n={n})")
    if n < 8:
        print(f"  too few compounds to correlate (n={n} < 8) — reported, not gated")
        return None
    rho = potency.spearman(xs, ys)
    p = potency.permutation_p(xs, ys, n_perm=N_PERM, seed=SEED, alternative="less")
    lo, hi = potency.bootstrap_ci(xs, ys, n_boot=N_BOOT, seed=SEED)
    print(f"  Spearman ρ = {rho:+.3f}   95% CI [{lo:+.3f}, {hi:+.3f}]   "
          f"one-sided perm p = {p:.4g}")
    if gate:
        passed = rho <= RHO_BAR and p < P_BAR
        print(f"  H1 bar: ρ ≤ {RHO_BAR} AND p < {P_BAR}  ->  "
              f"{'PASS' if passed else 'FAIL'}")
        return passed
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--enzyme-screen", default="work/potency_screen_enzyme")
    ap.add_argument("--wholecell-screen", default="work/potency_screen_wholecell")
    args = ap.parse_args()

    print("=" * 72)
    print("CONTINUOUS-POTENCY GATE (#81) — geometry vs measured pChEMBL")
    print("=" * 72)
    clean = verify_freeze()
    if not clean:
        print("\nRefusing to certify: a frozen input moved. Fix the drift and re-run.")

    fe = iron_position(RECEPTOR)
    print(f"\nreceptor iron    : ({fe[0]:.3f}, {fe[1]:.3f}, {fe[2]:.3f})")

    # ---- PRIMARY: enzyme (CHEMBL1780) ----
    e_pch, e_verdict = load_pchembl(ENZYME_TSV)
    e_dist = closest_n_fe(args.enzyme_screen, fe)
    e_cids = sorted(set(e_pch) & set(e_dist))
    missing = sorted(set(e_pch) - set(e_dist))
    if missing:
        print(f"\nWARNING: {len(missing)} enzyme compounds have no docked pose: {missing}")

    print("\n" + "-" * 72)
    print("PRIMARY — enzyme (CHEMBL1780), on-target binding potency")
    print("-" * 72)
    passed = correlate("H1: full in-domain set", e_cids, e_dist, e_pch, gate=True)
    novel = [c for c in e_cids if e_verdict[c] == "in-envelope-novel"]
    correlate("leakage-guarded: near-known azoles excluded", novel, e_dist, e_pch)
    bound = [c for c in e_cids if e_dist[c] <= FE_CUTOFF]
    correlate(f"iron-bound only (N–Fe ≤ {FE_CUTOFF} Å)", bound, e_dist, e_pch)

    # ---- SECONDARY: whole-cell ----
    w_pch, _ = load_pchembl(WHOLECELL_TSV)
    w_dist = closest_n_fe(args.wholecell_screen, fe)
    w_cids = sorted(set(w_pch) & set(w_dist))
    print("\n" + "-" * 72)
    print("SECONDARY — whole-cell C. albicans potency (confounded; does NOT gate)")
    print("-" * 72)
    correlate("whole-cell in-domain set", w_cids, w_dist, w_pch)

    print("\n" + "=" * 72)
    if passed and clean:
        print("GATE: PASS — H1 supported: geometry tracks on-target potency (see ρ/CI above).")
        return 0
    print("GATE: FAIL — H1 not supported at the pre-registered bar.")
    print("Per pre-registration: geometry enriches actives but does not demonstrably rank")
    print("continuous on-target potency. Do not subset-hunt these poses for a pass.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
