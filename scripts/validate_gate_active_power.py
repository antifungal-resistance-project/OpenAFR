"""Active-power gate — implements work/PREREGISTRATION_active_power.md exactly.

The Run-2 primary gate passed on 7 held-out azoles vs 348 property-matched decoys at
mode C AUC 0.794, but with a bootstrap 95% CI of 0.590-0.946 whose lower bound sits BELOW
the 0.70 bar (preprint Limitation #1). This gate re-runs the identical contrast — same
rigid 5TZ1, same mode C (minimum nitrogen-to-heme-iron distance), same property-matched
decoy pool, same frozen protocol/seed — with the ACTIVE side widened from 7 to a blind
N=300 sample of measured-active azoles, so the CI becomes tight and its LOWER BOUND can be
compared against the bar.

Primary (gates):   actives (N=300) vs the 350 property-matched decoys, mode C.
                   PASS bar (protocol.yaml): AUC point estimate >= 0.70.
                   Limitation #1 is RETIRED only if the bootstrap 95% CI lower bound >= 0.70.
Secondary (report only):
                   S1 product-claim contrast — the same actives vs the NON-azole verified
                      inactives (novel-chemotype triage; comparators 0.810 / 0.823).
                   S2 permutation p, bootstrap CI, the tip.
                   S3 EF@1/5/10%, BEDROC — non-gating (seed-fragile per preprint §3.2).
                   S4 domain-edge: AUC in heavy-atom halves (size-artifact check).

Pre-registered handling rule (from Run 2): an active with no usable pose is ranked LAST,
never dropped. EF is NOT part of the gate here (Run-2's EF@1% 12.68x did not survive
5-seed resampling; preprint §3.2).

Usage:
    python scripts/validate_gate_active_power.py [SCREEN_DIR] [IRON_RECEPTOR]
"""
import hashlib
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import protocol
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses
from openafr.scoring import report, bootstrap_auc_ci, permutation_p
from scripts.validate_gate_verified import azole_bearing

PREREG = "work/PREREGISTRATION_active_power.md"
PREREG_SHA = "9327563ac69ff087f17af7e79fb6ee8245dc8ff5fb7f8759a2a6ed9099d77d02"
SAMPLE = "data/ligands/active_power_sample.smi"
SAMPLE_SHA = "ba8d1f3cb05717017ababa321c8a74b4e0122c5eae94c91ceab92915dd3c37ae"
DECOYS = "data/ligands/decoys_holdout.smi"
INACTIVES = "data/ligands/verified_inactives.smi"


def _sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def check_frozen(path, expected):
    if not os.path.exists(path):
        print(f"WARNING: {path} missing — cannot verify it was fixed in advance.")
        return
    got = _sha(path)
    if got == expected:
        print(f"  frozen OK  {path}  (sha256 {got[:16]}...)")
    else:
        print(f"  !! MODIFIED SINCE FREEZE: {path}")
        print(f"     expected {expected[:16]}...  got {got[:16]}...  — any pass is not credible.")


def _names(smi_path):
    """Molecule names as prep_ligands writes them: the tab-delimited second column."""
    return {l.split("\t")[1].strip() for l in open(smi_path) if "\t" in l}


def _best_distance(screen, fe):
    """name -> best (minimum) nitrogen-to-iron distance over all poses, or None."""
    out = {}
    for f in sorted(os.listdir(screen)):
        if not f.endswith(".pdbqt"):
            continue
        best = None
        for _num, atoms in read_poses(os.path.join(screen, f)):
            d = min_nitrogen_iron_distance(atoms, fe)
            if d is not None and (best is None or d < best):
                best = d
        out[f[:-6]] = best
    return out


def _rows(actives, others, dist):
    """(name, key, is_active) for the actives + a given decoy/inactive pool. A member with
    no docked pose is included with key None so report() ranks it LAST, never dropped."""
    rows = []
    for name in sorted(actives | others):
        rows.append((name, dist.get(name), name in actives))
    return rows


def _ordered_flags(rows):
    ok = sorted((r for r in rows if r[1] is not None), key=lambda r: r[1])
    bad = [r for r in rows if r[1] is None]
    return [r[2] for r in ok + bad]


def _graded(label, rows, ef_fractions, bedroc_alpha, gating=False):
    report(label, rows, gating, ef_fractions, bedroc_alpha)
    flags = _ordered_flags(rows)
    lo, hi = bootstrap_auc_ci(flags)
    p = permutation_p(flags)
    print(f"  bootstrap 95% CI : {lo:.3f} - {hi:.3f}")
    print(f"  permutation p    : {p:.4g}  (20,000 shuffles)")
    return lo, hi


def _domain_edge(actives, others, dist, ef_fractions, bedroc_alpha):
    """S4 — split actives at the median heavy-atom count and report AUC in each half vs the
    full decoy pool, to confirm the signal is not carried by molecular size alone."""
    try:
        from rdkit import Chem
    except ImportError:
        print("\nS4 domain-edge: rdkit unavailable — skipped.")
        return
    heavy = {}
    for smi in (SAMPLE,):
        for l in open(smi):
            if "\t" not in l:
                continue
            s, name = l.split("\t")[0], l.split("\t")[1].strip()
            m = Chem.MolFromSmiles(s)
            if m is not None:
                heavy[name] = m.GetNumHeavyAtoms()
    sized = sorted((h, n) for n, h in heavy.items() if n in actives)
    if not sized:
        return
    med = sized[len(sized) // 2][0]
    small = {n for h, n in sized if h <= med}
    large = {n for h, n in sized if h > med}
    print(f"\nS4 domain-edge (median heavy-atoms = {med}):")
    for tag, half in (("<= median", small), ("> median", large)):
        if half:
            report(f"  actives {tag} ({len(half)}) vs decoys", _rows(half, others, dist),
                   False, ef_fractions, bedroc_alpha)


def run_gate(screen="work/screen_active_power", receptor="work/receptor_A.pdb"):
    print("=" * 70)
    print("ACTIVE-POWER GATE — WIDEN THE HELD-OUT ACTIVES TO RETIRE THE n=7 TIP")
    print("=" * 70)
    print("freeze checks:")
    check_frozen(PREREG, PREREG_SHA)
    check_frozen(SAMPLE, SAMPLE_SHA)
    protocol.verify()

    gate = protocol.load()["gate"]
    min_auc = gate["min_auc"]
    ef_fractions, bedroc_alpha = gate["ef_fractions"], gate["bedroc_alpha"]

    fe = iron_position(receptor)
    dist = _best_distance(screen, fe)

    actives = _names(SAMPLE)
    pmdecoys = _names(DECOYS)
    seen_act = {n for n in actives if n in dist}
    print(f"\nreceptor (iron)  : {receptor}")
    print(f"actives          : {len(actives)}  (docked: {len(seen_act)})")
    print(f"prop-matched dec : {len(pmdecoys)}  (docked: {len(pmdecoys & set(dist))})")
    if seen_act != actives:
        missing = actives - seen_act
        print(f"  note: {len(missing)} actives have no pose — ranked LAST per pre-registration.")
    print(f"pass bar         : AUC point >= {min_auc}; Limitation #1 RETIRED iff CI lower >= {min_auc}")

    # PRIMARY — actives vs property-matched decoys
    prim = _rows(actives, pmdecoys, dist)
    lo, hi = _graded("PRIMARY — actives (N=300) vs 350 property-matched decoys, mode C",
                     prim, ef_fractions, bedroc_alpha, gating=True)
    from openafr.scoring import metrics
    auc = metrics(_ordered_flags(prim))[0]

    # SECONDARY S1 — product-claim contrast (vs NON-azole verified inactives)
    if os.path.exists(INACTIVES):
        inact = _names(INACTIVES)
        azo = azole_bearing(INACTIVES)
        non_azole = inact - azo
        docked_non_azole = non_azole & set(dist)
        if docked_non_azole:
            _graded(f"S1 PRODUCT — actives vs NON-azole verified inactives "
                    f"({len(non_azole)}; docked {len(docked_non_azole)}) [comparators 0.810/0.823]",
                    _rows(actives, non_azole, dist), ef_fractions, bedroc_alpha)
        else:
            print("\nS1 product contrast: no NON-azole verified inactives docked in this "
                  "screen — skipped (dock verified_inactives.smi to enable).")

    # S4 — size-artifact check
    _domain_edge(actives, pmdecoys, dist, ef_fractions, bedroc_alpha)

    print("\n" + "=" * 70)
    print(f"PRIMARY mode C: AUC {auc:.3f} (need >= {min_auc})   "
          f"bootstrap 95% CI {lo:.3f} - {hi:.3f}")
    if auc < min_auc:
        print("GATE: FAIL — the Run-2 0.794 point estimate does NOT hold on a blind N=300 "
              "sample.\nThis is a material negative about the headline claim; report it as "
              "prominently as a pass (pre-registered FAIL branch).")
        return 1
    if lo >= min_auc:
        print("GATE: PASS — RETIRED. With n=7 sampling noise removed, the mode-C AUC clears "
              "the\n0.70 bar with its LOWER BOUND, not just its point estimate. Preprint §6.1 "
              "updates\nfrom 'a true AUC under the bar cannot be excluded' to this tight interval.")
        return 0
    print("GATE: PASS (point) but MARGINAL — the point estimate clears 0.70 yet the CI lower "
          f"bound\n({lo:.3f}) still sits below it even at N=300. Limitation #1 is QUANTIFIED, "
          "not retired;\nreport the narrowed interval honestly.")
    return 0


def main():
    screen = sys.argv[1] if len(sys.argv) > 1 else "work/screen_active_power"
    receptor = sys.argv[2] if len(sys.argv) > 2 else "work/receptor_A.pdb"
    return run_gate(screen, receptor)


if __name__ == "__main__":
    sys.exit(main())
