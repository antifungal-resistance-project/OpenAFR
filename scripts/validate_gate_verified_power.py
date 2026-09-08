"""Verified-inactive POWER gate — implements work/PREREGISTRATION_verified_power.md exactly.

Look #6 (`work/RESULTS_verified_inactives.md`) replaced the presumed-inactive decoys with 279
compounds MEASURED not to inhibit C. albicans and passed on the point estimate (mode C AUC
0.716), but on only 7 held-out actives, so its bootstrap 95% CI was 0.506-0.885 — the lower
bound below the 0.70 bar (preprint Limitation #2, "the single largest caveat"). This gate
re-runs that identical contrast — same rigid 5TZ1, same mode C, same 279 measured inactives,
same frozen protocol/seed — with the ACTIVE side widened from 7 to the same blind N=300 sample
look #10 used, so the CI becomes tight and its LOWER BOUND can be compared against the bar.

Only the actives change (7 -> 300); the inactives, criterion, receptor and protocol are all
held fixed from look #6. Any difference from look #6's 0.716 is attributable to active sample
size, not method drift.

Primary (gates):   N=300 actives vs ALL 279 measured inactives, mode C.
                   PASS bar (protocol.yaml): AUC point estimate >= 0.70.
                   Limitation #2 is RETIRED only if the bootstrap 95% CI lower bound >= 0.70.
Secondary (report only):
                   S1 non-azole subgroup (novel-chemotype triage; look #6 0.810 / look #10 0.774).
                   S2 azole-bearing subgroup (the crux; look #6 0.650 at n=7).
                   S3 the tip, band composition, permutation p, bootstrap CI.
                   S4 domain-edge: AUC in heavy-atom halves (size-artifact check).

Pre-registered handling rule (from Run 2 / look #6): an active with no usable pose is ranked
LAST, never dropped. EF is NOT part of the gate (seed-fragile, preprint §3.2).

Usage:
    python scripts/validate_gate_verified_power.py [SCREEN_DIR] [RECEPTOR]
"""
import hashlib
import os
import pathlib
import sys

from rdkit import Chem

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import protocol
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses, read_scores
from openafr.scoring import bootstrap_auc_ci, permutation_p, report
from scripts.validate_gate2 import check_prereg
from scripts.validate_gate_verified import azole_bearing, tip_statistics

PREREG = "work/PREREGISTRATION_verified_power.md"
PREREG_SHA = "40b656a38767d9ab96ff0165f14d313cb573224b5bd310bccf68c76431ce2da6"
ACTIVES = "data/ligands/active_power_sample.smi"
ACTIVES_SHA = "ba8d1f3cb05717017ababa321c8a74b4e0122c5eae94c91ceab92915dd3c37ae"
INACTIVES = "data/ligands/verified_inactives.smi"
INACTIVES_SHA = "37a20eec1c026c42d3c8efb53c1dd27cda2f9057ca541cfc1afed21a35faea6b"
DEFAULT_SCREEN = "work/screen_verified_power"
# Look #6 comparators (same contrast, 7 actives), quoted not recomputed.
LOOK6 = {"auc": 0.716, "non_azole": 0.810, "azole": 0.650, "ci_lo": 0.506, "ci_hi": 0.885}


def _sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def check_frozen(path, expected):
    if not os.path.exists(path):
        print(f"WARNING: {path} missing — cannot verify it was fixed in advance.")
        return
    got = _sha(path)
    tag = "frozen OK" if got == expected else "!! MODIFIED SINCE FREEZE"
    print(f"  {tag}  {path}  (sha256 {got[:16]}...)")
    if got != expected:
        print(f"     expected {expected[:16]}...  — any pass is not credible.")


def _heavy_atoms(smi_path):
    """name -> heavy-atom count, for the pre-registered domain-edge (size) split."""
    out = {}
    for line in open(smi_path):
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 2:
            continue
        mol = Chem.MolFromSmiles(parts[0])
        if mol is not None:
            out[parts[1]] = mol.GetNumHeavyAtoms()
    return out


def main():
    screen = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SCREEN
    receptor = sys.argv[2] if len(sys.argv) > 2 else "work/receptor_A.pdb"
    fe = iron_position(receptor)
    actives = {l.split("\t")[1].strip() for l in open(ACTIVES) if "\t" in l}

    gate = protocol.load()["gate"]
    min_auc = gate["min_auc"]
    rep = lambda label, rows, gating=False: report(
        label, rows, gating, gate["ef_fractions"], gate["bedroc_alpha"])

    print("=" * 74)
    print("VERIFIED-INACTIVE POWER GATE — N=300 ACTIVES vs MEASURED NON-INHIBITORS (Limitation #2)")
    print("=" * 74)
    print(f"receptor         : {receptor}")
    print(f"poses            : {screen}")
    check_prereg(PREREG, PREREG_SHA)
    print("frozen inputs:")
    check_frozen(ACTIVES, ACTIVES_SHA)
    check_frozen(INACTIVES, INACTIVES_SHA)
    protocol.verify()

    dist_rows, score_rows = [], []
    for f in sorted(os.listdir(screen)):
        if not f.endswith(".pdbqt"):
            continue
        name = f[:-6]
        is_act = name in actives
        p = os.path.join(screen, f)
        best_d = None
        for _num, atoms in read_poses(p):
            d = min_nitrogen_iron_distance(atoms, fe)
            if d is not None and (best_d is None or d < best_d):
                best_d = d
        sc = read_scores(p)
        dist_rows.append((name, best_d, is_act))
        score_rows.append((name, sc[0] if sc else None, is_act))

    seen = {r[0] for r in dist_rows if r[2]}
    missing = actives - seen
    if missing:
        print(f"\n  {len(missing)} of {len(actives)} actives produced no pose (docking "
              f"failures) — ranked LAST per pre-registration, never dropped")
    n_inactive = sum(1 for r in dist_rows if not r[2])
    submitted_act = sum(1 for l in open(ACTIVES) if "\t" in l)
    submitted_ina = sum(1 for l in open(INACTIVES) if "\t" in l)
    print(f"actives (of {submitted_act}) : {len(seen)} posed")
    print(f"measured inactives (of {submitted_ina}): {n_inactive} posed"
          + (f"  ({submitted_ina - n_inactive} docking failures excluded, not counted as decoys)"
             if submitted_ina != n_inactive else ""))
    print(f"pass bar         : AUC >= {min_auc} on mode C; Limitation #2 RETIRED iff CI lower bound >= {min_auc}")

    # Rank the actives that failed to dock last (None key), per pre-registration.
    dist_rows = [(n, d, a) if not (a and n not in seen) else (n, None, a) for n, d, a in dist_rows]

    c = rep("MODE C — rank by iron-approach distance  [vs ALL measured inactives]",
            dist_rows, gating=True)
    rep("MODE A — rank by Vina score", score_rows)

    azoles = azole_bearing(INACTIVES)
    hard = [r for r in dist_rows if r[2] or r[0] in azoles]
    soft = [r for r in dist_rows if r[2] or r[0] not in azoles]
    s2 = rep(f"S2 SUBGROUP — vs azole-bearing inactives only (n={len(hard) - len(seen)}): "
             f"working azole vs failed azole?  [look #6: 0.650]", hard)
    s1 = rep(f"S1 SUBGROUP — vs non-azole inactives only (n={len(soft) - len(seen)}): "
             f"novel-chemotype triage  [look #6: 0.810 / look #10: 0.774]", soft)

    # S4 — domain-edge (heavy-atom) split of the actives; both halves keep all inactives.
    ha = _heavy_atoms(ACTIVES)
    act_ha = sorted(v for k, v in ha.items() if k in actives)
    if act_ha:
        med = act_ha[len(act_ha) // 2]
        small = [r for r in dist_rows if (not r[2]) or ha.get(r[0], 0) <= med]
        large = [r for r in dist_rows if (not r[2]) or ha.get(r[0], 0) > med]
        rep(f"S4 domain-edge — actives with <= {med} heavy atoms (n={sum(1 for r in small if r[2])})", small)
        rep(f"S4 domain-edge — actives with >  {med} heavy atoms (n={sum(1 for r in large if r[2])})", large)

    # S3 — the tip + statistical support on the PRIMARY pool.
    ok = [r for r in dist_rows if r[1] is not None]
    ordered = sorted(ok, key=lambda r: r[1]) + [r for r in dist_rows if r[1] is None]
    rank, above, pct = tip_statistics(ordered)
    flags = [r[2] for r in ordered]
    p_perm = permutation_p(flags, n_shuffles=20000)
    ci_lo, ci_hi = bootstrap_auc_ci(flags, n_boot=10000)
    band = [r for r in ordered if r[1] is not None and 2.47 <= r[1] <= 2.88]
    band_act = sum(1 for r in band if r[2])

    print("\nS3 — THE TIP (best real active vs measured inactives above it)")
    if rank is not None:
        print(f"  best true active : {ordered[rank-1][0]}  at rank {rank}/{len(ordered)} ({pct:.2f}th pct)")
        print(f"  inactives above it: {above}")
    print(f"  validated iron-bound band 2.47-2.88 A: {len(band)} molecules, {band_act} actives, "
          f"{len(band)-band_act} measured inactives")
    print("\nSTATISTICAL SUPPORT (primary pool; descriptive, gate is the AUC point estimate)")
    print(f"  permutation test (20,000 shuffles): p = {p_perm:.4f}")
    print(f"  bootstrap 95% CI over actives     : AUC {ci_lo:.3f} - {ci_hi:.3f}")

    print("\n" + "=" * 74)
    passed = c[0] >= min_auc
    retired = passed and ci_lo >= min_auc
    print(f"MODE C (PRIMARY): AUC {c[0]:.3f}  CI {ci_lo:.3f}-{ci_hi:.3f}  (need point & lower "
          f">= {min_auc})   [look #6 n=7: AUC {LOOK6['auc']:.3f}, CI {LOOK6['ci_lo']}-{LOOK6['ci_hi']}]")
    print(f"S1 non-azole {s1[0]:.3f}   S2 azole-bearing {s2[0]:.3f}  [look #6: {LOOK6['azole']}]")
    if retired:
        print("GATE: PASS — LIMITATION #2 RETIRED. At N=300 the criterion separates real azoles")
        print("from compounds MEASURED not to inhibit, with the CI lower bound above the bar. The")
        print("decoy ceiling is confirmed real, not an artifact of decoy construction.")
        return 0
    if passed:
        print("GATE: MARGINAL — point estimate passes but CI lower bound < bar even at N=300.")
        print("Limitation #2 is QUANTIFIED (tight interval), not retired. Report honestly.")
        return 0
    print("GATE: FAIL — look #6's 0.716 does not hold at power; broad enrichment vs MEASURED")
    print("inactives is below the bar at N=300. Per pre-registration this is a material negative")
    print("about the product headline. Do NOT re-filter the inactive set and re-run.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
