"""Verified-inactive gate — implements work/PREREGISTRATION_verified_inactives.md.

The 6th look on the ledger, and the FIRST on an independent dataset. Every previous look
(run-2 baseline, Y132F transfer, consensus, reliability, axial) ranked the same 7 held-out
azoles against the same 348 *presumed*-inactive decoys. The stopping rule in the preprint
(§5) closed that dataset: "a further attempt requires an independent dataset — new actives
or a new decoy construction — not another feature on this one." This is the new decoy
construction: compounds MEASURED not to inhibit C. albicans (openafr/inactives.py).

Deliberately NOT a new criterion. The ranking key is mode C — minimum nitrogen-to-heme-iron
distance — exactly as in run 2, under the same frozen protocol, same seed 42, same receptor,
same 7 held-out actives. The only thing that changes anywhere in this run is where the
non-actives came from. That makes the run-2 baseline the direct comparator and any
difference attributable to decoy provenance rather than to method drift.

Two questions, both fixed in advance:
    H1 (gates)      does the criterion separate the actives from MEASURED inactives at all?
                    AUC >= 0.70 (protocol min_auc). EF is reported, never gated: with 7
                    actives it is a coarse ladder (preprint §4.2 recommendation 2).
    H2 (reported)   the tip. On presumed decoys, six outranked the best true active. The
                    pool-size-neutral statistic is the best active's PERCENTILE rank; the
                    raw count of inactives above it is reported alongside.

Two subgroups are pre-specified (not post-hoc slices): 58.4% of the verified inactives carry
an azole ring, because the only compounds anyone measures against C. albicans are antifungal
med-chem programs. The azole-bearing half asks the sharpest form of the question — can pose
geometry tell a working azole from a failed azole analogue? Both subgroups are reported with
AUC; neither gates.

Usage:
    python scripts/validate_gate_verified.py [SCREEN_DIR] [RECEPTOR]
"""
import os
import pathlib
import sys

from rdkit import Chem

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import protocol
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses, read_scores
from openafr.scoring import bootstrap_auc_ci, permutation_p, report
from scripts.validate_gate2 import ACTIVES, check_prereg

PREREG = "work/PREREGISTRATION_verified_inactives.md"
PREREG_SHA = "d2926f776149690689063da2ffa59235a41ad9bfafee860641ba0012cf2e2701"
DEFAULT_SCREEN = "work/screen_verified"
INACTIVES = "data/ligands/verified_inactives.smi"
# The azole warhead: imidazole / 1,2,4- and 1,2,3-triazole / tetrazole rings.
AZOLE_SMARTS = [Chem.MolFromSmarts(s) for s in ("c1ncnn1", "c1cncn1", "c1ccnn1", "c1cnnn1")]
# The run-2 baseline: same seed, same receptor, same actives, PRESUMED-inactive decoys.
# Published in work/RESULTS_run2_holdout.md; quoted here as the comparator, not recomputed.
# Its best active ranked 3rd of 355 (active ranks 3, 8, 29, 43, 49, 124, 274), so two
# presumed decoys outranked it -- the tip statistic this run reproduces on measured data.
RUN2 = {"auc": 0.794, "ef1": 12.68, "n_decoys": 348, "n_total": 355,
        "best_active_rank": 3, "above_best_active": 2}


def tip_statistics(ordered):
    """ordered: rows (name, key, is_active) already sorted best-first.

    Returns (rank of the best active, 1-based; how many inactives outrank it; its percentile
    among all ranked molecules). The percentile is the pool-size-neutral form — a raw count
    of decoys above the best active is not comparable between pools of different sizes.
    """
    for i, row in enumerate(ordered):
        if row[2]:
            return i + 1, i, (i + 1) / len(ordered) * 100.0
    return None, None, None


def azole_bearing(smi_path):
    """Names in a .smi whose molecule carries an azole ring. Pre-specified subgroup split."""
    out = set()
    for line in open(smi_path):
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 2:
            continue
        mol = Chem.MolFromSmiles(parts[0])
        if mol is not None and any(mol.HasSubstructMatch(q) for q in AZOLE_SMARTS):
            out.add(parts[1])
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

    print("=" * 70)
    print("VERIFIED-INACTIVE GATE — INDEPENDENT DECOY CONSTRUCTION (MEASURED NON-INHIBITORS)")
    print("=" * 70)
    print(f"receptor         : {receptor}")
    print(f"poses            : {screen}")
    check_prereg(PREREG, PREREG_SHA)
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
    if seen != actives:
        print(f"\n!! MISSING ACTIVES: {actives - seen} — results not valid")
    n_inactive = sum(1 for r in dist_rows if not r[2])
    print(f"held-out actives : {len(seen)}")
    print(f"verified inactives: {n_inactive}")
    # Docking failures are reported, never silently absorbed: a molecule that produced no
    # pose file is missing from the pool entirely, which shrinks the denominator.
    if os.path.exists(INACTIVES):
        submitted = sum(1 for l in open(INACTIVES) if "\t" in l)
        if submitted != n_inactive:
            print(f"  !! {submitted - n_inactive} of {submitted} inactives produced no pose "
                  f"(docking failures) — excluded from the pool, not counted as decoys")
    print(f"pass bar         : AUC >= {min_auc} on mode C (EF reported, not gating)")

    c = rep("MODE C — rank by iron-approach distance  [vs MEASURED inactives]",
            dist_rows, gating=True)
    rep("MODE A — rank by Vina score", score_rows)

    ok = [r for r in dist_rows if r[1] is not None]
    ordered = sorted(ok, key=lambda r: r[1]) + [r for r in dist_rows if r[1] is None]
    rank, above, pct = tip_statistics(ordered)

    if os.path.exists(INACTIVES):
        azoles = azole_bearing(INACTIVES)
        hard = [r for r in dist_rows if r[2] or r[0] in azoles]
        soft = [r for r in dist_rows if r[2] or r[0] not in azoles]
        rep(f"SUBGROUP (pre-specified) — vs azole-bearing inactives only "
            f"(n={len(hard) - len(seen)}): can geometry tell a working azole from a failed one?",
            hard)
        rep(f"SUBGROUP (pre-specified) — vs non-azole inactives only "
            f"(n={len(soft) - len(seen)})", soft)

    print("\nH2 — THE TIP (the question the presumed-inactive set could not answer)")
    if rank is None:
        print("  no active was rankable — H2 undefined")
    else:
        print(f"  best true active : {ordered[rank-1][0]}  at rank {rank}/{len(ordered)}"
              f"  ({pct:.2f}th percentile)")
        print(f"  inactives above it: {above}")
    print(f"  comparator (run 2, PRESUMED decoys): best active rank "
          f"{RUN2['best_active_rank']}/{RUN2['n_total']} "
          f"({RUN2['best_active_rank'] / RUN2['n_total'] * 100:.2f}th percentile), "
          f"{RUN2['above_best_active']} decoys above it")
    print("  top 10 by iron-approach distance:")
    for i, r in enumerate(ordered[:10], 1):
        print(f"    {i:2d}. {r[0]:<44} {r[1]:.3f} A   {'ACTIVE' if r[2] else 'inactive'}")

    # The same statistical support run 2 reported, so the two runs are readable side by
    # side. Descriptive only -- the gate is the AUC point estimate against the frozen bar.
    flags = [r[2] for r in ordered]
    p_perm = permutation_p(flags, n_shuffles=20000)
    ci_lo, ci_hi = bootstrap_auc_ci(flags, n_boot=10000)
    act_ranks = [i + 1 for i, r in enumerate(ordered) if r[2]]
    print("\nSTATISTICAL SUPPORT (descriptive; does not gate)")
    print(f"  permutation test (20,000 label shuffles): p = {p_perm:.4f}")
    print(f"  bootstrap 95% CI over actives           : AUC {ci_lo:.3f} - {ci_hi:.3f}")
    print(f"  active rank positions (of {len(ordered)})           : "
          + ", ".join(str(r) for r in act_ranks))

    print("\n" + "=" * 70)
    passed = c[0] >= min_auc
    print(f"MODE C: AUC {c[0]:.3f} (need >= {min_auc})   EF@5% {c[1][1]:.2f}x   "
          f"BEDROC {c[2]:.3f}   [run-2 comparator AUC {RUN2['auc']:.3f}]")
    if passed:
        print("GATE: PASS — the geometric criterion separates real azoles from compounds")
        print("MEASURED not to inhibit. The separation published against presumed-inactive")
        print("decoys is not an artifact of decoy construction.")
        return 0
    print("GATE: FAIL — the criterion does not separate real azoles from measured")
    print("non-inhibitors. Per pre-registration this is the stronger negative: the ceiling")
    print("is not a property of how the decoys were built, and the shortlist claim must be")
    print("withdrawn or re-scoped. Do NOT re-filter the inactive set and re-run.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
