"""Ensemble within-azole gate — implements work/PREREGISTRATION_ensemble.md.

The within-azole question, one last permitted way. mode C (0.650) and mode F (0.590) both
failed to rank working azoles above measured-inactive azole analogues on the SINGLE rigid 5TZ1
conformer. Both FAILs pre-committed the same escape clause: the ceiling may be a property of the
one rigid pose, and the only permitted next attack is a real experimental conformer ENSEMBLE.

This grader is that attack, and it is disciplined exactly as the prereg demands: the criterion
is mode C, UNCHANGED (minimum ligand-nitrogen-to-heme-iron distance); the ensemble is EXHAUSTIVE
(all three C. albicans CYP51 crystal structures — 5TZ1, 5FSA, 5V5Z — so no conformer is chosen);
the aggregation is ensemble-MINIMUM, fixed in advance ("does there EXIST a conformer in which
this molecule coordinates the iron like a real azole?"). Zero free parameters beyond the frozen
protocol. The gate pool is the identical 7 held-out azoles + azole-bearing measured inactives on
which mode C scored 0.650 — a clean head-to-head.

Self-docking is leakage-free on the gate by construction (none of the 7 held-out actives was
co-crystallised in any conformer). The n=15 powered secondary re-introduces posaconazole (5FSA)
and itraconazole (5V5Z); their self-pair conformer is excluded there.

Usage:
    python scripts/validate_gate_ensemble.py    # reads the three screen dirs below
"""
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import protocol
from openafr.ensemble import ensemble_rows
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses
from openafr.scoring import bootstrap_auc_ci, metrics, permutation_p, report
from scripts.prep_receptor import sha256
from scripts.validate_gate2 import ACTIVES, check_prereg
from scripts.validate_gate_verified import azole_bearing, tip_statistics

PREREG = "work/PREREGISTRATION_ensemble.md"
PREREG_SHA = "8041b38ef37dde7a9a45243c9d6f5702ed568624826417f92c5a193549542c4e"
INACTIVES = "data/ligands/verified_inactives.smi"
TRAIN_ACTIVES = "data/ligands/actives.smi"

# The EXHAUSTIVE ensemble: id -> (screen dir, receptor anchor, pinned anchor sha256).
# Order matters for the marginal-conformer increment (5TZ1 first, the project's base receptor).
CONFORMERS = [
    ("5TZ1", "work/screen_verified", "work/receptor_A.pdb",
     "a1a8410868d1b9ee2835117cf6327ba4721ca3830a4d5689b7e9ea0ee1eee22f"),
    ("5FSA", "work/screen_5FSA", "work/receptor_5FSA_A.pdb",
     "9531a1531234072dae1f35df89a8411479393a91e170ff8d8e2b56af6bb777e3"),
    ("5V5Z", "work/screen_5V5Z", "work/receptor_5V5Z_A.pdb",
     "93a261010162782c625b36480568eb1b939b28414d9290817c2b4009e25103cd"),
]
# Self-docking guard: a co-crystal azole is not scored against the conformer solved with it.
SELF_PAIR = {"posaconazole": {"5FSA"}, "itraconazole": {"5V5Z"}}

MODE_C_SINGLE = 0.650        # look #6, single-conformer 5TZ1, this exact pool
MODE_F_SINGLE = 0.590        # look #7, channel engagement, this exact pool
NOVEL_SINGLE = 0.810         # look #6, single-conformer 5TZ1, non-azole inactives
GUARD_FLOOR = 0.78           # the novel-chemotype guard must stay ~>= this to call a real win


def best_nfe_by_molecule(screen, fe):
    """{name: min-over-poses N-Fe distance or None} for every pose file in a screen dir."""
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


def rows_to_key(rows, actives, exclude_names=None):
    """(name, distance) -> scoring rows (name, key, is_active); lower distance ranks first."""
    exclude_names = exclude_names or set()
    return [(n, d, n in actives) for n, d in rows if n not in exclude_names]


def auc_of(key_rows, gate):
    """AUC only, on (name, key, is_active) rows (unrankable placed last)."""
    ok = [r for r in key_rows if r[1] is not None]
    ordered = sorted(ok, key=lambda r: r[1]) + [r for r in key_rows if r[1] is None]
    return metrics([r[2] for r in ordered], gate["ef_fractions"], gate["bedroc_alpha"])[0]


def main():
    gate = protocol.load()["gate"]
    min_auc = gate["min_auc"]
    rep = lambda label, rows, gating=False: report(
        label, rows, gating, gate["ef_fractions"], gate["bedroc_alpha"])

    print("=" * 70)
    print("ENSEMBLE WITHIN-AZOLE GATE — mode C over the 3 C. albicans CYP51 conformers")
    print("=" * 70)
    check_prereg(PREREG, PREREG_SHA)
    protocol.verify()

    # Hash-gate every receptor anchor: what is graded must be what was pinned in the prereg.
    for cid, _screen, anchor, want in CONFORMERS:
        if not os.path.exists(anchor):
            print(f"!! {anchor} missing — run scripts/prep_receptor_ensemble.py")
            return 1
        got = sha256(anchor)
        ok = got == want
        print(f"receptor {cid:5}: {anchor}  sha256 {got[:16]}...  "
              + ("[pinned]" if ok else f"!! MOVED (expected {want[:16]}...)"))
        if not ok:
            print("   Any 'pass' from this run is not credible — the receptor changed.")
            return 1

    holdout = {l.split("\t")[1].strip() for l in open(ACTIVES) if "\t" in l}
    train = {l.split("\t")[1].strip() for l in open(TRAIN_ACTIVES) if "\t" in l}
    all_act = holdout | train
    azoles = azole_bearing(INACTIVES)

    # Per-conformer best-pose mode C for every molecule.
    per_conf = {}
    for cid, screen, anchor, _sha in CONFORMERS:
        if not os.path.isdir(screen):
            print(f"!! screen dir {screen} missing — dock this conformer first")
            return 1
        per_conf[cid] = best_nfe_by_molecule(screen, iron_position(anchor))
    names = sorted(set().union(*[set(d) for d in per_conf.values()]))
    for cid in per_conf:
        miss = set(names) - set(per_conf[cid])
        if miss:
            print(f"  !! {cid}: {len(miss)} molecules absent (docking failures): "
                  f"{sorted(miss)[:5]}{'...' if len(miss) > 5 else ''}")
    missing_act = holdout - set(names)
    if missing_act:
        print(f"\n!! MISSING HELD-OUT ACTIVES: {missing_act} — results not valid")

    def ens_rows(mode="min", conformer_ids=None, exclude=None):
        sub = {c: per_conf[c] for c in (conformer_ids or per_conf)}
        return ensemble_rows(names, sub, mode, exclude=exclude)

    # ---- PRIMARY (gates): ensemble-min, 7 held-out + azole-bearing inactives ----
    pool = holdout | azoles
    prim_rows = [(n, d) for n, d in ens_rows("min") if n in pool]
    c = rep(f"ENSEMBLE-MIN mode C  [7 held-out azoles vs azole-bearing inactives, n_act={len(holdout)}]",
            rows_to_key(prim_rows, holdout), gating=True)
    print(f"  single-conformer comparators on THIS pool: mode C {MODE_C_SINGLE:.3f} (look #6), "
          f"mode F {MODE_F_SINGLE:.3f} (look #7)")
    primary_auc = c[0]

    # ---- SECONDARY 1: per-conformer AUC on the gate pool (5TZ1 must ~= 0.650) ----
    print("\nPER-CONFORMER AUC on the gate pool (secondary; 5TZ1 alone must ~= 0.650 sanity)")
    for cid, _s, _a, _h in CONFORMERS:
        rows = [(n, d) for n, d in ens_rows("min", [cid]) if n in pool]
        a = auc_of(rows_to_key(rows, holdout), gate)
        flag = "  <- must ~= 0.650" if cid == "5TZ1" else ""
        print(f"  {cid:5}: AUC {a:.3f}{flag}")

    # ---- SECONDARY 2: ensemble-MEAN (diagnostic, never promoted) ----
    mean_rows = [(n, d) for n, d in ens_rows("mean") if n in pool]
    mean_auc = auc_of(rows_to_key(mean_rows, holdout), gate)
    print(f"\nENSEMBLE-MEAN mode C AUC (diagnostic, does NOT gate): {mean_auc:.3f}")

    # ---- SECONDARY 3: marginal-conformer increment ----
    print("\nMARGINAL-CONFORMER INCREMENT (ensemble-min, cumulative)")
    for k in range(1, len(CONFORMERS) + 1):
        ids = [c[0] for c in CONFORMERS[:k]]
        rows = [(n, d) for n, d in ens_rows("min", ids) if n in pool]
        print(f"  {{{', '.join(ids)}}}: AUC {auc_of(rows_to_key(rows, holdout), gate):.3f}")

    # ---- SECONDARY 4: n=15 powered form, self-pair excluded ----
    pool15 = all_act | azoles
    excl_rows = [(n, d) for n, d in ens_rows("min", exclude=SELF_PAIR) if n in pool15]
    a15 = auc_of(rows_to_key(excl_rows, all_act), gate)
    raw_rows = [(n, d) for n, d in ens_rows("min") if n in pool15]
    a15_raw = auc_of(rows_to_key(raw_rows, all_act), gate)
    print(f"\nn=15 POWERED (add 8 training azoles; posaconazole/itraconazole self-pair EXCLUDED)")
    print(f"  self-pair excluded (headline): AUC {a15:.3f}   |   raw (with self-pairs): {a15_raw:.3f}")

    # ---- SECONDARY 5: novel-chemotype GUARD (must not break the 0.810 that works) ----
    nonaz = {n for n in names if n not in azoles and n not in all_act}
    guard_pool7 = holdout | nonaz
    guard7_rows = [(n, d) for n, d in ens_rows("min") if n in guard_pool7]
    guard7 = auc_of(rows_to_key(guard7_rows, holdout), gate)
    guard_pool15 = all_act | nonaz
    guard15_rows = [(n, d) for n, d in ens_rows("min", exclude=SELF_PAIR) if n in guard_pool15]
    guard15 = auc_of(rows_to_key(guard15_rows, all_act), gate)
    print(f"\nNOVEL-CHEMOTYPE GUARD — ensemble-min vs NON-azole inactives (comparator {NOVEL_SINGLE:.3f})")
    print(f"  n_act=7 : AUC {guard7:.3f}   |   n_act=15: AUC {guard15:.3f}")
    print(f"  (an ensemble that lifts within-azole but drops this below ~{GUARD_FLOOR} is NOT a win)")

    # ---- statistical support on the primary pool ----
    ok = [r for r in rows_to_key(prim_rows, holdout) if r[1] is not None]
    ordered = sorted(ok, key=lambda r: r[1]) + [r for r in rows_to_key(prim_rows, holdout) if r[1] is None]
    flags = [r[2] for r in ordered]
    rank, above, pct = tip_statistics(ordered)
    print("\nSTATISTICAL SUPPORT (descriptive; the gate is the AUC point estimate)")
    print(f"  permutation test (20,000 shuffles): p = {permutation_p(flags, n_shuffles=20000):.4f}")
    ci_lo, ci_hi = bootstrap_auc_ci(flags, n_boot=10000)
    print(f"  bootstrap 95% CI over actives     : AUC {ci_lo:.3f} - {ci_hi:.3f}")
    if rank is not None:
        print(f"  best active                       : {ordered[rank-1][0]} at rank "
              f"{rank}/{len(ordered)} ({pct:.2f}th pct); {above} inactives above it")

    # ---- verdict, per the pre-committed interpretation ----
    print("\n" + "=" * 70)
    print(f"ENSEMBLE-MIN: AUC {primary_auc:.3f} (need >= {min_auc})   "
          f"[single-conformer C {MODE_C_SINGLE:.3f} / F {MODE_F_SINGLE:.3f}]   "
          f"guard(n=7) {guard7:.3f} [comparator {NOVEL_SINGLE:.3f}]")
    if primary_auc >= min_auc and guard7 >= GUARD_FLOOR:
        print("GATE: PASS — a real experimental conformer ensemble lifts within-azole ranking")
        print("above the single-rigid-conformer ceiling WITHOUT breaking novel-chemotype triage.")
        print("The mode C/F FAILs were induced-fit artefacts. Within-class triage re-opens (n=7,")
        print("wide CI); a separately pre-registered independent-active confirmation is required")
        print("before any within-class wet-lab claim.")
        return 0
    if primary_auc >= min_auc:
        print("GATE: within-azole cleared, but the NOVEL-CHEMOTYPE GUARD collapsed — the ensemble")
        print("just made coordination easier for everything. Per pre-registration this is NOT an")
        print("improvement; the single rigid conformer remains the product receptor.")
        return 1
    if primary_auc > MODE_C_SINGLE:
        print("GATE: PARTIAL — the ensemble helps but does not clear the bar. Three crystallographic")
        print("conformers are not enough flexibility. Per pre-registration the next step, if pursued,")
        print("is an MD-generated ensemble (separately registered) — NOT re-tuning this one.")
        return 1
    print("GATE: FAIL — the within-azole ceiling survives a complete, exhaustive experimental")
    print("conformer ensemble. The limitation is not captured by crystallographic receptor")
    print("flexibility. This closes the crystallographic-ensemble line; the defensible product")
    print("remains novel-chemotype triage (AUC 0.810).")
    return 1


if __name__ == "__main__":
    sys.exit(main())
