"""Ensemble confirmation gate — implements work/PREREGISTRATION_ensemble_confirm.md.

Look #8 (work/RESULTS_ensemble.md) passed within-azole at ensemble-min mode C AUC 0.750, but
on only the 7 held-out azoles (bootstrap CI 0.683-0.825 straddles the bar). The ensemble
pre-registration's PASS branch requires the payment before any within-class claim: a
"separately pre-registered confirmation on an INDEPENDENT active set."

This is that confirmation. Everything is held fixed from look #8 — the exhaustive 3-conformer
C. albicans CYP51 ensemble (hash-pinned), the mode C criterion, ensemble-MINIMUM aggregation,
the frozen protocol, and the SAME frozen azole-bearing verified-inactive pool. The ONLY thing
that changes is the actives: an independent set of measured-active azoles built from the ChEMBL
C. albicans cache by a frozen rule (scripts/make_verified_actives.py), novelty-filtered < 0.70
against every molecule look #8 used and against the 5TZ1 co-crystal ligand. If ensemble-min
still separates these DIFFERENT, measured-active azoles from measured-inactive azole analogues
at AUC >= 0.70, the look #8 result was not a fluke of those particular 7 molecules.

Usage:
    python scripts/validate_gate_ensemble_confirm.py
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
from scripts.validate_gate2 import check_prereg
from scripts.validate_gate_ensemble import (
    CONFORMERS, MODE_C_SINGLE, MODE_F_SINGLE, NOVEL_SINGLE, GUARD_FLOOR,
    best_nfe_by_molecule, rows_to_key, auc_of)
from scripts.validate_gate_verified import azole_bearing, tip_statistics

PREREG = "work/PREREGISTRATION_ensemble_confirm.md"
PREREG_SHA = "7754050e94523a9c2f5d2f9626708cd5bd02f12122f32e5d549049eccfb05c3b"
INACTIVES = "data/ligands/verified_inactives.smi"
ACTIVES = "data/ligands/verified_actives_sample.smi"   # the frozen N=200 seed-42 test set
ENSEMBLE_MIN_LOOK8 = 0.750                          # the result this run confirms or overturns


def main():
    gate = protocol.load()["gate"]
    min_auc = gate["min_auc"]
    rep = lambda label, rows, gating=False: report(
        label, rows, gating, gate["ef_fractions"], gate["bedroc_alpha"])

    print("=" * 70)
    print("ENSEMBLE CONFIRMATION GATE — independent measured-active azole set")
    print("=" * 70)
    check_prereg(PREREG, PREREG_SHA)
    protocol.verify()
    for cid, _screen, anchor, want in CONFORMERS:
        got = sha256(anchor)
        ok = got == want
        print(f"receptor {cid:5}: sha256 {got[:16]}...  " + ("[pinned]" if ok else "!! MOVED"))
        if not ok:
            print("   Any 'pass' from this run is not credible — the receptor changed.")
            return 1

    actives = {l.split("\t")[1].strip() for l in open(ACTIVES) if "\t" in l}
    azoles = azole_bearing(INACTIVES)
    print(f"independent actives: {len(actives)}  (from {ACTIVES})")
    print(f"azole-bearing inactives (frozen, same as look #8): {len(azoles)}")

    per_conf = {}
    for cid, screen, anchor, _sha in CONFORMERS:
        if not os.path.isdir(screen):
            print(f"!! screen dir {screen} missing — dock this conformer first")
            return 1
        per_conf[cid] = best_nfe_by_molecule(screen, iron_position(anchor))
    names = sorted(set().union(*[set(d) for d in per_conf.values()]))
    missing_act = actives - set(names)
    if missing_act:
        print(f"\n!! {len(missing_act)} independent actives have no pose (docking failures): "
              f"{sorted(missing_act)[:5]}{'...' if len(missing_act) > 5 else ''}")
        print("   They are ranked LAST (never dropped), which is conservative for the gate.")

    def ens_rows(mode="min", conformer_ids=None):
        sub = {c: per_conf[c] for c in (conformer_ids or per_conf)}
        return ensemble_rows(names, sub, mode)

    pool = actives | azoles
    prim = [(n, d) for n, d in ens_rows("min") if n in pool]
    c = rep(f"ENSEMBLE-MIN mode C  [independent actives vs azole-bearing inactives, "
            f"n_act={len(actives & set(names))}]", rows_to_key(prim, actives), gating=True)
    print(f"  comparators: look #8 ensemble-min {ENSEMBLE_MIN_LOOK8:.3f} (7 held-out azoles); "
          f"single-conformer mode C {MODE_C_SINGLE:.3f}, mode F {MODE_F_SINGLE:.3f}")
    primary_auc = c[0]

    print("\nPER-CONFORMER AUC on this pool (secondary)")
    for cid, _s, _a, _h in CONFORMERS:
        rows = [(n, d) for n, d in ens_rows("min", [cid]) if n in pool]
        print(f"  {cid:5}: AUC {auc_of(rows_to_key(rows, actives), gate):.3f}")

    mean_rows = [(n, d) for n, d in ens_rows("mean") if n in pool]
    print(f"\nENSEMBLE-MEAN mode C AUC (diagnostic, does NOT gate): "
          f"{auc_of(rows_to_key(mean_rows, actives), gate):.3f}")

    print("\nMARGINAL-CONFORMER INCREMENT (ensemble-min, cumulative)")
    for k in range(1, len(CONFORMERS) + 1):
        ids = [cc[0] for cc in CONFORMERS[:k]]
        rows = [(n, d) for n, d in ens_rows("min", ids) if n in pool]
        print(f"  {{{', '.join(ids)}}}: AUC {auc_of(rows_to_key(rows, actives), gate):.3f}")

    nonaz = {n for n in names if n not in azoles and n not in actives}
    guard_rows = [(n, d) for n, d in ens_rows("min") if n in (actives | nonaz)]
    guard = auc_of(rows_to_key(guard_rows, actives), gate)
    print(f"\nNOVEL-CHEMOTYPE GUARD — independent actives vs NON-azole inactives "
          f"(comparator {NOVEL_SINGLE:.3f}): AUC {guard:.3f}")

    ok = [r for r in rows_to_key(prim, actives) if r[1] is not None]
    ordered = sorted(ok, key=lambda r: r[1]) + [r for r in rows_to_key(prim, actives) if r[1] is None]
    flags = [r[2] for r in ordered]
    rank, above, pct = tip_statistics(ordered)
    print("\nSTATISTICAL SUPPORT (descriptive; the gate is the AUC point estimate)")
    print(f"  permutation test (20,000 shuffles): p = {permutation_p(flags, n_shuffles=20000):.4f}")
    ci_lo, ci_hi = bootstrap_auc_ci(flags, n_boot=10000)
    print(f"  bootstrap 95% CI over actives     : AUC {ci_lo:.3f} - {ci_hi:.3f}")
    if rank is not None:
        print(f"  best active                       : {ordered[rank-1][0]} at rank "
              f"{rank}/{len(ordered)} ({pct:.2f}th pct); {above} inactives above it")

    print("\n" + "=" * 70)
    print(f"ENSEMBLE-MIN (independent actives): AUC {primary_auc:.3f} (need >= {min_auc})   "
          f"[look #8 {ENSEMBLE_MIN_LOOK8:.3f}]   guard {guard:.3f}")
    if primary_auc >= min_auc and guard >= GUARD_FLOOR:
        print("GATE: PASS — the within-azole lift REPLICATES on an independent measured-active")
        print("azole set. The look #8 result was not a fluke of the 7 held-out azoles. Within-class")
        print("triage is confirmed at the ensemble level (still a triage AUC, not top-of-list")
        print("enrichment — read the EF/BEDROC and the tip before any wet-lab prioritisation).")
        return 0
    if primary_auc >= min_auc:
        print("GATE: within-azole replicates but the NOVEL-CHEMOTYPE GUARD dropped — read both.")
        return 1
    print("GATE: FAIL — the within-azole lift does NOT replicate on an independent active set.")
    print("The look #8 AUC 0.750 was consistent with a 7-molecule sampling fluctuation (its CI")
    print("straddled the bar). No within-class claim is earned; novel-chemotype triage (0.810)")
    print("remains the defensible product.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
