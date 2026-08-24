"""Channel-engagement gate — implements work/PREREGISTRATION_channel_engagement.md.

The within-azole question look #6 left open. Look #6 showed the iron-distance criterion
(mode C) scores AUC 0.650 against azole-bearing measured inactives — below the bar. It cannot
tell a working azole from a failed azole analogue, because every azole coordinates the iron.
This is one orthogonal, mechanistically-fixed feature at that gap, under the same one-look
discipline the axial pre-registration used to reuse already-docked poses.

mode F — channel engagement. Azole SAR turns on the N1 tail threading the CYP51 substrate-
access channel, a signal the N-Fe distance is blind to. Per molecule:
  * channel-shell = protein residues within 4.5 A of the co-crystal azole VT1 (label-blind:
    the residues the real bound drug touches, defined from the crystal, using NO test labels);
  * measured in the COORDINATING pose (min N-Fe, the same pose mode C selects) to remove any
    over-poses aggregation degree of freedom;
  * feature = count of distinct channel-shell residues the ligand contacts (<= 4.5 A), ranked
    descending. Zero free parameters beyond the repo-wide 4.5 A cutoff.

Gate: AUC >= protocol min_auc on the matched pool (7 held-out azole actives + the azole-bearing
verified inactives) — the identical pool where mode C scored 0.650. AUC is threshold-free, so
there is no cutoff to tune. Everything else (n=15 powered form, r(F,C), r(F,size), the size-
residualized AUC, EF/BEDROC, permutation p, bootstrap CI, the tip) is reported, never gates.

Usage:
    python scripts/validate_gate_channel.py [SCREEN_DIR] [RECEPTOR]
"""
import math
import os
import pathlib
import sys

from rdkit import Chem

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import mapping, protocol
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses
from openafr.scoring import bootstrap_auc_ci, metrics, permutation_p, report
from scripts.validate_gate2 import ACTIVES, check_prereg
from scripts.validate_gate_verified import RUN2, azole_bearing, tip_statistics

PREREG = "work/PREREGISTRATION_channel_engagement.md"
PREREG_SHA = "62592a101ad893f193a0239d47dc04348967bfcf4138184f5db299d6ce516548"
DEFAULT_SCREEN = "work/screen_verified"
INACTIVES = "data/ligands/verified_inactives.smi"
TRAIN_ACTIVES = "data/ligands/actives.smi"     # the 8 training azoles (n=15 powered secondary)
REF_LIGAND = "work/ref_VT1.pdb"
CUTOFF = mapping.CONTACT_CUTOFF                 # 4.5 A, fixed repo-wide
MODE_C_WITHIN_AZOLE = 0.650                     # look #6 comparator on this exact pool


def channel_shell(receptor, ref_ligand=REF_LIGAND, cutoff=CUTOFF):
    """Protein residue positions within `cutoff` of the co-crystal ligand — label-blind.

    Defined from the crystal structure alone (the residues the real bound drug touches); no
    active/inactive label is consulted. read_residues reads protein ATOM records only, so the
    heme HETATM group is excluded automatically — this is the protein channel, not the
    coordination shell.
    """
    residues = mapping.read_residues(receptor)
    lig = mapping.read_ligand_atoms(ref_ligand)
    shell = {}
    for pos, (_resname, atoms) in residues.items():
        d = mapping._min_dist(atoms, lig)
        if d is not None and d <= cutoff:
            shell[pos] = atoms
    return shell


def coordinating_pose(path, fe):
    """The (atoms of the) pose with the smallest N-Fe distance — the same pose mode C uses.

    Returns None when the molecule has no nitrogen-bearing pose (no coordinating pose to
    select); such a molecule is unrankable and placed last, exactly as under mode C.
    """
    best_d, best_atoms = None, None
    for _num, atoms in read_poses(path):
        d = min_nitrogen_iron_distance(atoms, fe)
        if d is not None and (best_d is None or d < best_d):
            best_d, best_atoms = d, atoms
    return best_atoms, best_d


def engagement(pose_atoms, shell, cutoff=CUTOFF):
    """Distinct channel-shell residues the pose contacts (any ligand heavy atom <= cutoff)."""
    n = 0
    for _pos, res_atoms in shell.items():
        if any(math.dist((lx, ly, lz), rxyz) <= cutoff
               for (_ln, lx, ly, lz) in pose_atoms for (_rn, rxyz) in res_atoms):
            n += 1
    return n


def heavy_atom_count(pose_atoms):
    """Ligand heavy-atom count (read_poses already drops hydrogens)."""
    return len(pose_atoms)


def _pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return num / (dx * dy) if dx and dy else float("nan")


def _residualize(feat, size):
    """Residuals of `feat` after an OLS fit on `size` (removes the linear size component)."""
    n = len(feat)
    ms, mf = sum(size) / n, sum(feat) / n
    sxx = sum((s - ms) ** 2 for s in size)
    if not sxx:
        return list(feat)
    b = sum((s - ms) * (f - mf) for s, f in zip(size, feat)) / sxx
    a = mf - b * ms
    return [f - (a + b * s) for f, s in zip(feat, size)]


def collect(screen, actives, fe, shell):
    """Rows (name, engagement, is_active, heavy_atoms, nfe) for every pose file in `screen`.

    engagement is None (unrankable, placed last) when there is no coordinating pose.
    """
    rows = []
    for f in sorted(os.listdir(screen)):
        if not f.endswith(".pdbqt"):
            continue
        name = f[:-6]
        pose, nfe = coordinating_pose(os.path.join(screen, f), fe)
        eng = engagement(pose, shell) if pose is not None else None
        ha = heavy_atom_count(pose) if pose is not None else None
        rows.append((name, eng, name in actives, ha, nfe))
    return rows


def _rank_key(rows):
    """Map rows to scoring.report's (name, key, is_active): higher engagement ranks first,
    so the sort key is NEGATIVE engagement; unrankable (None) stays None -> placed last."""
    return [(r[0], (-r[1] if r[1] is not None else None), r[2]) for r in rows]


def main():
    screen = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SCREEN
    receptor = sys.argv[2] if len(sys.argv) > 2 else "work/receptor_A.pdb"
    fe = iron_position(receptor)
    holdout = {l.split("\t")[1].strip() for l in open(ACTIVES) if "\t" in l}
    train = {l.split("\t")[1].strip() for l in open(TRAIN_ACTIVES) if "\t" in l} \
        if os.path.exists(TRAIN_ACTIVES) else set()

    gate = protocol.load()["gate"]
    min_auc = gate["min_auc"]
    rep = lambda label, rows, gating=False: report(
        label, rows, gating, gate["ef_fractions"], gate["bedroc_alpha"])

    print("=" * 70)
    print("CHANNEL-ENGAGEMENT GATE — THE WITHIN-AZOLE QUESTION (mode F, orthogonal feature)")
    print("=" * 70)
    print(f"receptor         : {receptor}")
    print(f"poses            : {screen}")
    check_prereg(PREREG, PREREG_SHA)
    protocol.verify()

    shell = channel_shell(receptor)
    print(f"channel-shell    : {len(shell)} protein residues within {CUTOFF:.1f} A of VT1 "
          f"(label-blind, from the crystal)")

    if not os.path.exists(INACTIVES):
        print(f"!! {INACTIVES} missing — cannot select the azole-bearing subgroup")
        return 1
    azoles = azole_bearing(INACTIVES)                 # azole-bearing verified inactives

    all_actives = holdout | train
    rows = collect(screen, all_actives, fe, shell)
    have = {r[0] for r in rows}
    missing = holdout - have
    if missing:
        print(f"\n!! MISSING HELD-OUT ACTIVES: {missing} — results not valid")

    # Primary pool: 7 held-out azole actives + azole-bearing inactives (matched to look #6).
    prim = [r for r in rows if (r[0] in holdout) or (r[0] in azoles)]
    c = rep(f"MODE F — channel engagement  [vs azole-bearing inactives, n_act={len(holdout)}]",
            _rank_key(prim), gating=True)
    print(f"  mode-C comparator on THIS pool (look #6): AUC {MODE_C_WITHIN_AZOLE:.3f}")

    # ---- secondaries (reported, never gate) ----
    # n=15 powered form: add the 8 training azoles (held-out w.r.t. THIS feature).
    if train:
        pool15 = [r for r in rows if (r[0] in all_actives) or (r[0] in azoles)]
        n_act15 = sum(1 for r in pool15 if r[2])
        rep(f"MODE F — powered form (n_act={n_act15}: +8 training azoles)", _rank_key(pool15))

    # Orthogonality + size diagnostics on the primary pool (rankable rows only).
    ok = [r for r in prim if r[1] is not None and r[4] is not None]
    eng = [r[1] for r in ok]
    ha = [r[3] for r in ok]
    nfe = [r[4] for r in ok]
    print("\nORTHOGONALITY / SIZE DIAGNOSTICS (primary pool, descriptive)")
    print(f"  r(mode F, mode C = N-Fe)  : {_pearson(eng, nfe):+.3f}   "
          f"(near 0 = genuinely new signal)")
    print(f"  r(mode F, heavy-atom count): {_pearson(eng, ha):+.3f}")

    # Size-adjusted AUC: rank by engagement residualized on heavy-atom count.
    resid = _residualize(eng, ha)
    resid_rows = [(ok[i][0], -resid[i], ok[i][2]) for i in range(len(ok))]
    ordered_r = sorted(resid_rows, key=lambda r: r[1])
    auc_resid = metrics([r[2] for r in ordered_r], gate["ef_fractions"], gate["bedroc_alpha"])[0]
    print(f"  size-residualized AUC     : {auc_resid:.3f}   "
          f"(collapse here = mode F was mostly re-measuring size)")

    # Statistical support + the tip, exactly as look #6.
    ordered = sorted(ok, key=lambda r: -r[1]) + [r for r in prim if r[1] is None]
    flags = [r[2] for r in ordered]
    rank, above, pct = tip_statistics([(r[0], r[1], r[2]) for r in ordered])
    print("\nSTATISTICAL SUPPORT (descriptive; does not gate)")
    print(f"  permutation test (20,000 shuffles): p = {permutation_p(flags, n_shuffles=20000):.4f}")
    ci_lo, ci_hi = bootstrap_auc_ci(flags, n_boot=10000)
    print(f"  bootstrap 95% CI over actives     : AUC {ci_lo:.3f} - {ci_hi:.3f}")
    if rank is not None:
        print(f"  best active                       : {ordered[rank-1][0]} at rank "
              f"{rank}/{len(ordered)} ({pct:.2f}th pct); {above} inactives above it")

    print("\n" + "=" * 70)
    passed = c[0] >= min_auc
    print(f"MODE F: AUC {c[0]:.3f} (need >= {min_auc})   size-residualized AUC {auc_resid:.3f}   "
          f"[mode-C on this pool {MODE_C_WITHIN_AZOLE:.3f}]")
    if passed:
        print("GATE: PASS — channel engagement separates working azoles from failed azole")
        print("analogues where iron-distance could not. Read the size-residualized AUC before")
        print("claiming a size-INDEPENDENT fit signal (see the pre-committed interpretation).")
        return 0
    print("GATE: FAIL — an orthogonal, mechanistically-motivated tail-fit feature also cannot")
    print("rank within the azole class on a single rigid conformer. Per pre-registration this is")
    print("the decisive evidence for Option B (re-scope to novel-chemotype triage, AUC 0.810).")
    print("Do NOT re-define the feature and re-run; the next attack needs ensemble/flexible poses.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
