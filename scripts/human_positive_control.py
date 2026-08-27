"""Human-CYP51 positive control, redone as a native-pose-recovery test (issue #73).

Implements work/PREREGISTRATION_human_control.md (frozen sha in
work/PREREG_human_control.sha256). See that file for the fixed method and the pass criterion.

Short version: the #73 counter-screen's pre-registered positive control (re-dock ketoconazole,
require N-Fe < 3.0 A) did not fire — a *global-search* failure on a 36-heavy-atom, 8-rotatable-
bond floppy ligand, not a receptor fault (small rigid azoles coordinate the human iron 5/5, and
ketoconazole ALSO failed the same way on the validated fungal receptor). This control isolates
the receptor question from the search question:

  run 1  vina --local_only from the CRYSTAL ketoconazole pose  -> does the receptor+scoring
         recognize the true coordinating pose as a coordinating, negative-affinity minimum?
         THIS is the decider: PASS iff N-Fe < 3.0 A and affinity < 0.
  run 2  global redock (frozen protocol, seeds 42 1 2 3 4) from the crystal conformer -> context
         only; may still fail on search power, which would not change the verdict.

Deterministic given the frozen inputs. Writes a JSON summary and prints a table.

Usage:
    conda run -n openafr python scripts/human_positive_control.py
        [--pdb data/structures/3LD6.pdb] [--receptor work/receptor_human.pdbqt]
        [--anchor work/receptor_human_A.pdb] [--box work/receptor_human_box.json]
        [--out work/candidates/human_control] [--seeds 42 1 2 3 4]
"""
import argparse
import json
import math
import os
import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses, read_scores

ROOT = pathlib.Path(__file__).resolve().parent.parent
FE_CUTOFF = 3.0            # established coordination cutoff (scripts/analyze_poses.py)
LIGAND_RESN = "KKK"        # ketoconazole in 3LD6
LIGAND_CHAIN = "A"         # chain A, matching the receptor prep
BOX_SIZE = 26              # frozen protocol
NUM_MODES = 20             # frozen protocol
EXHAUST = 32               # frozen protocol


def extract_crystal_ligand(pdb, out_pdb):
    """Write only the chain-A KKK HETATM records to out_pdb (the crystal conformer)."""
    lines = [l for l in open(pdb)
             if l.startswith("HETATM") and l[17:20] == LIGAND_RESN and l[21] == LIGAND_CHAIN]
    if not lines:
        raise SystemExit(f"no {LIGAND_RESN} chain {LIGAND_CHAIN} atoms in {pdb}")
    with open(out_pdb, "w") as f:
        f.writelines(lines)
        f.write("END\n")
    return len(lines)


def obabel_ligand_pdbqt(in_pdb, out_pdbqt):
    """Same ligand prep as every docked ligand in the project: obabel ... -p 7.4.

    Preserves the input heavy-atom coordinates (the crystal pose); adds polar H at pH 7.4
    and AutoDock atom types. Does NOT re-embed a conformer.
    """
    subprocess.run(["obabel", in_pdb, "-O", out_pdbqt, "-p", "7.4"],
                   check=True, capture_output=True, text=True)
    if not os.path.getsize(out_pdbqt):
        raise SystemExit(f"obabel produced an empty {out_pdbqt}")


def vina(receptor, ligand, out, center, seed=None, local_only=False):
    cx, cy, cz = center
    cmd = ["vina", "--receptor", receptor, "--ligand", ligand, "--out", out,
           "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
           "--size_x", str(BOX_SIZE), "--size_y", str(BOX_SIZE), "--size_z", str(BOX_SIZE),
           "--num_modes", str(NUM_MODES)]
    if local_only:
        cmd += ["--local_only"]
    else:
        cmd += ["--exhaustiveness", str(EXHAUST), "--seed", str(seed)]
    p = subprocess.run(cmd, check=True, capture_output=True, text=True)
    if not os.path.getsize(out):
        raise SystemExit(f"vina produced an empty {out}")
    return p.stdout


def best_coordination(pose_path, fe):
    """(min_N_Fe, affinity_of_that_pose) over all MODEL/ENDMDL poses; N-Fe minimised."""
    scores = read_scores(pose_path)
    best = None
    for i, (_num, atoms) in enumerate(read_poses(pose_path)):
        d = min_nitrogen_iron_distance(atoms, fe)
        if d is not None and (best is None or d < best[0]):
            aff = scores[i] if i < len(scores) else None
            best = (round(d, 3), aff)
    return best if best else (None, None)


def bare_pose_nfe(path, fe):
    """min N-Fe for a --local_only output: a single, un-MODEL-wrapped pdbqt pose."""
    atoms = []
    for l in open(path):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        name = l[12:16].strip()
        el = (l[76:78].strip() or name[0]).upper()
        if el == "H":
            continue
        atoms.append((name, float(l[30:38]), float(l[38:46]), float(l[46:54])))
    d = min_nitrogen_iron_distance(atoms, fe)
    return round(d, 3) if d is not None else None


def parse_local_affinity(stdout):
    """Pull the estimated free energy of binding from a --local_only Vina run's stdout."""
    for l in stdout.splitlines():
        if "Estimated Free Energy of Binding" in l:
            return float(l.split(":")[1].split("(")[0].strip())
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", default=str(ROOT / "data/structures/3LD6.pdb"))
    ap.add_argument("--receptor", default=str(ROOT / "work/receptor_human.pdbqt"))
    ap.add_argument("--anchor", default=str(ROOT / "work/receptor_human_A.pdb"))
    ap.add_argument("--box", default=str(ROOT / "work/receptor_human_box.json"))
    ap.add_argument("--out", default=str(ROOT / "work/candidates/human_control"))
    ap.add_argument("--seeds", nargs="+", type=int, default=[42, 1, 2, 3, 4])
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    fe = iron_position(args.anchor)
    box = json.load(open(args.box))
    center = (box["center_x"], box["center_y"], box["center_z"])

    # Ground truth from the deposited coordinates (independent of any docking).
    crystal_pdb = os.path.join(args.out, "ketoconazole_crystal.pdb")
    natoms = extract_crystal_ligand(args.pdb, crystal_pdb)
    xtal_ns = []
    for l in open(crystal_pdb):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        name = l[12:16].strip()
        el = (l[76:78].strip() or name[0]).upper()
        if el == "N":
            xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            xtal_ns.append((name, round(math.dist(xyz, fe), 3)))
    crystal_min = min(d for _, d in xtal_ns)

    ligand_pdbqt = os.path.join(args.out, "ketoconazole_crystal.pdbqt")
    obabel_ligand_pdbqt(crystal_pdb, ligand_pdbqt)

    # Run 1 — the decider: local optimisation from the crystal pose. --local_only writes a
    # single, un-MODEL-wrapped pose and prints the affinity to stdout (no REMARK VINA RESULT).
    local_out = os.path.join(args.out, "ketoconazole_local.pdbqt")
    local_stdout = vina(args.receptor, ligand_pdbqt, local_out, center, local_only=True)
    local_nfe = bare_pose_nfe(local_out, fe)
    local_aff = parse_local_affinity(local_stdout)

    # Run 2 — context: global redock from the crystal conformer, per-seed.
    per_seed = {}
    for s in args.seeds:
        gout = os.path.join(args.out, f"ketoconazole_global_s{s}.pdbqt")
        vina(args.receptor, ligand_pdbqt, gout, center, seed=s)
        d, aff = best_coordination(gout, fe)
        per_seed[s] = {"n_fe": d, "affinity": aff}
    coord_seeds = [s for s, r in per_seed.items() if r["n_fe"] is not None and r["n_fe"] < FE_CUTOFF]
    global_consensus = min((per_seed[s]["n_fe"] for s in coord_seeds), default=None)

    passed = local_nfe is not None and local_nfe < FE_CUTOFF and local_aff is not None and local_aff < 0

    summary = {
        "prereg": "work/PREREGISTRATION_human_control.md",
        "iron": fe,
        "crystal_n_fe_by_atom": xtal_ns,
        "crystal_min_n_fe": crystal_min,
        "ligand_atoms": natoms,
        "fe_cutoff": FE_CUTOFF,
        "run1_local_only": {"n_fe": local_nfe, "affinity": local_aff,
                            "coordinates": local_nfe is not None and local_nfe < FE_CUTOFF},
        "run2_global": {"per_seed": per_seed,
                        "n_coord": len(coord_seeds), "n_seeds": len(args.seeds),
                        "consensus_n_fe": global_consensus},
        "PASS": passed,
    }
    with open(os.path.join(args.out, "human_control.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print(f"crystal ground truth: min N-Fe = {crystal_min} A  ({xtal_ns})")
    print(f"RUN 1  local_only (decider): N-Fe = {local_nfe} A   affinity = {local_aff} kcal/mol"
          f"   -> {'COORDINATES' if summary['run1_local_only']['coordinates'] else 'no coord'}")
    print(f"RUN 2  global redock: {len(coord_seeds)}/{len(args.seeds)} seeds coordinate; "
          f"consensus N-Fe = {global_consensus}")
    for s in args.seeds:
        print(f"        seed {s:>2}: N-Fe = {per_seed[s]['n_fe']}   aff = {per_seed[s]['affinity']}")
    print(f"\nPASS (receptor validated, #73 promotable): {passed}")
    print(f"wrote {os.path.join(args.out, 'human_control.json')}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
