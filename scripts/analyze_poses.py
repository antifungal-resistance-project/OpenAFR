"""Score docked poses with a metal-coordination constraint.

Pipeline stage: validate/ (seed)

Why this exists
---------------
AutoDock Vina's scoring function ignores partial charges, so it cannot model the
azole-nitrogen -> heme-iron coordination bond that defines how these drugs work.
Empirically (see work/RESULTS_redock_VT1.md) Vina DOES generate the correct pose
but misranks it, behind an iron-free pose, by less than 0.1 kcal/mol.

So we never rank on docking score alone. A pose only counts if it makes the iron
bond. Chemistry first, score second.

    all poses
        |
        +--> N-Fe distance < CUTOFF ? --no--> rejected (cannot be a real azole pose)
        |                              |
        |                             yes
        |                              v
        +------------------> best score among iron-bound poses = the answer
"""
import math
import os
import sys

FE_CUTOFF = 3.0  # angstroms; a real azole-iron coordination bond is 2.0-2.5 A


def read_poses(path):
    """Yield (model_number, [(atom_name, x, y, z), ...]) for each pose, heavy atoms only."""
    cur, num = [], None
    for line in open(path):
        if line.startswith("MODEL"):
            num, cur = int(line.split()[1]), []
        elif line.startswith("ENDMDL"):
            yield num, cur
        elif line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            element = (line[76:78].strip() or name[0]).upper()
            if element == "H":
                continue
            cur.append((name, float(line[30:38]), float(line[38:46]), float(line[46:54])))


def read_scores(path):
    return [float(l.split()[3]) for l in open(path) if "REMARK VINA RESULT" in l]


def iron_position(receptor_pdb):
    for l in open(receptor_pdb):
        if l.startswith("HETATM") and l[76:78].strip().upper() == "FE":
            return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    raise SystemExit(f"No heme iron found in {receptor_pdb} — wrong structure?")


def analyze(pdbqt, fe):
    """Return (best_score, best_constrained_score, best_fe_distance, n_poses, n_passing)."""
    scores = read_scores(pdbqt)
    poses = list(read_poses(pdbqt))
    if not poses:
        return None
    best_fe, constrained, all_fe = None, [], []
    for (num, atoms), score in zip(poses, scores):
        nitrogens = [a for a in atoms if a[0].startswith("N")]
        if not nitrogens:
            continue
        d = min(math.dist(a[1:], fe) for a in nitrogens)
        all_fe.append(d)
        if d < FE_CUTOFF:
            constrained.append((score, num, d))
    if all_fe:
        best_fe = min(all_fe)
    return {
        "best_score": scores[0] if scores else None,
        "n_poses": len(poses),
        "n_passing": len(constrained),
        "best_fe": best_fe,
        "constrained": min(constrained) if constrained else None,
    }


def main():
    docked_dir = sys.argv[1] if len(sys.argv) > 1 else "work/docked"
    receptor = sys.argv[2] if len(sys.argv) > 2 else "work/receptor_A.pdb"
    fe = iron_position(receptor)

    files = sorted(f for f in os.listdir(docked_dir) if f.endswith(".pdbqt"))
    print(f"Heme iron at ({fe[0]:.2f}, {fe[1]:.2f}, {fe[2]:.2f})")
    print(f"Metal-coordination cutoff: N-Fe < {FE_CUTOFF} A\n")
    print(f"{'ligand':16}{'top score':>10}{'poses':>7}{'iron-bound':>12}"
          f"{'closest N-Fe':>14}{'constrained score':>19}")
    print("-" * 78)

    passed, failed = [], []
    for f in files:
        name = f[:-6]
        r = analyze(os.path.join(docked_dir, f), fe)
        if r is None:
            print(f"{name:16}{'NO POSES — docking failed':>50}")
            failed.append(name)
            continue
        cs = f"{r['constrained'][0]:.2f}" if r["constrained"] else "none"
        print(f"{name:16}{r['best_score']:>10.2f}{r['n_poses']:>7}{r['n_passing']:>12}"
              f"{r['best_fe']:>14.2f}{cs:>19}")
        (passed if r["n_passing"] else failed).append(name)

    print()
    print(f"{len(passed)}/{len(files)} ligands produced at least one iron-bound pose")
    if failed:
        print(f"no iron-bound pose: {', '.join(failed)}")


if __name__ == "__main__":
    main()
