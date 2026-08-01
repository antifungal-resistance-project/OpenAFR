"""Parsers for AutoDock Vina output (.pdbqt) and receptor structures (.pdb).

One canonical copy of the three parsers that every gate/analysis script needs.
All coordinate readers skip hydrogens: docking is done on the polar/heavy-atom
skeleton and the iron-coordination geometry is measured between heavy atoms only.
"""
import math


def read_poses(path):
    """Yield (model_number, [(atom_name, x, y, z), ...]) for each docked pose.

    Heavy atoms only. Element is taken from PDB columns 77-78 when present, else
    inferred from the first character of the atom name.
    """
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
    """Return the Vina affinity (kcal/mol) of every pose, in file order (best first)."""
    return [float(l.split()[3]) for l in open(path) if "REMARK VINA RESULT" in l]


def iron_position(receptor_pdb):
    """Return (x, y, z) of the heme iron. Raise SystemExit if the receptor has none."""
    for l in open(receptor_pdb):
        if l.startswith("HETATM") and l[76:78].strip().upper() == "FE":
            return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    raise SystemExit(f"no heme iron found in {receptor_pdb} — wrong structure?")


def min_nitrogen_iron_distance(atoms, fe):
    """Closest approach of any nitrogen in one pose to the heme iron, or None.

    None means the pose has no nitrogen at all and so cannot be a real azole pose.
    """
    ns = [a for a in atoms if a[0].startswith("N")]
    if not ns:
        return None
    return min(math.dist(a[1:], fe) for a in ns)
