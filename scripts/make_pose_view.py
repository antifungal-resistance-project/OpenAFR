"""Export a molecule's best iron-coordinating pose for 3D viewing.

Given a docked screen (one dir, or a consensus parent with per-seed s<seed>/ subdirs) and
a molecule name, this finds the single pose with the smallest nitrogen-to-heme-iron distance,
writes it as a standalone PDB alongside the receptor, and emits a PyMOL script that renders a
labelled scene (pocket, heme, iron, the ligand, and the measured N-Fe distance).

Uses only the pinned env (python + obabel); it does NOT import or require PyMOL. The emitted
.pml is plain text you can open in any PyMOL/ChimeraX, or render headless with the optional
`pymolview` env (see scripts/render_views.sh).

    python scripts/make_pose_view.py efinaconazole \
        --screen work/screen_consensus_wt --receptor work/receptor_A.pdb --out work/views/efinaconazole_wt
"""
import argparse
import glob
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses


def pose_files(screen, name):
    """All docked .pdbqt files for `name`: consensus s<seed>/ subdirs, or a flat dir."""
    seeded = [d for d in glob.glob(os.path.join(screen, "s*"))
              if os.path.isdir(d) and os.path.basename(d)[1:].isdigit()]
    dirs = seeded or [screen]
    return [os.path.join(d, f"{name}.pdbqt") for d in dirs
            if os.path.exists(os.path.join(d, f"{name}.pdbqt"))]


def best_pose(files, fe):
    """(file, model_number, distance) of the single closest-approach pose across files."""
    best = None
    for f in files:
        for num, atoms in read_poses(f):
            d = min_nitrogen_iron_distance(atoms, fe)
            if d is not None and (best is None or d < best[2]):
                best = (f, num, d)
    return best


def extract_model(pdbqt, model_num, out_pdbqt):
    """Write just MODEL `model_num` from a multi-model docked pdbqt to its own file."""
    keep, on = [], False
    for line in open(pdbqt):
        if line.startswith("MODEL"):
            on = int(line.split()[1]) == model_num
        if on:
            keep.append(line)
        if line.startswith("ENDMDL") and on:
            break
    with open(out_pdbqt, "w") as f:
        f.writelines(keep)


PML = """# Auto-generated scene — {name} best iron-coordinating pose (N-Fe = {dist:.3f} A)
load {receptor}, receptor
load {ligand}, ligand
hide everything
bg_color white
# protein pocket
show cartoon, receptor
color grey80, receptor and polymer
set cartoon_transparency, 0.5, receptor
# heme + iron
show sticks, receptor and resn HEM
color salmon, receptor and resn HEM
show spheres, receptor and elem Fe
color orange, receptor and elem Fe
set sphere_scale, 0.5, elem Fe
# the docked molecule
show sticks, ligand
util.cbay ligand
# the coordinating contact: closest ligand N to the iron
select ironN, (ligand and elem N) within 4.0 of (receptor and elem Fe)
distance nfe, ironN, (receptor and elem Fe), 4.0
color black, nfe
# frame it on the iron
orient (receptor and elem Fe) expand 8
zoom (receptor and elem Fe), 10
set ray_shadows, 0
ray 1400, 1050
png {png}, dpi=150
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("name")
    ap.add_argument("--screen", required=True)
    ap.add_argument("--receptor", default="work/receptor_A.pdb")
    ap.add_argument("--out", required=True, help="output path prefix")
    args = ap.parse_args()

    fe = iron_position(args.receptor)
    files = pose_files(args.screen, args.name)
    if not files:
        sys.exit(f"no docked poses for {args.name} under {args.screen}")
    bp = best_pose(files, fe)
    if bp is None:
        sys.exit(f"{args.name} has no nitrogen-bearing pose to view")
    src, model, dist = bp

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    lig_pdb = f"{args.out}.pdb"
    with tempfile.NamedTemporaryFile("w", suffix=".pdbqt", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        extract_model(src, model, tmp_path)
        subprocess.run(["obabel", tmp_path, "-O", lig_pdb], capture_output=True)
    finally:
        os.unlink(tmp_path)

    pml_path = f"{args.out}.pml"
    with open(pml_path, "w") as f:
        f.write(PML.format(name=args.name, dist=dist,
                           receptor=os.path.abspath(args.receptor),
                           ligand=os.path.abspath(lig_pdb),
                           png=os.path.abspath(f"{args.out}.png")))
    print(f"{args.name}: best pose N-Fe {dist:.3f} A  (from {os.path.relpath(src)}, model {model})")
    print(f"  ligand PDB : {lig_pdb}")
    print(f"  PyMOL scene: {pml_path}  ->  {args.out}.png")


if __name__ == "__main__":
    main()
