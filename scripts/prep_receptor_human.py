"""3LD6.pdb -> work/receptor_human.pdbqt: the HUMAN CYP51 counter-screen receptor.

Pipeline stage: receptor prep — the off-target selectivity target (issue #73).

Why this exists
---------------
The whole method scores how closely a nitrogen reaches the FUNGAL heme iron. But humans
have their OWN CYP51 (lanosterol 14alpha-demethylase), and a molecule tuned to coordinate
the fungal heme iron can just as easily coordinate the human one — that is a core antifungal
toxicity mechanism (hepatotoxicity, drug interactions). A candidate is only interesting if it
reaches the fungal iron MORE readily than the human iron. This builds the human receptor so
the identical criterion can be run against it and a per-candidate selectivity margin computed.

Chain of custody mirrors scripts/prep_receptor.py exactly, so the two receptors are prepared
the same way and the only thing that differs between the fungal and human scores is the
protein:

    data/structures/3LD6.pdb          human CYP51 + ketoconazole (RCSB, chains A+B)
        |  keep chain A protein + chain-A HEM; drop chain B, waters, glycans (GLC),
        |  and the bound ketoconazole (KKK)
        v
    (chain-A atom set)  ==  work/receptor_human_A.pdb   <- tracked anchor the grader reads
        |                                                  the human heme-iron position from
        |  obabel -xr -p 7.4  (rigid receptor, polar H at pH 7.4, AutoDock atom types)
        v
    work/receptor_human.pdbqt         the docking input

The docking BOX for this receptor is centered on the crystallographic ketoconazole (KKK)
centroid — the human active site — exactly as the fungal box is centered on 5TZ1's bound VT1.
That center is derived here from the structure (never hand-copied) and written to
work/receptor_human_box.json, which scripts/screen_human.sh forwards to screen.sh as a
box-center override. Size, exhaustiveness, num_modes and seed stay from the frozen protocol.

Ground truth to verify against: 3LD6's own LINK record ties the ketoconazole N2 to the heme
iron at 2.41 A — so a correctly prepared receptor + a re-docked ketoconazole should reproduce
an iron-bound pose near that distance.
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_IN = os.path.join(ROOT, "data", "structures", "3LD6.pdb")
DEFAULT_OUT = os.path.join(ROOT, "work", "receptor_human.pdbqt")
ANCHOR = os.path.join(ROOT, "work", "receptor_human_A.pdb")
BOX_JSON = os.path.join(ROOT, "work", "receptor_human_box.json")

CHAIN = "A"            # crystal has two copies (A, B); use chain A, as the fungal prep does
HETERO_KEEP = {"HEM"}  # keep heme; drop waters (HOH), glycans (GLC), ketoconazole (KKK)
LIGAND = "KKK"         # ketoconazole — its chain-A centroid defines the docking box center


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def iron_xyz(path):
    """(x, y, z) of the heme iron, or None. Reads both the anchor PDB (HETATM FE) and
    obabel's pdbqt (ATOM FE)."""
    for l in open(path):
        if l[:6].strip() in ("ATOM", "HETATM") and l[12:16].strip() == "FE":
            return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return None


def ligand_centroid(pdb_path, chain=CHAIN, resname=LIGAND):
    """Centroid (x, y, z) of the named ligand's atoms in the given chain — the box center."""
    xs = ys = zs = n = 0
    for l in open(pdb_path):
        if (l[:6].strip() == "HETATM" and l[21] == chain
                and l[17:20].strip() == resname):
            xs += float(l[30:38]); ys += float(l[38:46]); zs += float(l[46:54]); n += 1
    if not n:
        sys.exit(f"ERROR: no {resname} atoms in chain {chain} — cannot place the docking box")
    return (xs / n, ys / n, zs / n, n)


def extract_chain(pdb_path):
    """Original 3LD6 lines for chain-A protein + chain-A heme, in file order."""
    kept = []
    for line in open(pdb_path):
        rec = line[:6].strip()
        if rec == "ATOM" and line[21] == CHAIN:
            kept.append(line)
        elif rec == "HETATM" and line[21] == CHAIN and line[17:20].strip() in HETERO_KEEP:
            kept.append(line)
    return kept


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdb_in", nargs="?", default=DEFAULT_IN)
    ap.add_argument("pdbqt_out", nargs="?", default=DEFAULT_OUT)
    ap.add_argument("--ph", default="7.4")
    args = ap.parse_args()

    if shutil.which("obabel") is None:
        sys.exit("ERROR: obabel not on PATH. Run inside the pinned env: conda activate openafr")
    if not os.path.exists(args.pdb_in):
        sys.exit(f"ERROR: input structure not found: {args.pdb_in}")

    print(f"input : {os.path.relpath(args.pdb_in, ROOT)}  (sha256 {sha256(args.pdb_in)[:12]})")

    kept = extract_chain(args.pdb_in)
    if not kept:
        sys.exit(f"ERROR: no chain-{CHAIN} atoms found in {args.pdb_in}")
    n_hem = sum(1 for l in kept if l[17:20].strip() == "HEM")
    print(f"  extracted chain {CHAIN}: {len(kept)} atoms "
          f"({len(kept) - n_hem} protein, {n_hem} heme)")

    # Write/verify the tracked anchor (the grader reads the human iron position from it).
    if os.path.exists(ANCHOR):
        ref = [l for l in open(ANCHOR) if l[:6].strip() in ("ATOM", "HETATM")]
        if len(ref) != len(kept):
            sys.exit(f"ERROR: extracted {len(kept)} atoms but anchor {ANCHOR} has {len(ref)}")
        print(f"  anchor check: matches existing {os.path.relpath(ANCHOR, ROOT)}")
    else:
        with open(ANCHOR, "w") as fh:
            fh.writelines(kept)
            fh.write("END\n")
        print(f"  wrote anchor: {os.path.relpath(ANCHOR, ROOT)} (sha256 {sha256(ANCHOR)[:12]})")

    # Box center = crystallographic ketoconazole centroid (the human active site).
    cx, cy, cz, nlig = ligand_centroid(args.pdb_in)
    with open(BOX_JSON, "w") as fh:
        json.dump({"center_x": round(cx, 3), "center_y": round(cy, 3),
                   "center_z": round(cz, 3), "source": f"{LIGAND} chain {CHAIN} centroid "
                   f"({nlig} atoms) in {os.path.basename(args.pdb_in)}"}, fh, indent=2)
    print(f"  box center: ({cx:.3f}, {cy:.3f}, {cz:.3f}) from {LIGAND} "
          f"-> {os.path.relpath(BOX_JSON, ROOT)}")

    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as tmp:
        tmp.writelines(kept)
        tmp.write("END\n")
        tmp_path = tmp.name
    try:
        os.makedirs(os.path.dirname(os.path.abspath(args.pdbqt_out)), exist_ok=True)
        proc = subprocess.run(["obabel", tmp_path, "-O", args.pdbqt_out, "-xr", "-p", args.ph],
                              capture_output=True, text=True)
    finally:
        os.unlink(tmp_path)

    if proc.returncode != 0 or not os.path.exists(args.pdbqt_out) or os.path.getsize(args.pdbqt_out) == 0:
        sys.exit(f"ERROR: obabel conversion failed\n{proc.stderr}")

    # Name the REMARK after the real structure so the receptor hashes reproducibly.
    lines = open(args.pdbqt_out).readlines()
    for i, line in enumerate(lines):
        if line.startswith("REMARK") and "Name =" in line:
            lines[i] = f"REMARK  Name = {os.path.basename(args.pdb_in)}\n"
            break
    with open(args.pdbqt_out, "w") as fh:
        fh.writelines(lines)

    fe = iron_xyz(args.pdbqt_out)
    if fe is None:
        sys.exit("ERROR: no heme FE atom in the output receptor — conversion dropped the iron")
    fe_ref = iron_xyz(ANCHOR)
    if fe_ref is not None and max(abs(a - b) for a, b in zip(fe, fe_ref)) > 1e-3:
        sys.exit(f"ERROR: heme iron moved during conversion: {fe} vs anchor {fe_ref}")

    print(f"output: {os.path.relpath(args.pdbqt_out, ROOT)}  (sha256 {sha256(args.pdbqt_out)[:12]})")
    print(f"  heme Fe at ({fe[0]:.3f}, {fe[1]:.3f}, {fe[2]:.3f})  [matches anchor]")
    print("done. Verify by re-docking ketoconazole: expect an iron-bound pose near the "
          "crystallographic 2.41 A N-Fe (3LD6 LINK record).")


if __name__ == "__main__":
    main()
