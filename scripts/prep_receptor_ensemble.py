"""Extra C. albicans CYP51 conformers -> aligned docking receptors, built by script.

Pipeline stage: receptor prep for the ENSEMBLE within-azole look
(work/PREREGISTRATION_ensemble.md). The single rigid 5TZ1 conformer could not rank within
the azole class (mode C 0.650, mode F 0.590). The one permitted next attack is a *real
experimental conformer ensemble* — every C. albicans CYP51 crystal structure in the PDB,
which is exactly three: 5TZ1 (already the project receptor), 5FSA, 5V5Z. This script turns
the two new deposits into docking receptors in the SAME coordinate frame as 5TZ1, so the
hash-frozen docking box (protocol.yaml, unchanged) lands in each pocket.

Chain of custody, each link checkable (mirrors scripts/prep_receptor.py for 5TZ1):

    data/structures/<ID>.pdb        raw RCSB deposit (chains, waters, bound azole, heme)
        |  keep chain A protein + chain-A HEM; drop chain B, waters, the co-crystal azole
        |  Kabsch-superpose chain-A CA onto the 5TZ1 anchor (work/receptor_A.pdb) by resSeq
        v
    work/receptor_<ID>_A.pdb        tracked frozen anchor; the gate reads its heme-iron
        |                           position, and its sha256 is pinned in the prereg
        |  obabel -xr -p 7.4  (rigid receptor, polar H at pH 7.4, AutoDock atom types)
        v
    work/receptor_<ID>.pdbqt        the docking input scripts/screen.sh reads

Two things make this honest rather than a free knob:
  * the conformer set is EXHAUSTIVE (all three C. albicans CYP51 structures), so there is
    no conformer to select;
  * the superposition only places the shared box — the N-Fe distance the gate measures is
    computed in each conformer's own frame, against that conformer's own iron.

The alignment is deterministic (numpy SVD Kabsch on the shared-resSeq CA set), so the anchor
and pdbqt are byte-reproducible and hash-pinnable exactly like 5TZ1's.

Usage:
    python scripts/prep_receptor_ensemble.py 5FSA 5V5Z
"""
import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from scripts.prep_receptor import normalize_name_remark, sha256  # reuse, don't duplicate

ANCHOR_5TZ1 = os.path.join(ROOT, "work", "receptor_A.pdb")
CHAIN = "A"
HETERO_KEEP = {"HEM"}          # keep the heme cofactor; drop waters + the co-crystal azole


def ca_coords(lines):
    """{resSeq+iCode -> (x,y,z)} for chain-A CA atoms; first altloc wins on duplicates."""
    out = {}
    for l in lines:
        if l[:6].strip() == "ATOM" and l[21] == CHAIN and l[12:16].strip() == "CA":
            key = l[22:27]                        # resSeq (22-25) + insertion code (26)
            if key not in out:                    # first altloc only, deterministic
                out[key] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return out


def kabsch(P, Q):
    """Rotation R and translation t mapping P onto Q (least-squares). Returns (R, t, rmsd)."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = Q.mean(0) - R @ P.mean(0)
    aligned = (R @ P.T).T + t
    rmsd = float(np.sqrt(((aligned - Q) ** 2).sum(1).mean()))
    return R, t, rmsd


def extract_chain_a(pdb_path):
    """Chain-A protein ATOM + chain-A HEM HETATM, original RCSB lines (co-crystal azole and
    waters dropped). Raises if no chain-A heme iron is present."""
    kept = []
    for line in open(pdb_path):
        rec = line[:6].strip()
        if rec == "ATOM" and line[21] == CHAIN:
            kept.append(line)
        elif rec == "HETATM" and line[21] == CHAIN and line[17:20].strip() in HETERO_KEEP:
            kept.append(line)
    if not any(l[:6].strip() == "HETATM" and l[12:16].strip() == "FE" for l in kept):
        sys.exit(f"ERROR: no chain-{CHAIN} heme FE in {pdb_path} — wrong chain/structure?")
    return kept


def apply_transform(lines, R, t):
    """Rewrite only the x/y/z columns (31-54) of each ATOM/HETATM line; everything else
    (atom name, resname, occupancy, b-factor, element) is preserved byte-for-byte."""
    out = []
    for l in lines:
        xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        x, y, z = R @ xyz + t
        out.append(f"{l[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{l[54:]}")
    return out


def prep_one(conf_id, ph="7.4"):
    raw = os.path.join(ROOT, "data", "structures", f"{conf_id}.pdb")
    anchor_out = os.path.join(ROOT, "work", f"receptor_{conf_id}_A.pdb")
    pdbqt_out = os.path.join(ROOT, "work", f"receptor_{conf_id}.pdbqt")
    if not os.path.exists(raw):
        sys.exit(f"ERROR: {raw} not found — download it from RCSB first")

    print(f"\n=== {conf_id} ===")
    print(f"input : data/structures/{conf_id}.pdb  (sha256 {sha256(raw)[:12]})")

    kept = extract_chain_a(raw)
    n_hem = sum(1 for l in kept if l[17:20].strip() == "HEM")
    print(f"  extracted chain {CHAIN}: {len(kept)} atoms ({len(kept)-n_hem} protein, {n_hem} heme)")

    # Superpose this conformer's CA onto the 5TZ1 anchor, over residues present in BOTH.
    ref_ca = ca_coords([l for l in open(ANCHOR_5TZ1)])
    mov_ca = ca_coords(kept)
    shared = sorted(set(ref_ca) & set(mov_ca))
    if len(shared) < 100:
        sys.exit(f"ERROR: only {len(shared)} shared CA residues with 5TZ1 — cannot align")
    P = np.array([mov_ca[k] for k in shared])
    Q = np.array([ref_ca[k] for k in shared])
    R, t, rmsd = kabsch(P, Q)
    print(f"  aligned onto 5TZ1 frame: {len(shared)} shared CA, backbone RMSD {rmsd:.3f} A")

    aligned = apply_transform(kept, R, t)
    with open(anchor_out, "w") as fh:
        fh.writelines(aligned)
        fh.write("END\n")
    print(f"  anchor: work/receptor_{conf_id}_A.pdb  (sha256 {sha256(anchor_out)})")

    # obabel needs a file; convert the aligned chain-A PDB to a rigid receptor PDBQT.
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as tmp:
        tmp.writelines(aligned)
        tmp.write("END\n")
        tmp_path = tmp.name
    try:
        proc = subprocess.run(["obabel", tmp_path, "-O", pdbqt_out, "-xr", "-p", ph],
                              capture_output=True, text=True)
    finally:
        os.unlink(tmp_path)
    if proc.returncode != 0 or not os.path.exists(pdbqt_out) or os.path.getsize(pdbqt_out) == 0:
        sys.exit(f"ERROR: obabel conversion failed\n{proc.stderr}")
    normalize_name_remark(pdbqt_out, raw)

    # A receptor is useless if the conversion dropped the iron every gate distance is measured to.
    if not any(l[:6].strip() in ("ATOM", "HETATM") and l[12:16].strip() == "FE"
               for l in open(pdbqt_out)):
        sys.exit("ERROR: no heme FE atom in the output receptor — conversion dropped the iron")
    print(f"  output: work/receptor_{conf_id}.pdbqt  (sha256 {sha256(pdbqt_out)})")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("conformers", nargs="+", help="PDB IDs, e.g. 5FSA 5V5Z")
    ap.add_argument("--ph", default="7.4", help="protonation pH for obabel (default: 7.4)")
    args = ap.parse_args()
    if not os.path.exists(ANCHOR_5TZ1):
        sys.exit(f"ERROR: {ANCHOR_5TZ1} missing — run scripts/prep_receptor.py first")
    import shutil
    if shutil.which("obabel") is None:
        sys.exit("ERROR: obabel not on PATH. Run inside the pinned env: conda activate openafr")
    for conf_id in args.conformers:
        prep_one(conf_id, args.ph)
    print("\ndone. Pin the anchor sha256s in work/PREREGISTRATION_ensemble.md, then freeze it.")


if __name__ == "__main__":
    main()
