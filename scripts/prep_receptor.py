"""5TZ1.pdb -> work/receptor.pdbqt: the docking receptor, built by script.

Pipeline stage: receptor prep (the manual step that used to block a fresh clone).

The run-2 result was docked against a receptor a human prepared by hand and never
captured, so `work/receptor.pdbqt` (gitignored) did not exist on a clean checkout
and `run_screen2.sh` had nothing to dock into. This script reproduces that file
deterministically from the raw RCSB structure.

Chain of custody, each link checkable:

    data/structures/5TZ1.pdb          raw RCSB deposit (chains A+B, waters, bound VT1)
        |  keep chain A protein + chain-A HEM; drop chain B, waters, the ligand
        v
    (chain-A atom set)  ==  work/receptor_A.pdb   <-- tracked frozen anchor the gate
        |                                             reads the heme-iron position from
        |  obabel -xr -p 7.4  (rigid receptor, polar H at pH 7.4, AutoDock atom types)
        v
    work/receptor.pdbqt               the docking input run_screen2.sh reads

The middle equality is asserted, not assumed: the extracted atoms must match
work/receptor_A.pdb exactly (same atoms, order, coordinates). That is what keeps
what gets *docked* tied to the same receptor the gate *grades* against — the
concern that this whole step exists to close.

Verify a rebuilt receptor is the right one by re-docking VT-1161 and checking the
top score and iron-bound pose against work/RESULTS_redock_VT1.md:
    top affinity ~= -12.1 kcal/mol; an iron-bound pose at N-Fe ~= 2.6 A.
Vina ignores receptor partial charges (RESULTS_redock_VT1.md, finding 3), so the
all-zero charges obabel writes here do not affect the score — atom types and
geometry do, and those are what this script pins.
"""
import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_IN = os.path.join(ROOT, "data", "structures", "5TZ1.pdb")
DEFAULT_OUT = os.path.join(ROOT, "work", "receptor.pdbqt")
ANCHOR = os.path.join(ROOT, "work", "receptor_A.pdb")

CHAIN = "A"          # crystal has two copies (A, B); run-2 used chain A
HETERO_KEEP = {"HEM"}  # keep the heme cofactor; drop waters (HOH) and the bound ligand (VT1)


def atom_key(line):
    """Identity of an atom independent of its serial number (which gets renumbered).

    (record, name, altloc, resname, resSeq, x, y, z, element) — the fields that must
    be identical for two files to be the same receptor.
    """
    return (
        line[:6].strip(), line[12:16].strip(), line[16], line[17:20].strip(),
        line[22:26].strip(), line[30:38].strip(), line[38:46].strip(),
        line[46:54].strip(), line[76:78].strip(),
    )


def iron_xyz(path):
    """(x, y, z) of the heme iron, or None. Keys on atom name so it reads both the
    anchor PDB (HETATM FE) and obabel's pdbqt (ATOM FE, AutoDock type in cols 77-78)."""
    for l in open(path):
        if l[:6].strip() in ("ATOM", "HETATM") and l[12:16].strip() == "FE":
            return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return None


def extract_chain_a(pdb_path):
    """Return the original 5TZ1 lines for chain-A protein + chain-A heme, in file order."""
    kept = []
    for line in open(pdb_path):
        rec = line[:6].strip()
        if rec == "ATOM" and line[21] == CHAIN:
            kept.append(line)
        elif rec == "HETATM" and line[21] == CHAIN and line[17:20].strip() in HETERO_KEEP:
            kept.append(line)
    return kept


def verify_against_anchor(kept):
    """Hard-check the extracted atoms match the tracked work/receptor_A.pdb exactly."""
    if not os.path.exists(ANCHOR):
        print(f"  WARNING: {os.path.relpath(ANCHOR, ROOT)} not found; skipping anchor check")
        return
    ref = [line for line in open(ANCHOR) if line[:6].strip() in ("ATOM", "HETATM")]
    if len(kept) != len(ref):
        sys.exit(f"ERROR: extracted {len(kept)} atoms but anchor has {len(ref)} — "
                 "chain-A selection does not match work/receptor_A.pdb")
    for i, (a, b) in enumerate(zip(kept, ref)):
        if atom_key(a) != atom_key(b):
            sys.exit(f"ERROR: atom {i + 1} differs from work/receptor_A.pdb\n"
                     f"  extracted: {atom_key(a)}\n  anchor:    {atom_key(b)}")
    print(f"  anchor check: {len(kept)} atoms match work/receptor_A.pdb exactly")


def normalize_name_remark(pdbqt_path, source_pdb):
    """Rewrite obabel's `REMARK  Name = <scratch path>` to name the real input structure.

    obabel copies its input filename into the first REMARK, and the input here is a
    randomly-named temp file — so two byte-identical receptors hashed differently on every
    rebuild, and the receptor could not be hash-pinned the way protocol.yaml and the
    pre-registrations are. The atom records were always identical; only this one line moved.
    """
    lines = open(pdbqt_path).readlines()
    for i, line in enumerate(lines):
        if line.startswith("REMARK") and "Name =" in line:
            lines[i] = f"REMARK  Name = {os.path.basename(source_pdb)}\n"
            break
    with open(pdbqt_path, "w") as fh:
        fh.writelines(lines)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdb_in", nargs="?", default=DEFAULT_IN,
                    help="raw structure (default: data/structures/5TZ1.pdb)")
    ap.add_argument("pdbqt_out", nargs="?", default=DEFAULT_OUT,
                    help="docking receptor to write (default: work/receptor.pdbqt)")
    ap.add_argument("--ph", default="7.4", help="protonation pH for obabel (default: 7.4)")
    ap.add_argument("--mutant", action="store_true",
                    help="pdb_in is an already-extracted (chain A + heme) mutant receptor, "
                         "e.g. from scripts/mutate_receptor.py. Skips the wild-type "
                         "atom-identity check (a mutant is meant to differ) but still "
                         "enforces that the heme iron matches work/receptor_A.pdb.")
    args = ap.parse_args()

    if shutil.which("obabel") is None:
        sys.exit("ERROR: obabel not on PATH. Run inside the pinned env: conda activate openafr")
    if not os.path.exists(args.pdb_in):
        sys.exit(f"ERROR: input structure not found: {args.pdb_in}")

    print(f"input : {os.path.relpath(args.pdb_in, ROOT)}  (sha256 {sha256(args.pdb_in)[:12]})")

    kept = extract_chain_a(args.pdb_in)
    if not kept:
        sys.exit(f"ERROR: no chain-{CHAIN} atoms found in {args.pdb_in}")
    n_hem = sum(1 for l in kept if l[17:20].strip() == "HEM")
    print(f"  extracted chain {CHAIN}: {len(kept)} atoms ({len(kept) - n_hem} protein, {n_hem} heme)")
    if args.mutant:
        print("  --mutant: skipping wild-type atom-identity check (iron check still enforced below)")
    else:
        verify_against_anchor(kept)

    # obabel needs a file; write the extracted chain-A PDB to a temp file (original
    # RCSB lines, unmodified) and convert it to a rigid receptor PDBQT.
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as tmp:
        tmp.writelines(kept)
        tmp.write("END\n")
        tmp_path = tmp.name
    try:
        os.makedirs(os.path.dirname(os.path.abspath(args.pdbqt_out)), exist_ok=True)
        # -xr: rigid receptor (no flexible torsions).  -p: add polar H at the given pH.
        # A "Failed to kekulize aromatic bonds" warning on the heme is expected and benign.
        proc = subprocess.run(
            ["obabel", tmp_path, "-O", args.pdbqt_out, "-xr", "-p", args.ph],
            capture_output=True, text=True,
        )
    finally:
        os.unlink(tmp_path)

    if proc.returncode != 0 or not os.path.exists(args.pdbqt_out) or os.path.getsize(args.pdbqt_out) == 0:
        sys.exit(f"ERROR: obabel conversion failed\n{proc.stderr}")

    normalize_name_remark(args.pdbqt_out, args.pdb_in)

    # The receptor is useless for this pipeline if the heme iron did not survive the
    # conversion, or moved: every gate distance is measured to it.
    fe = iron_xyz(args.pdbqt_out)
    if fe is None:
        sys.exit("ERROR: no heme FE atom in the output receptor — conversion dropped the iron")
    fe_ref = iron_xyz(ANCHOR) if os.path.exists(ANCHOR) else None
    if fe_ref is not None and max(abs(a - b) for a, b in zip(fe, fe_ref)) > 1e-3:
        sys.exit(f"ERROR: heme iron moved during conversion: {fe} vs anchor {fe_ref}")

    print(f"output: {os.path.relpath(args.pdbqt_out, ROOT)}  (sha256 {sha256(args.pdbqt_out)[:12]})")
    print(f"  heme Fe at ({fe[0]:.3f}, {fe[1]:.3f}, {fe[2]:.3f})"
          + ("  [matches anchor]" if fe_ref is not None else ""))
    print("done. Verify with a VT-1161 redock against work/RESULTS_redock_VT1.md.")


if __name__ == "__main__":
    main()
