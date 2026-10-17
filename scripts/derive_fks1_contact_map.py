#!/usr/bin/env python3
"""Derive the FKS1 structural contact map from the 9WZU cryo-EM structure.

Downloads (or reuses) the 9WZU mmCIF file, aligns the S. cerevisiae Fks1 sequence
around each C. auris FKS1 hot-spot window to derive position offsets, computes the
minimum heavy-atom distance from each hot-spot residue to the enfumafungin ligand
(same binding site as echinocandins — see PREREGISTRATION_fks1_structural_proximity.md),
and writes data/structures/fks1_contact_map.json.

Usage:
    python scripts/derive_fks1_contact_map.py [--pdb /path/to/9WZU.cif]

The --pdb flag allows supplying a pre-downloaded file (avoids the network).
Otherwise the script fetches 9WZU.cif from RCSB and verifies its SHA-256.

This script is meant to be run ONCE to produce the pinned contact map; the map
ships in the repo so runtime code (openafr/fks1_structural.py) never hits the
network. Re-run only when the prereg explicitly calls for a new structure.
"""
import argparse
import hashlib
import json
import math
import os
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

# ── pinned structure ──────────────────────────────────────────────────────────
PDB_ID = "9WZU"
PDB_SHA256 = "1de4c46b84f4caa958a82caa52aab1fff034d0d1b339ea10d62df5119dc45492"
PDB_URL = f"https://files.rcsb.org/download/{PDB_ID}.cif"
DRUG_COMP_ID = "A1E06"   # Enfumafungin (same FKS1 binding site as echinocandins)

# Contact-tier cutoffs (Angstroms; pre-registered in PREREGISTRATION_fks1_structural_proximity.md)
CUTOFFS = {"drug_contact": 5.0, "pocket_adjacent": 8.0, "pocket_proximate": 12.0}

OUT_PATH = os.path.join(ROOT, "data", "structures", "fks1_contact_map.json")

# ── residue codes ─────────────────────────────────────────────────────────────
AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E",
    "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


# ── mmCIF parser ──────────────────────────────────────────────────────────────
def parse_atom_site(cif_path):
    """Parse ATOM/HETATM records from the atom_site loop of an mmCIF file.

    Returns a list of dicts with keys that mirror the _atom_site column names.
    Hydrogen atoms are excluded.
    """
    cols, rows, in_loop = [], [], False
    with open(cif_path) as f:
        for raw in f:
            line = raw.rstrip()
            if line.startswith("loop_"):
                cols = []
                in_loop = True
                continue
            if in_loop and line.startswith("_atom_site."):
                cols.append(line.strip().split(".")[-1])
                continue
            if cols and (line.startswith("ATOM") or line.startswith("HETATM")):
                parts = line.split()
                if len(parts) >= len(cols):
                    row = dict(zip(cols, parts))
                    if row.get("type_symbol", "").upper() == "H":
                        continue
                    rows.append(row)
                continue
            if cols and line and not line.startswith("_") and not line.startswith("#"):
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    pass  # non-record lines inside the loop reset nothing; keep reading
    return rows


# ── Smith–Waterman local alignment (stdlib only) ─────────────────────────────
def sw_align(seq1, seq2, match=2, mismatch=-1, gap=-2):
    """Return aligned (i1, i2) index pairs (0-based) for the best local alignment.

    Insertions / deletions are represented as (None, i2) or (i1, None).
    """
    n, m = len(seq1), len(seq2)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    T = [[None] * (m + 1) for _ in range(n + 1)]
    best, bp = 0, (0, 0)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq1[i - 1] == seq2[j - 1] else mismatch
            d, u, l = H[i-1][j-1] + s, H[i-1][j] + gap, H[i][j-1] + gap
            mx = max(0, d, u, l)
            H[i][j] = mx
            if mx > 0:
                T[i][j] = "D" if mx == d else ("U" if mx == u else "L")
            if mx > best:
                best, bp = mx, (i, j)
    pairs = []
    i, j = bp
    while i > 0 and j > 0 and H[i][j] > 0:
        t = T[i][j]
        if t == "D":
            pairs.append((i - 1, j - 1))
            i -= 1
            j -= 1
        elif t == "U":
            pairs.append((i - 1, None))
            i -= 1
        else:
            pairs.append((None, j - 1))
            j -= 1
    pairs.reverse()
    return pairs


# ── contact-map derivation ────────────────────────────────────────────────────
def derive(cif_path):
    from openafr.fks1_caller import FKS1_WINDOWS, load_reference

    rows = parse_atom_site(cif_path)

    # Drug coordinates (enfumafungin heavy atoms, auth_asym_id not required since
    # label_comp_id is unique to the drug in 9WZU)
    drug_coords = [
        (float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]))
        for r in rows
        if r.get("label_comp_id") == DRUG_COMP_ID
    ]
    assert drug_coords, f"No {DRUG_COMP_ID} atoms found in {cif_path}"
    print(f"  Drug atoms ({DRUG_COMP_ID}): {len(drug_coords)}")

    # Sc residue → list of (atom_name, x, y, z, resname)
    from collections import defaultdict
    sc_atoms = defaultdict(list)
    for r in rows:
        if r.get("group_PDB") == "ATOM":
            pos = int(r["auth_seq_id"])
            sc_atoms[pos].append((
                r["label_atom_id"],
                float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]),
                AA3TO1.get(r.get("auth_comp_id", ""), "?"),
            ))

    # C. auris FKS1 reference protein (for alignment)
    _, ca_prot = load_reference()

    # For each hot-spot window, find the Ca→Sc offset via local sequence alignment
    # on a 50-residue window centred on the window.
    def compute_offset(ca_start, ca_end, sc_search_start, sc_search_end):
        ca_seq = ca_prot[ca_start - 1: ca_end + 4]
        sc_seq = "".join(
            sc_atoms[p][0][4] if sc_atoms.get(p) else "?"
            for p in range(sc_search_start, sc_search_end + 1)
        )
        pairs = sw_align(ca_seq, sc_seq)
        offsets = []
        for ci, si in pairs:
            if ci is not None and si is not None:
                ca_abs = ca_start + ci
                sc_abs = sc_search_start + si
                if ca_start <= ca_abs <= ca_end:
                    offsets.append(sc_abs - ca_abs)
        if not offsets:
            raise RuntimeError(f"No alignment found for Ca {ca_start}-{ca_end}")
        from collections import Counter
        most_common = Counter(offsets).most_common(1)[0][0]
        return most_common

    # Search windows (generous padding around each hot-spot)
    hs1_w = FKS1_WINDOWS["HS1"]
    hs2_w = FKS1_WINDOWS["HS2"]
    hs3_w = FKS1_WINDOWS["HS3"]

    offsets = {}
    for wname, w, sc_s, sc_e in [
        ("HS1", hs1_w, 620, 670),
        ("HS3", hs3_w, 680, 720),
        ("HS2", hs2_w, 1340, 1390),
    ]:
        off = compute_offset(w.start_res, w.end_res, sc_s, sc_e)
        offsets[wname] = off
        print(f"  {wname}: Ca→Sc offset = +{off}")

    # Compute minimum heavy-atom distance for each hot-spot residue
    def min_dist(sc_pos):
        atoms = sc_atoms.get(sc_pos, [])
        if not atoms:
            return None
        return min(
            math.dist((a[1], a[2], a[3]), d)
            for a in atoms
            for d in drug_coords
        )

    def tier(d):
        if d is None:
            return "no_structure"
        if d <= CUTOFFS["drug_contact"]:
            return "drug_contact"
        if d <= CUTOFFS["pocket_adjacent"]:
            return "pocket_adjacent"
        if d <= CUTOFFS["pocket_proximate"]:
            return "pocket_proximate"
        return "distal"

    residues = {}
    for wname, w in FKS1_WINDOWS.items():
        off = offsets[wname]
        for ca_pos in range(w.start_res, w.end_res + 1):
            sc_pos = ca_pos + off
            d = min_dist(sc_pos)
            t = tier(d)
            ca_aa = ca_prot[ca_pos - 1]
            sc_atoms_here = sc_atoms.get(sc_pos, [])
            sc_aa = sc_atoms_here[0][4] if sc_atoms_here else "?"
            residues[str(ca_pos)] = {
                "ca_aa": ca_aa,
                "sc_pos": sc_pos,
                "sc_aa": sc_aa,
                "min_dist_to_drug_A": round(d, 3) if d is not None else None,
                "tier": t,
                "window": wname,
            }
            flag = " ←" if t in ("drug_contact", "pocket_adjacent") else ""
            print(f"    Ca{ca_pos}({ca_aa})=Sc{sc_pos}({sc_aa})  {d:.2f}Å  {t}{flag}"
                  if d else f"    Ca{ca_pos}({ca_aa})=Sc{sc_pos}(?)  N/A  {t}")

    # Sanity checks: known panel markers must be in drug_contact or pocket_adjacent
    panel_checks = {635: "F", 639: "S", 691: "W", 1354: "R"}
    for pos, expected_aa in panel_checks.items():
        entry = residues.get(str(pos))
        assert entry, f"Panel residue Ca{pos} missing from contact map"
        assert entry["ca_aa"] == expected_aa, \
            f"Ca{pos} expected {expected_aa}, got {entry['ca_aa']} — alignment error?"
        t = entry["tier"]
        assert t in ("drug_contact", "pocket_adjacent"), \
            f"Panel residue Ca{pos}({expected_aa}) landed in unexpected tier {t!r}; " \
            f"distance {entry['min_dist_to_drug_A']}Å — alignment may be wrong"
    print("  Sanity checks passed: F635, S639, W691, R1354 all in drug_contact or pocket_adjacent")

    return {
        "source_pdb": PDB_ID,
        "source_cif_sha256": PDB_SHA256,
        "source_title": "Cryo-EM structure of Fks1 in complex with enfumafungin (9WZU, 2.97A)",
        "drug_comp_id": DRUG_COMP_ID,
        "drug_name": "Enfumafungin (same FKS1 binding site as echinocandins)",
        "contact_cutoffs_A": CUTOFFS,
        "alignment_offsets_sc_minus_ca": offsets,
        "note": (
            "Minimum heavy-atom distance from each C. auris FKS1 hot-spot residue "
            "(via S. cerevisiae alignment, offset as above) to any enfumafungin heavy atom. "
            "NOT a phenotype prediction. See PREREGISTRATION_fks1_structural_proximity.md."
        ),
        "residues": residues,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", help="Pre-downloaded 9WZU.cif path (skip network fetch)")
    args = ap.parse_args()

    cif_path = args.pdb
    if cif_path is None:
        tmp = os.path.join(os.path.dirname(OUT_PATH), f"{PDB_ID}.cif")
        print(f"Downloading {PDB_URL} → {tmp}")
        urllib.request.urlretrieve(PDB_URL, tmp)
        cif_path = tmp

    # Verify SHA-256
    print(f"Verifying SHA-256 of {cif_path}...")
    h = hashlib.sha256(open(cif_path, "rb").read()).hexdigest()
    if h != PDB_SHA256:
        sys.exit(f"SHA-256 mismatch!\n  expected: {PDB_SHA256}\n  got:      {h}")
    print(f"  OK ({h[:16]}…)")

    print("Deriving contact map…")
    contact_map = derive(cif_path)

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w") as f:
        json.dump(contact_map, f, indent=2)

    # Pin the output SHA-256
    out_sha = hashlib.sha256(open(OUT_PATH, "rb").read()).hexdigest()
    print(f"\nWrote {OUT_PATH}")
    print(f"Contact map SHA-256: {out_sha}")
    pin_path = os.path.join(ROOT, "work", "PREREG_fks1_structural_proximity.sha256")
    print(f"\nAdd to {pin_path}:")
    print(f"  fks1_contact_map.json  {out_sha}")
    print(f"  9WZU.cif               {PDB_SHA256}")


if __name__ == "__main__":
    main()
