"""work/receptor_A.pdb + mutation spec -> a mutant heavy-atom receptor PDB.

Pipeline stage: receptor mutation (the *C. auris* resistance-pocket port). See
`.context/auris_port_scope.md` for why this exists — the gate and repurposing screen
dock against the *C. albicans* wild-type 5TZ1 pocket, but the mission is the resistant
*C. auris* pocket, for which no experimental structure exists.

The approach (Path A): graft *C. auris* ERG11 resistance point mutations onto the
already-validated 5TZ1 chain-A scaffold, changing only the mutated side chains and
nothing else — same heme, same iron, same coordinate frame, same docking box. That keeps
the entire frozen-protocol chain of custody intact and isolates the resistance mutation
as the only variable.

DETERMINISTIC BY CONSTRUCTION. The pinned environment has no side-chain repacker
(no PyMOL/FoldX/SCWRL), and guessing a rotamer would break the same reproducibility
contract as the hash-frozen protocol. So this script *only* supports mutations that are a
pure atom **deletion** from the wild-type side chain — where the mutant side chain is a
strict subset of the wild-type atoms in identical positions, so the result is unique with
no rotamer choice:

    Y132F  (Tyr -> Phe)   delete the para-hydroxyl OH; the phenyl ring is unchanged
    V125A  (Val -> Ala)   delete CG1, CG2

These cover the dominant global *C. auris* resistance mutation (Y132F, clades I & IV).
Build-class mutations that *add* atoms — K143R (Lys->Arg, clade I) and F126L
(Phe->Leu, clade III) — are REFUSED, not approximated: they need a pinned repacker, which
is a documented follow-up (`.context/auris_port_scope.md`, Path A step 2). Refusing is the
honest behaviour; a silently guessed side chain is worse than none.

The output is a heavy-atom PDB (chain A + heme), identical to the wild-type anchor
everywhere except the mutated residues. Convert it to a docking receptor with the same
obabel path the wild type uses:

    python scripts/mutate_receptor.py Y132F -o work/receptor_auris_Y132F.pdb
    python scripts/prep_receptor.py --mutant work/receptor_auris_Y132F.pdb \\
        work/receptor_auris_Y132F.pdbqt

The `--mutant` flag makes prep_receptor skip the wild-type atom-identity check (the point
of a mutant is that it differs) while still enforcing that the heme iron did not move.
"""
import argparse
import hashlib
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_ANCHOR = os.path.join(ROOT, "work", "receptor_A.pdb")

# One-letter -> three-letter for the residues that appear in supported specs.
AA1TO3 = {
    "A": "ALA", "F": "PHE", "K": "LYS", "L": "LEU", "R": "ARG", "V": "VAL", "Y": "TYR",
}

# Canonical heavy-atom names per residue — used to *assert* the wild-type residue is
# intact before mutating (assert, don't assume, same as prep_receptor's anchor check).
HEAVY_ATOMS = {
    "ALA": {"N", "CA", "C", "O", "CB"},
    "VAL": {"N", "CA", "C", "O", "CB", "CG1", "CG2"},
    "PHE": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
}

# Explicit whitelist of deterministic deletion mutations: (wt3, mut3) -> atoms to remove.
# Whitelisted, not inferred by name-subset: PHE's atom names are a subset of TYR's *and*
# LEU's names look like a subset of PHE's, but PHE->LEU is NOT a valid deletion (deleting
# the far ring atoms leaves an aromatic-geometry CG/CD, not a leucine isopropyl). Only
# pairs where every retained atom keeps both its identity and geometry belong here.
SUPPORTED = {
    ("TYR", "PHE"): {"OH"},          # Phe = Tyr minus the para-hydroxyl
    ("VAL", "ALA"): {"CG1", "CG2"},  # Ala = Val minus the two gamma methyls
}


def parse_spec(spec):
    """'Y132F' -> ('TYR', 132, 'PHE'). Errors are fatal and specific."""
    s = spec.strip().upper()
    if len(s) < 3 or not s[0].isalpha() or not s[-1].isalpha():
        sys.exit(f"ERROR: bad mutation spec {spec!r} — expected e.g. Y132F")
    wt1, mut1, mid = s[0], s[-1], s[1:-1]
    if not mid.isdigit():
        sys.exit(f"ERROR: bad mutation spec {spec!r} — residue number not numeric")
    if wt1 not in AA1TO3 or mut1 not in AA1TO3:
        sys.exit(f"ERROR: unknown amino-acid letter in spec {spec!r}")
    return AA1TO3[wt1], int(mid), AA1TO3[mut1]


def read_atoms(path):
    """All ATOM/HETATM lines of a PDB, in file order (this is chain A + heme already)."""
    if not os.path.exists(path):
        sys.exit(f"ERROR: receptor anchor not found: {path}")
    return [l for l in open(path) if l[:6].strip() in ("ATOM", "HETATM")]


def resseq(line):
    return line[22:26].strip()


def atname(line):
    return line[12:16].strip()


def resname(line):
    return line[17:20].strip()


def apply_mutation(lines, wt3, pos, mut3):
    """Return new lines with residue `pos` mutated wt3 -> mut3 by atom deletion.

    Asserts the target residue is the expected, intact wild type before touching it, and
    refuses any mutation not on the deterministic whitelist.
    """
    if (wt3, mut3) not in SUPPORTED:
        sys.exit(
            f"ERROR: {wt3}{pos}{mut3} is not a deterministic deletion mutation.\n"
            f"  Supported (atom-deletion only): "
            + ", ".join(f"{a}->{b}" for a, b in SUPPORTED)
            + ".\n  Build-class mutations that add atoms (e.g. K143R, F126L) need a pinned\n"
            "  side-chain repacker — see .context/auris_port_scope.md (Path A, step 2)."
        )

    target = [l for l in lines if resseq(l) == str(pos) and l[:6].strip() == "ATOM"]
    if not target:
        sys.exit(f"ERROR: no ATOM records for residue {pos} in chain A")
    found = {resname(l) for l in target}
    if found != {wt3}:
        sys.exit(f"ERROR: residue {pos} is {sorted(found)}, spec expects {wt3} "
                 "— wrong structure or wrong numbering")
    have = {atname(l) for l in target}
    if have != HEAVY_ATOMS[wt3]:
        missing = HEAVY_ATOMS[wt3] - have
        extra = have - HEAVY_ATOMS[wt3]
        sys.exit(f"ERROR: {wt3}{pos} is not an intact wild-type residue "
                 f"(missing {sorted(missing)}, unexpected {sorted(extra)}) — refusing to "
                 "mutate a residue that already carries hydrogens or alt-locs")

    remove = SUPPORTED[(wt3, mut3)]
    out = []
    for l in lines:
        if resseq(l) == str(pos) and l[:6].strip() == "ATOM":
            if atname(l) in remove:
                continue                      # deleted atom
            out.append(l[:17] + mut3.ljust(3) + l[20:])  # rename residue, keep coords
        else:
            out.append(l)
    return out


def iron_xyz(lines):
    for l in lines:
        if l[:6].strip() in ("ATOM", "HETATM") and atname(l) == "FE":
            return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return None


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mutations", nargs="+", help="mutation specs, e.g. Y132F V125A")
    ap.add_argument("-o", "--out", required=True, help="mutant receptor PDB to write")
    ap.add_argument("--anchor", default=DEFAULT_ANCHOR,
                    help="wild-type heavy-atom receptor (default: work/receptor_A.pdb)")
    args = ap.parse_args()

    wt = read_atoms(args.anchor)
    fe0 = iron_xyz(wt)
    if fe0 is None:
        sys.exit(f"ERROR: no heme FE in anchor {args.anchor} — not a valid receptor")
    print(f"anchor: {os.path.relpath(args.anchor, ROOT)}  "
          f"(sha256 {sha256(args.anchor)[:12]}, {len(wt)} atoms)")

    specs = [parse_spec(s) for s in args.mutations]
    positions = [p for _, p, _ in specs]
    if len(positions) != len(set(positions)):
        sys.exit("ERROR: two specs target the same residue number")

    lines = wt
    for wt3, pos, mut3 in specs:
        lines = apply_mutation(lines, wt3, pos, mut3)
        print(f"  applied {wt3}{pos}{mut3}: removed {sorted(SUPPORTED[(wt3, mut3)])}")

    # Guarantee nothing but the intended residues changed, and the iron is untouched.
    changed_pos = {str(p) for p in positions}
    for a, b in _zip_nonmutated(wt, lines, changed_pos):
        if a != b:
            sys.exit("ERROR: a non-mutated atom line changed — refusing to write")
    fe1 = iron_xyz(lines)
    if fe1 != fe0:
        sys.exit(f"ERROR: heme iron moved: {fe0} -> {fe1}")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as f:
        f.writelines(lines)
        f.write("END\n")
    print(f"output: {os.path.relpath(args.out, ROOT)}  (sha256 {sha256(args.out)[:12]}, "
          f"{len(lines)} atoms, {len(wt) - len(lines)} removed)")
    print(f"  heme Fe unchanged at ({fe0[0]:.3f}, {fe0[1]:.3f}, {fe0[2]:.3f})")
    print("done. Convert with: python scripts/prep_receptor.py --mutant "
          f"{os.path.relpath(args.out, ROOT)} <out.pdbqt>")


def _zip_nonmutated(wt, mut, changed_pos):
    """Pairs of (wt_line, mut_line) for every residue not in changed_pos, aligned by the
    mutant's own order. Atoms of mutated residues are skipped in both lists."""
    wt_keep = [l for l in wt if resseq(l) not in changed_pos]
    mut_keep = [l for l in mut if resseq(l) not in changed_pos]
    if len(wt_keep) != len(mut_keep):
        sys.exit("ERROR: non-mutated atom count changed — a mutation touched the wrong "
                 "residue")
    return zip(wt_keep, mut_keep)


if __name__ == "__main__":
    main()
