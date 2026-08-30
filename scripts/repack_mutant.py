"""work/receptor_A.pdb + an ADD-ATOM mutation spec -> a mutant heavy-atom receptor PDB,
with the new side chain built by a pinned side-chain repacker (PyMOL).

READ scripts/mutate_receptor.py FIRST. That script builds the *deletion-class* mutants
(Y132F, V125A) with ZERO modeling: the mutant side chain is a strict subset of the wild-type
atoms in identical positions, so the result is unique and no rotamer is ever chosen. It
REFUSES the build-class mutations that *add* atoms — K143R (Lys->Arg, clade I) and F126L
(Phe->Leu, clade III) — because placing the new atoms requires choosing a rotamer, and a
silently guessed side chain would break the reproducibility contract the whole project rests
on.

This script is the deliberately-separate, honestly-labeled home for those add-atom mutants.
It crosses that line ON PURPOSE and documents the cost: unlike a deletion, the resulting
receptor carries a **modeled rotamer**, not an observed geometry. To keep that honest and
reproducible it is:

  - PINNED. The repacker is PyMOL's mutagenesis wizard (`pymol-open-source=3.1.0`, pinned in
    the separate build-time env `environment-repack.yml` — it does not co-solve with the main
    env's rdkit/vina pins, and is build-time only, so it lives in its own env). No hand-placed
    atoms.
  - DETERMINISTIC. For the fixed backbone frame the wizard enumerates the backbone-dependent
    rotamer library and this script selects the **minimum-strain** rotamer (fewest steric
    clashes with the fixed pocket) — a single, reproducible choice, not a judgement call.
    Re-running the build on the same anchor yields a byte-identical PDB (verified in
    tests/test_repack_mutant.py).
  - ISOLATED + VERIFIED. Only the mutated residue's side chain may change: every other atom
    must stay bit-identical to the wild-type anchor and the heme iron must not move, or the
    build is REFUSED, not written. So the resistance mutation remains the only variable versus
    the wild-type candidate dock, exactly as for the deletion mutants.

What this is NOT: it is not an induced-fit / flexible-receptor model (the backbone and every
other side chain are frozen), and a RETAINED verdict downstream still means only "coordination
geometry not destroyed," never "beats resistance" (see work/PREREGISTRATION_candidate_resistance.md).
The modeled-rotamer mutants are labeled as such wherever they surface so they are never
conflated with the exact-geometry Y132F/V125A mutants.

Run under a PyMOL-bearing interpreter (`conda activate openafr-repack`, or any `pymol -cq`):

    pymol -cq scripts/repack_mutant.py -- K143R -o work/receptor_auris_K143R.pdb
    pymol -cq scripts/repack_mutant.py -- F126L -o work/receptor_auris_F126L.pdb
    python scripts/prep_receptor.py --mutant work/receptor_auris_K143R.pdb \
        work/receptor_auris_K143R.pdbqt
"""
import argparse
import hashlib
import os
import sys

# Resolve the repo root from THIS file's path. `pymol -cq script.py` does not set __file__ to
# the script, so fall back to argv[0] (the script path pymol passes) before __file__.
_self = next((p for p in (sys.argv[0] if sys.argv else "", __file__)
              if p and p.endswith("repack_mutant.py")), __file__)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(_self)))
DEFAULT_ANCHOR = os.path.join(ROOT, "work", "receptor_A.pdb")

AA1TO3 = {"F": "PHE", "K": "LYS", "L": "LEU", "R": "ARG"}

# Canonical heavy-atom names per residue — used to ASSERT the wild-type residue is intact
# before repacking (assert, don't assume; same discipline as mutate_receptor.py) and to check
# the repacked residue came out with the expected chemistry.
HEAVY_ATOMS = {
    "LYS": {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"},
    "ARG": {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
    "PHE": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "LEU": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"},
}

# The build-class (add-atom) mutations this repacker OWNS. These are exactly the specs
# mutate_receptor.py refuses. Kept explicit so a typo can never silently repack the wrong pair.
SUPPORTED = {("LYS", "ARG"), ("PHE", "LEU")}

# Non-mutated atoms must match the anchor to this tolerance (A). PyMOL rewrites coordinates,
# but the mutagenesis wizard perturbs nothing outside the target side chain, so the observed
# delta is 0.000; a strict bound catches any accidental global move.
COORD_TOL = 1e-3


def parse_spec(spec):
    """'K143R' -> ('LYS', 143, 'ARG'). Errors are fatal and specific."""
    s = spec.strip().upper()
    if len(s) < 3 or not s[0].isalpha() or not s[-1].isalpha() or not s[1:-1].isdigit():
        sys.exit(f"ERROR: bad mutation spec {spec!r} — expected e.g. K143R")
    wt1, mut1, pos = s[0], s[-1], int(s[1:-1])
    if wt1 not in AA1TO3 or mut1 not in AA1TO3:
        sys.exit(f"ERROR: unknown amino-acid letter in spec {spec!r}")
    wt3, mut3 = AA1TO3[wt1], AA1TO3[mut1]
    if (wt3, mut3) not in SUPPORTED:
        sys.exit(
            f"ERROR: {spec} is not an add-atom mutation this repacker owns.\n"
            f"  Supported: " + ", ".join(f"{a}->{b}" for a, b in sorted(SUPPORTED)) + ".\n"
            "  Deletion-class mutations (Y132F, V125A) belong in scripts/mutate_receptor.py,\n"
            "  which builds them with no rotamer modeling at all.")
    return wt3, pos, mut3


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def read_residue_atoms(path, pos):
    """{atom_name: (x,y,z)} for chain-A residue `pos`, plus its residue name(s)."""
    coords, resn = {}, set()
    for l in open(path):
        if l[:6].strip() == "ATOM" and l[21] == "A" and l[22:26].strip() == str(pos):
            coords[l[12:16].strip()] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            resn.add(l[17:20].strip())
    return coords, resn


def assert_intact_wildtype(anchor, pos, wt3):
    """Refuse to repack unless residue `pos` is the expected, intact wild-type residue."""
    coords, resn = read_residue_atoms(anchor, pos)
    if not coords:
        sys.exit(f"ERROR: no chain-A ATOM records for residue {pos} in {anchor}")
    if resn != {wt3}:
        sys.exit(f"ERROR: residue {pos} is {sorted(resn)}, spec expects {wt3} "
                 "— wrong structure or wrong numbering")
    have = set(coords)
    if have != HEAVY_ATOMS[wt3]:
        sys.exit(f"ERROR: {wt3}{pos} is not an intact wild-type residue "
                 f"(missing {sorted(HEAVY_ATOMS[wt3] - have)}, "
                 f"unexpected {sorted(have - HEAVY_ATOMS[wt3])}) — refusing to repack a "
                 "residue that already carries hydrogens, alt-locs or a prior mutation")


def repack(anchor, pos, mut3, out):
    """Run the PyMOL mutagenesis wizard: mutate chain-A residue `pos` to `mut3`, choosing the
    minimum-strain backbone-dependent rotamer, and write a heavy-atom PDB. Import is local so
    the module parses (and its pure helpers stay testable) without PyMOL installed."""
    from pymol import cmd
    cmd.reinitialize()
    cmd.set("retain_order", 1)          # keep the anchor's atom order in the output
    cmd.load(anchor, "rec")
    cmd.remove("hydro")                 # heavy-atom receptor, matching the anchor
    cmd.wizard("mutagenesis")
    cmd.refresh_wizard()
    w = cmd.get_wizard()
    w.set_mode(mut3)
    w.do_select(f"chain A and resi {pos}")
    strains = list(getattr(w, "bump_scores", []) or [])
    if not strains:
        cmd.set_wizard()
        sys.exit(f"ERROR: PyMOL loaded no rotamers for residue {pos}->{mut3}")
    best = strains.index(min(strains)) + 1     # minimum-strain rotamer (1-indexed state)
    cmd.frame(best)
    w.apply()
    cmd.set_wizard()
    cmd.remove("hydro")                 # the wizard adds H to score strain; drop them
    cmd.sort("rec")
    cmd.save(out, "rec")
    return best, min(strains), len(strains)


def iron_xyz(path):
    for l in open(path):
        if l[:6].strip() in ("ATOM", "HETATM") and l[12:16].strip() == "FE":
            return (round(float(l[30:38]), 3), round(float(l[38:46]), 3),
                    round(float(l[46:54]), 3))
    return None


def verify_isolated(anchor, out, pos, mut3):
    """Refuse the build unless ONLY residue `pos` changed and the iron is unmoved."""
    # Every non-mutated atom must be bit-identical (to COORD_TOL) to the anchor.
    def nonmut(path):
        d = {}
        for l in open(path):
            if l[:6].strip() in ("ATOM", "HETATM") and l[22:26].strip() != str(pos):
                d[(l[21], l[22:26].strip(), l[12:16].strip())] = (
                    float(l[30:38]), float(l[38:46]), float(l[46:54]))
        return d
    a, b = nonmut(anchor), nonmut(out)
    if set(a) != set(b):
        sys.exit("ERROR: the repack added/removed atoms outside the target residue — refusing")
    worst = max((max(abs(x - y) for x, y in zip(a[k], b[k])) for k in a), default=0.0)
    if worst > COORD_TOL:
        sys.exit(f"ERROR: a non-mutated atom moved by {worst:.4f} A (> {COORD_TOL}) — refusing")
    # The mutated residue must have come out as the expected residue with its full heavy set.
    coords, resn = read_residue_atoms(out, pos)
    if resn != {mut3} or set(coords) != HEAVY_ATOMS[mut3]:
        sys.exit(f"ERROR: repacked residue {pos} is {sorted(resn)} with atoms "
                 f"{sorted(coords)} — expected intact {mut3}")
    if iron_xyz(out) != iron_xyz(anchor):
        sys.exit(f"ERROR: heme iron moved: {iron_xyz(anchor)} -> {iron_xyz(out)}")
    return worst


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mutation", help="add-atom mutation spec, e.g. K143R or F126L")
    ap.add_argument("-o", "--out", required=True, help="mutant receptor PDB to write")
    ap.add_argument("--anchor", default=DEFAULT_ANCHOR,
                    help="wild-type heavy-atom receptor (default: work/receptor_A.pdb)")
    args = ap.parse_args(argv)

    wt3, pos, mut3 = parse_spec(args.mutation)
    if not os.path.exists(args.anchor):
        sys.exit(f"ERROR: anchor not found: {args.anchor}")
    print(f"anchor: {os.path.relpath(args.anchor, ROOT)}  (sha256 {sha256(args.anchor)[:12]})")
    assert_intact_wildtype(args.anchor, pos, wt3)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    state, strain, n = repack(args.anchor, pos, mut3, args.out)
    print(f"  repacked {wt3}{pos}{mut3}: PyMOL rotamer {state}/{n}, "
          f"minimum strain {strain:.3f} (backbone-dependent library)")

    worst = verify_isolated(args.anchor, args.out, pos, mut3)
    fe = iron_xyz(args.out)
    print(f"output: {os.path.relpath(args.out, ROOT)}  (sha256 {sha256(args.out)[:12]})")
    print(f"  only residue {pos} changed (worst non-mutated atom delta {worst:.4f} A); "
          f"heme Fe unchanged at {fe}")
    print("  MODELED-ROTAMER mutant — labeled as such downstream; not an exact-geometry build.")
    print("done. Convert with: python scripts/prep_receptor.py --mutant "
          f"{os.path.relpath(args.out, ROOT)} <out.pdbqt>")
    return 0


# `pymol -cq script.py -- ...` execs this module with __name__ == "pymol" (not "__main__"),
# so trigger on both: run as a plain script under the pinned env, or under `pymol -cq`.
if __name__ in ("__main__", "pymol"):
    sys.exit(main())
