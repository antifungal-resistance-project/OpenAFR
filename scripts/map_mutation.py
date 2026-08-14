"""Mutation -> pocket mapping CLI for the C. auris azole early-warning track (issue #23).

Turns a flagged ERG11 substitution into a structural verdict against the modeled
azole pocket: where the residue is, whether it contacts the drug, how far it sits
from the catalytic iron, and whether the mutant can be built and docked. All method
logic lives in openafr.mapping; this is the CLI + the wet-lab-facing wiring.

Three subcommands
-----------------
  # Map ad-hoc tokens (the three canonical C. auris ERG11 mutations + a non-residue):
  python scripts/map_mutation.py map Y132F K143R F126L TR34

  # Map the emergence detector's own output (feeds off scripts/detect_emergence.py demo):
  python scripts/map_mutation.py from-emergence

  # For deletion-buildable mutations, actually build the mutant receptor (reusing the
  # C. auris port from #13) and print the make_pose_view / render_views commands:
  python scripts/map_mutation.py map Y132F --build

This is the structural half of the loop and where the moat starts: a genomic feed can
tell you Y132F is rising, but only this step says Y132F reshapes the drug-contact
surface of the pocket while K143R sits in the second shell.
"""
import argparse
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from openafr import mapping  # noqa: E402
from openafr.pdbqt import iron_position  # noqa: E402


def _print_verdict(m):
    if not m.get("canonical", False):
        print(f"  {m['token']:10s} UNMAPPABLE -- {m['reason']}")
        return
    label = f"{m['wt']}{m['position']}{m['mut']}"
    if not m.get("mappable", False):
        print(f"  {m['token']:10s} ({label})  UNMAPPABLE -- {m['reason']}")
        return
    contact = "POCKET-CONTACT" if m["pocket_contact"] else "second-shell"
    build = "buildable(deletion)" if m["buildable"] else "needs-repacker"
    print(f"  {m['token']:10s} {contact:14s} drug {m['dist_to_ligand']:.2f} A  "
          f"iron {m['dist_to_iron']:.2f} A  [{build}]")
    print(f"      -> {m['structural_verdict']}")


def _maybe_build(m, args):
    """For a deletion-buildable, mappable mutation, build the mutant receptor and show
    the reuse of make_pose_view / render_views. Refusals stay refusals."""
    if not (args.build and m.get("mappable") and m.get("buildable")):
        return
    token = m["token"]
    out_pdb = ROOT / "work" / f"receptor_auris_{token}.pdb"
    print(f"\n  --build {token}: grafting onto the 5TZ1 chain-A port (#13)")
    r = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "mutate_receptor.py"), token,
         "-o", str(out_pdb)],
        capture_output=True, text=True)
    sys.stdout.write("      " + r.stdout.replace("\n", "\n      ").rstrip() + "\n")
    if r.returncode != 0:
        sys.stderr.write(r.stderr)
        print(f"      build FAILED -- refused, not approximated (see above)")
        return
    print(f"      view the mutant pocket with the existing tooling:")
    print(f"        python scripts/prep_receptor.py --mutant {out_pdb.relative_to(ROOT)} "
          f"work/receptor_auris_{token}.pdbqt")
    print(f"        python scripts/make_pose_view.py <azole> --receptor "
          f"{out_pdb.relative_to(ROOT)} \\")
    print(f"            --screen <docked_screen> --out work/views/{token}")
    print(f"        bash scripts/render_views.sh work/views/{token}")


def _map_tokens(tokens, args):
    residues = mapping.read_residues(args.receptor)
    ligand = mapping.read_ligand_atoms(args.ligand)
    fe = iron_position(args.receptor)
    print(f"receptor {pathlib.Path(args.receptor).name} "
          f"(residues {min(residues)}-{max(residues)}), "
          f"pocket = heavy-atom contact <= {mapping.CONTACT_CUTOFF} A of "
          f"{pathlib.Path(args.ligand).name}\n")
    verdicts = [mapping.map_mutation(t, residues, ligand, fe, args.contact_cutoff)
                for t in tokens]
    for m in verdicts:
        _print_verdict(m)
        _maybe_build(m, args)
    n_map = sum(1 for m in verdicts if m.get("mappable"))
    n_contact = sum(1 for m in verdicts if m.get("pocket_contact"))
    print(f"\n{len(verdicts)} token(s): {n_map} mapped, {n_contact} pocket-contact, "
          f"{len(verdicts) - n_map} unmappable (reported honestly).")


def cmd_map(args):
    _map_tokens(args.tokens, args)


def cmd_from_emergence(args):
    """Map the tokens the emergence demo flags, closing the detect -> map loop."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from detect_emergence import _demo_records  # reuse the #22 fixture
    from openafr import emergence as em
    result = em.detect_emergence(_demo_records(), "2024-12-31", window_days=180,
                                 min_count=3, min_delta=0.05)
    tokens = [s["mutation"] for s in result["signals"]]
    print(f"emergence demo flagged {len(tokens)} substitution(s): "
          f"{', '.join(tokens)}\n")
    _map_tokens(tokens, args)


def _add_struct_args(p):
    p.add_argument("--receptor", default=str(mapping.DEFAULT_RECEPTOR),
                   help="modeled receptor PDB (default work/receptor_A.pdb)")
    p.add_argument("--ligand", default=str(mapping.DEFAULT_LIGAND),
                   help="reference ligand defining the pocket (default work/ref_VT1.pdb)")
    p.add_argument("--contact-cutoff", type=float, default=mapping.CONTACT_CUTOFF,
                   help="heavy-atom contact distance for pocket-lining (default 4.5 A)")
    p.add_argument("--build", action="store_true",
                   help="build the mutant receptor for deletion-buildable mutations")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("map", help="map ad-hoc substitution tokens")
    p.add_argument("tokens", nargs="+", help="e.g. Y132F K143R F126L TR34")
    _add_struct_args(p)
    p.set_defaults(func=cmd_map)

    p = sub.add_parser("from-emergence",
                       help="map the tokens the emergence demo flags")
    _add_struct_args(p)
    p.set_defaults(func=cmd_from_emergence)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
