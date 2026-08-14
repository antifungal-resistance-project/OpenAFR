"""Structural so-what CLI: azole-fit consequence of a flagged mutation (issue #24).

Turns a flagged ERG11 substitution into the actionable fit verdict the alert layer
(#25) renders: does fluconazole still fit the pocket, and how much worse. All method
logic lives in openafr.structural (which builds on openafr.mapping, #23); this is the
CLI + the human-facing wiring.

Two subcommands
---------------
  # Estimate ad-hoc tokens (the canonical C. auris ERG11 mutations + a non-residue):
  python scripts/structural_sowhat.py estimate Y132F K143R F126L V125A TR34

  # Estimate the emergence detector's own output, closing detect -> map -> estimate:
  python scripts/structural_sowhat.py from-emergence

Each verdict states its evidence tier explicitly -- measured (a real #13 dock),
estimate (deterministic direction, dock for the number), or no confident call -- so a
non-bioinformatician can see not just the answer but how much to trust it.
"""
import argparse
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from openafr import mapping, structural  # noqa: E402
from openafr.pdbqt import iron_position  # noqa: E402

_TIER = {"empirical-dock": "MEASURED", "geometric-estimate": "ESTIMATE",
         "no-call": "NO-CALL"}


def _print_verdict(v):
    tier = _TIER[v["evidence_class"]]
    token = v["token"]
    print(f"  {token:10s} [{tier:8s} conf={v['confidence']:8s} dir={v['direction']}]")
    print(f"      {v['fluconazole_fit_verdict']}")
    if v["lost_contacts"]:
        lc = ", ".join(f"{c['atom']}@{c['dist_to_drug']:.2f}A" for c in v["lost_contacts"])
        print(f"      lost drug contacts (geometry): {lc}")
    if v["empirical"]:
        e = v["empirical"]
        print(f"      dock: gate {e['gate']}  AUC {e['auc']} (wt {e['auc_wildtype']})  "
              f"EF@1% {e['ef1_wildtype']}x -> {e['ef1_mutant']}x  [{e['ligand_set']}]")
    for c in v["caveats"]:
        print(f"      - caveat: {c}")


def _run(tokens, args):
    residues = mapping.read_residues(args.receptor)
    ligand = mapping.read_ligand_atoms(args.ligand)
    fe = iron_position(args.receptor)
    print(f"receptor {pathlib.Path(args.receptor).name}, "
          f"pocket = heavy-atom contact <= {args.contact_cutoff} A of "
          f"{pathlib.Path(args.ligand).name}\n")
    verdicts = []
    for t in tokens:
        mv = mapping.map_mutation(t, residues, ligand, fe, args.contact_cutoff)
        v = structural.estimate_fit_consequence(mv, residues, ligand, args.contact_cutoff)
        _print_verdict(v)
        verdicts.append(v)
    n = len(verdicts)
    measured = sum(1 for v in verdicts if v["confidence"] == "measured")
    estimate = sum(1 for v in verdicts if v["confidence"] == "estimate")
    nocall = sum(1 for v in verdicts if v["confidence"] == "none")
    print(f"\n{n} token(s): {measured} measured, {estimate} estimated, "
          f"{nocall} no-confident-call (honesty bar held).")


def cmd_estimate(args):
    _run(args.tokens, args)


def cmd_from_emergence(args):
    """Estimate the tokens the #22 emergence demo flags -- the full loop, no hand-off."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from detect_emergence import _demo_records
    from openafr import emergence as em
    result = em.detect_emergence(_demo_records(), "2024-12-31", window_days=180,
                                 min_count=3, min_delta=0.05)
    tokens = [s["mutation"] for s in result["signals"]]
    print(f"emergence demo flagged {len(tokens)} substitution(s): "
          f"{', '.join(tokens)}\n")
    _run(tokens, args)


def _add_struct_args(p):
    p.add_argument("--receptor", default=str(mapping.DEFAULT_RECEPTOR),
                   help="modeled receptor PDB (default work/receptor_A.pdb)")
    p.add_argument("--ligand", default=str(mapping.DEFAULT_LIGAND),
                   help="reference ligand defining the pocket (default work/ref_VT1.pdb)")
    p.add_argument("--contact-cutoff", type=float, default=mapping.CONTACT_CUTOFF,
                   help="heavy-atom contact distance for pocket-lining (default 4.5 A)")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("estimate", help="estimate fit consequence for ad-hoc tokens")
    p.add_argument("tokens", nargs="+", help="e.g. Y132F K143R F126L V125A TR34")
    _add_struct_args(p)
    p.set_defaults(func=cmd_estimate)

    p = sub.add_parser("from-emergence",
                       help="estimate the tokens the emergence demo flags")
    _add_struct_args(p)
    p.set_defaults(func=cmd_from_emergence)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
