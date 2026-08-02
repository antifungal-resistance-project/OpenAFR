"""Rank a screened library by the validated ranking criterion.

The screen stage (scripts/screen.sh) produces docked poses; this turns them into
a ranked candidate list. Molecules are ranked by the ONE criterion the gate
proved works on held-out actives: mode C, closest approach of a ligand nitrogen
to the heme iron (smaller = better). Vina's own affinity score is carried along
as a secondary column only — on this system it ranked worse than random, so it
does not decide order.

Per the pre-registration handling rule, a molecule that produced no usable pose
(no nitrogen near the iron) is ranked LAST, never dropped. Dropping unrankable
molecules is the easiest way to flatter a candidate list.

This does NOT re-run the gate. Screening novel molecules is only justified once
scripts/validate_gate2.py has PASSED under the frozen protocol — a rank here is a
hypothesis for a wet lab, never a hit.

Usage: python scripts/rank_candidates.py SCREEN_DIR [RECEPTOR_PDB] [-o OUT.tsv]
  SCREEN_DIR    directory of docked *.pdbqt poses (e.g. work/screen)
  RECEPTOR_PDB  receptor with the heme iron       (default work/receptor_A.pdb)
"""
import argparse
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance, read_poses, read_scores


def rank(screen_dir, receptor_pdb):
    """Return rows (name, best_dist_or_None, best_score_or_None) sorted best-first.

    best_dist is the minimum nitrogen-to-iron distance over all poses of a
    molecule. Rows with no usable distance sort last, preserving the gate's
    'unrankable ranks last' rule.
    """
    fe = iron_position(receptor_pdb)
    rows = []
    for f in sorted(os.listdir(screen_dir)):
        if not f.endswith(".pdbqt"):
            continue
        path = os.path.join(screen_dir, f)
        best_d = None
        for _num, atoms in read_poses(path):
            d = min_nitrogen_iron_distance(atoms, fe)
            if d is not None and (best_d is None or d < best_d):
                best_d = d
        sc = read_scores(path)
        rows.append((f[:-6], best_d, sc[0] if sc else None))
    # smaller distance first; None (unrankable) sorted to the very end.
    rows.sort(key=lambda r: (r[1] is None, r[1] if r[1] is not None else 0.0))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("screen_dir", help="directory of docked *.pdbqt poses")
    ap.add_argument("receptor_pdb", nargs="?", default="work/receptor_A.pdb",
                    help="receptor PDB with the heme iron (default work/receptor_A.pdb)")
    ap.add_argument("-o", "--out", help="write the ranked TSV here (default: stdout only)")
    args = ap.parse_args()

    rows = rank(args.screen_dir, args.receptor_pdb)

    lines = ["rank\tligand\tN_Fe_dist_A\tvina_kcal_mol"]
    for i, (name, d, sc) in enumerate(rows, 1):
        ds = f"{d:.3f}" if d is not None else "no_pose"
        ss = f"{sc:.1f}" if sc is not None else "NA"
        lines.append(f"{i}\t{name}\t{ds}\t{ss}")
    text = "\n".join(lines)

    print(text)
    print(f"\n# {len(rows)} molecules ranked by mode C (N-Fe distance, validated criterion).",
          file=sys.stderr)
    print("# Screening is only justified after scripts/validate_gate2.py PASSES.",
          file=sys.stderr)
    print("# A high rank is a hypothesis for a wet lab, never a hit.", file=sys.stderr)
    if args.out:
        pathlib.Path(args.out).write_text(text + "\n")
        print(f"# wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
