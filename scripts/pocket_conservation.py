"""Cross-species CYP51 pocket-conservation map (issue #76).

Docking happens against ONE species' CYP51 -- the *C. albicans* 5TZ1 chain-A
stand-in every gate and the frozen protocol share. This script answers the
transferability question that a single-species docking score cannot: for each
residue that lines the modeled azole pocket, is it CONSERVED in the pathogens the
mission targets (*C. auris*, *Aspergillus fumigatus*), or is it species-specific?

It is a SHARED map, not a per-candidate axis and not a gate: it never re-docks, never
promotes, and produces one pocket-conservation table plus the decisive
heme-coordination-shell read that governs EVERY candidate's N->Fe anchor. Which
individual residues a given candidate's pose touches needs the per-candidate pose
files (gitignored, regenerated off-box) -- documented as the refinement, not faked here.

Criterion + inputs are frozen in work/PREREGISTRATION_pocket_conservation.md
(sha work/PREREG_pocket_conservation.sha256). No free parameter: the 4.5 A pocket
cutoff is openafr.mapping's existing definition, the substitution scores are stock
BLOSUM62 with BLAST-default gaps, the verdict set is the three orthologs named in #76.

Usage:
    conda run -n openafr python scripts/pocket_conservation.py [--json]
Writes: work/candidates/pocket_conservation.json  and prints a readable table.
"""
import argparse
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from openafr import conservation as cz

RECEPTOR = os.path.join(ROOT, "work", "receptor_A.pdb")
LIGAND = os.path.join(ROOT, "work", "ref_VT1.pdb")

# Pinned sequences (committed with provenance + SHA-256).
ALBICANS = os.path.join(ROOT, "data", "cyp51_orthologs", "P10613.fasta")   # ERG11, structure frame
ASP_A = os.path.join(ROOT, "data", "cyp51_orthologs", "Q4WNT5.fasta")      # A. fumigatus CYP51A
ASP_B = os.path.join(ROOT, "data", "cyp51_orthologs", "Q96W81.fasta")      # A. fumigatus CYP51B
AURIS_CDS = os.path.join(ROOT, "data", "earlywarning", "erg11_reference", "erg11_cds.fasta")

OUT = os.path.join(ROOT, "work", "candidates", "pocket_conservation.json")


def load_sequences():
    alb = cz.read_fasta(ALBICANS)[1]
    aspA = cz.read_fasta(ASP_A)[1]
    aspB = cz.read_fasta(ASP_B)[1]
    auris = cz.translate(cz.read_fasta(AURIS_CDS)[1])   # from the pinned B8441 CDS
    return alb, auris, aspA, aspB


def build():
    alb, auris, aspA, aspB = load_sequences()
    result = cz.build_map(
        RECEPTOR, LIGAND, alb,
        verdict_orthologs={"auris": auris, "aspergillus_A": aspA},
        supplementary={"aspergillus_B": aspB},
    )
    result["meta"] = {
        "receptor": os.path.relpath(RECEPTOR, ROOT),
        "ligand": os.path.relpath(LIGAND, ROOT),
        "contact_cutoff_A": cz.CONTACT_CUTOFF,
        "coord_shell_A": cz.COORD_SHELL_A,
        "reference": "C. albicans ERG11 (UniProt P10613, 5TZ1 frame)",
        "verdict_orthologs": {
            "auris": "C. auris ERG11 (NCBI XM_085597798.1, B8441; translated CDS)",
            "aspergillus_A": "A. fumigatus CYP51A (UniProt Q4WNT5)",
        },
        "supplementary_orthologs": {
            "aspergillus_B": "A. fumigatus CYP51B (UniProt Q96W81) -- reported, not scored",
        },
    }
    return result


def _fmt_table(result):
    lines = []
    lines.append(f"{'res':>6} {'shell':>5} {'dFe':>6} {'dLig':>6} "
                 f"{'aur':>4} {'aspA':>4} {'aspB':>4}  verdict")
    for r in result["pocket"]:
        lines.append(
            f"{r['ref_aa']}{r['pos']:<5} "
            f"{'  *  ' if r['coord_shell'] else '     '} "
            f"{r['dist_to_iron']:>6.2f} {r['dist_to_ligand']:>6.2f} "
            f"{r['orthologs']['auris']:>4} {r['orthologs']['aspergillus_A']:>4} "
            f"{r['supplementary']['aspergillus_B']:>4}  {r['verdict']}"
            + ("  <-- transferability risk" if r["transferability_risk"] else ""))
    s = result["summary"]
    lines.append("")
    lines.append(f"pocket residues: {s['n_pocket']}  "
                 f"(IDENTICAL {s['n_identical']} / CONSERVATIVE {s['n_conservative']} "
                 f"/ DIVERGENT {s['n_divergent']})")
    lines.append(f"coordination shell (<= {cz.COORD_SHELL_A:.0f} A to iron): "
                 f"{s['n_coord_shell']} residues; divergent in shell: "
                 f"{s['coord_shell_divergent'] or 'none'}")
    lines.append(f"divergent pocket positions: {s['divergent_positions'] or 'none'}")
    lines.append("pocket-region identity to C. albicans is the point, not global identity; "
                 "global % identity: "
                 + ", ".join(f"{k} {v}%" for k, v in s["identity_to_albicans"].items()))
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", action="store_true", help="print the full JSON instead of a table")
    args = ap.parse_args()
    result = build()
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2, sort_keys=True)
        f.write("\n")
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print(_fmt_table(result))
        print(f"\nwrote {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
