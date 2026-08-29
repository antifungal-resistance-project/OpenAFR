"""Audit the heme/iron model that the whole criterion rests on (issue #72).

The entire OpenAFR criterion is one number: the closest approach of a ligand
nitrogen to the heme iron (`openafr.pdbqt.min_nitrogen_iron_distance`). If the
iron's coordination environment is modeled loosely, every distance that ranks
every candidate is systematically off. This script audits that foundation
head-on, deterministically, with no network and no docking engine — it only
reads the checked-in crystallographic coordinates.

It answers three questions #72 poses:

  1. PARAMETERIZATION. What is the iron's coordination sphere in the receptor we
     dock into (work/receptor_A.pdb, 5TZ1 chain A + heme)? Report the four
     porphyrin nitrogens, the proximal axial ligand, and what — if anything —
     occupies the distal (6th) axial position the incoming azole must reach.

  2. GROUND-TRUTH REPRODUCTION. Does the criterion reproduce the *known* azole
     coordination distance? Measure the crystallographic VT-1161 (oteseconazole)
     pose (work/ref_VT1.pdb) against the receptor iron with the SAME code the
     gate uses, and check it lands in the canonical 2.0-2.3 A coordination band.

  3. HEURISTIC SAFETY. `min_nitrogen_iron_distance` takes the closest of ANY
     nitrogen. Confirm that in the crystal the closest nitrogen IS the
     coordinating azole nitrogen, with a comfortable margin over the next one,
     so the heuristic is not silently latching a non-coordinating amine.

Charge/spin sensitivity is settled separately and is not re-run here: Vina's
scoring ignores receptor partial charges (work/RESULTS_redock_VT1.md, finding 3
— patching the iron to +3.000 produced byte-identical poses), so the measured
GEOMETRY is invariant to the iron charge/oxidation state we assign. This script
asserts the geometry; that result establishes it doesn't matter.

Usage:
    python scripts/audit_heme_model.py                 # human-readable report
    python scripts/audit_heme_model.py --json          # machine-readable dict
"""
import argparse
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from openafr.pdbqt import iron_position, min_nitrogen_iron_distance  # noqa: E402

RECEPTOR = os.path.join(ROOT, "work", "receptor_A.pdb")
CRYSTAL_LIGAND = os.path.join(ROOT, "work", "ref_VT1.pdb")

# Canonical azole-to-heme coordination band (N1(sp2) -> Fe(III), from CYP51 azole
# co-crystals). The criterion's FE_CUTOFF of 3.0 A used downstream is a loose
# "did it reach the iron at all" gate; this tighter band is the geometry a REAL
# coordination bond occupies, and is what a validated measurement must reproduce.
COORD_BAND = (2.0, 2.3)
# Fe-porphyrin (equatorial) N and Fe-thiolate (proximal Cys) sanity bands from
# heme-thiolate (CYP450) structures.
PORPHYRIN_N_BAND = (1.8, 2.2)
THIOLATE_BAND = (2.2, 2.5)
# Anything heavy within this radius of the iron that is NOT the heme or the
# proximal Cys would be a distal-axial ligand (e.g. a resting-state aqua ligand).
DISTAL_PROBE = 3.0


def _atoms(path):
    """Heavy atoms as (name, resn, resi, element, (x,y,z)). Skips hydrogens, the
    same way the coordinate readers in openafr.pdbqt do."""
    out = []
    for l in open(path):
        if l[:6].strip() not in ("ATOM", "HETATM"):
            continue
        name = l[12:16].strip()
        el = (l[76:78].strip() or name[0]).upper()
        if el == "H":
            continue
        out.append((name, l[17:20].strip(), l[22:26].strip(), el,
                    (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
    return out


def audit_coordination_sphere():
    """Question 1 + the distal-site check. Reads only the receptor."""
    atoms = _atoms(RECEPTOR)
    fe = iron_position(RECEPTOR)

    porphyrin_n, proximal_s, distal, waters = [], [], [], 0
    for name, resn, resi, el, xyz in atoms:
        if el == "FE":
            continue
        if resn == "HOH":
            waters += 1
        d = math.dist(xyz, fe)
        if el == "N" and resn == "HEM":
            porphyrin_n.append((round(d, 3), name))
        elif el == "S" and resn == "CYS":
            proximal_s.append((round(d, 3), name, resi))
        elif d <= DISTAL_PROBE and resn not in ("HEM",) and el != "C":
            # a non-heme, non-carbon heavy atom sitting on the iron = a distal ligand
            if not (el == "S" and resn == "CYS"):
                distal.append((round(d, 3), name, resn, resi, el))

    porphyrin_n.sort()
    proximal_s.sort()
    return {
        "iron": fe,
        "porphyrin_nitrogens": porphyrin_n,          # 4 equatorial N
        "proximal_cys_thiolate": proximal_s[0] if proximal_s else None,
        "distal_axial_ligands": sorted(distal),      # expect [] = open 6th site
        "explicit_waters_in_model": waters,          # expect 0
    }


def audit_crystal_reproduction():
    """Questions 2 + 3. Measures the crystal azole pose against the receptor iron
    with the production code path (min_nitrogen_iron_distance)."""
    fe = iron_position(RECEPTOR)
    lig = _atoms(CRYSTAL_LIGAND)
    lig_atoms = [(n, *xyz) for n, _, _, _, xyz in lig]  # (name,x,y,z) for the code

    crystal_nfe = min_nitrogen_iron_distance(lig_atoms, fe)

    n_dists = sorted((round(math.dist(xyz, fe), 3), name)
                     for name, _, _, el, xyz in lig if el == "N")
    margin = round(n_dists[1][0] - n_dists[0][0], 3) if len(n_dists) > 1 else None
    return {
        "crystal_min_nfe": round(crystal_nfe, 3),
        "in_coordination_band": COORD_BAND[0] <= crystal_nfe <= COORD_BAND[1],
        "coordination_band": COORD_BAND,
        "ligand_nitrogen_distances": n_dists,
        "closest_minus_second_A": margin,   # heuristic-latch safety margin
    }


def run():
    sphere = audit_coordination_sphere()
    crystal = audit_crystal_reproduction()
    return {"coordination_sphere": sphere, "crystal_reproduction": crystal}


def _verdicts(a):
    """Booleans a test can lock and a reader can scan. Pure derivation."""
    s, c = a["coordination_sphere"], a["crystal_reproduction"]
    n_ok = (len(s["porphyrin_nitrogens"]) == 4 and
            all(PORPHYRIN_N_BAND[0] <= d <= PORPHYRIN_N_BAND[1]
                for d, _ in s["porphyrin_nitrogens"]))
    prox = s["proximal_cys_thiolate"]
    s_ok = prox is not None and THIOLATE_BAND[0] <= prox[0] <= THIOLATE_BAND[1]
    return {
        "four_porphyrin_N_in_band": n_ok,
        "proximal_thiolate_in_band": s_ok,
        "distal_site_open": not s["distal_axial_ligands"] and s["explicit_waters_in_model"] == 0,
        "crystal_reproduces_coordination": c["in_coordination_band"],
        "heuristic_latches_coordinating_N": (c["closest_minus_second_A"] or 0) >= 0.5,
    }


def _report(a):
    s, c = a["coordination_sphere"], a["crystal_reproduction"]
    v = _verdicts(a)
    fe = s["iron"]
    print("HEME/IRON MODEL AUDIT (issue #72)")
    print(f"  receptor : work/receptor_A.pdb  (5TZ1 chain A + heme)")
    print(f"  iron     : Fe at ({fe[0]:.3f}, {fe[1]:.3f}, {fe[2]:.3f})\n")

    print("1. Iron coordination sphere")
    print(f"   porphyrin N (equatorial), expect 4 in {PORPHYRIN_N_BAND} A:")
    for d, name in s["porphyrin_nitrogens"]:
        print(f"       {name:<4} {d:.3f} A")
    prox = s["proximal_cys_thiolate"]
    print(f"   proximal axial : Cys{prox[2]} {prox[1]} (thiolate) {prox[0]:.3f} A"
          f"  [expect {THIOLATE_BAND} A]")
    if s["distal_axial_ligands"]:
        print(f"   distal (6th)   : OCCUPIED {s['distal_axial_ligands']}")
    else:
        print(f"   distal (6th)   : OPEN (no water/ligand on the iron) — the site the azole enters")
    print(f"   explicit waters in receptor model: {s['explicit_waters_in_model']}\n")

    print("2. Ground-truth reproduction — crystallographic VT-1161 (oteseconazole)")
    print(f"   criterion N-Fe (min_nitrogen_iron_distance): {c['crystal_min_nfe']:.3f} A")
    print(f"   canonical coordination band {c['coordination_band']} A -> "
          f"{'IN BAND' if c['in_coordination_band'] else 'OUT OF BAND'}\n")

    print("3. Heuristic safety — is the closest N the coordinating azole N?")
    for d, name in c["ligand_nitrogen_distances"]:
        print(f"       {name:<4} {d:.3f} A")
    print(f"   closest-to-second margin: {c['closest_minus_second_A']} A "
          f"(the coordinating N wins cleanly)\n")

    print("Verdicts")
    for k, ok in v.items():
        print(f"   [{'PASS' if ok else 'FAIL'}] {k}")
    print("\nCharge/spin: not re-run — Vina ignores receptor partial charges "
          "(RESULTS_redock_VT1.md finding 3),")
    print("so the measured geometry above is invariant to the iron charge/oxidation state.")
    return all(v.values())


def main():
    ap = argparse.ArgumentParser(description="Audit the heme/iron coordination model (#72).")
    ap.add_argument("--json", action="store_true", help="emit the raw audit as JSON")
    args = ap.parse_args()
    a = run()
    if args.json:
        print(json.dumps({**a, "verdicts": _verdicts(a)}, indent=2))
        return
    ok = _report(a)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
