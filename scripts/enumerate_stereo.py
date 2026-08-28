"""Stereochemistry enumeration + ligand prep for the per-candidate stereo axis (#85).

Deterministic, network-free, docking-free. For every shortlist molecule
(`work/candidates/shortlist_focus.smi`) it reports, straight from the supplied SMILES:

  - tetrahedral stereocenters (assigned + unassigned), via RDKit
    FindMolChiralCenters(includeUnassigned=True)
  - unspecified double-bond (E/Z) stereo
  - the number of distinct stereoisomers implied by the UNASSIGNED centers
    (EnumerateStereoisomers, onlyUnassigned=True, unique=True), capped
  - stereo_explicit = the supplied SMILES pins a SINGLE isomer (zero unassigned
    tetrahedral centers AND zero unassigned stereo bonds). When true, the shortlist
    rank already refers to one defined isomer and needs no isomer docking.

For the molecules that are NOT stereo-explicit it writes one embedded 3D SDF per
enumerated isomer into --ligand-dir (named <candidate>__isoNN.sdf, plus a
stereo_map.json recording each isomer's isomeric SMILES + CIP labels), so the
consensus screen can dock every isomer under the frozen protocol. The grader
scripts/candidate_stereo.py reads those per-isomer dock trees.

This file produces NO nitrogen-to-iron geometry; it is safe to run before the
pre-registration is frozen (the enumeration is a deterministic function of the
input SMILES, not a gradeable result).

Usage:
    python scripts/enumerate_stereo.py \
        [--smi work/candidates/shortlist_focus.smi] \
        [--ligand-dir work/candidates/stereo_ligands] \
        [--audit-out work/candidates/stereo_audit.json] \
        [--max-isomers 16] [--json]
"""
import argparse
import json
import os
import pathlib
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import (EnumerateStereoisomers,
                                               StereoEnumerationOptions)

ROOT = pathlib.Path(__file__).resolve().parent.parent


def read_smi(path):
    """[(smiles, name)] from a whitespace/tab SMILES file (name = last field)."""
    out = []
    for line in open(path):
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        # SMILES may itself contain no whitespace; name is the final tab/space field.
        parts = line.rsplit("\t", 1) if "\t" in line else line.rsplit(None, 1)
        if len(parts) != 2:
            continue
        smi, name = parts
        out.append((smi.strip(), name.strip()))
    return out


def unassigned_stereo_bonds(mol):
    """Count double bonds with perceivable but unspecified E/Z stereo."""
    Chem.FindPotentialStereoBonds(mol)
    n = 0
    for b in mol.GetBonds():
        if b.GetStereo() == Chem.BondStereo.STEREOANY:
            n += 1
    return n


def audit_one(smi, name, max_isomers):
    """Deterministic stereo audit + enumerated isomers for one molecule."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f"unparseable SMILES for {name}: {smi}")
    centers = Chem.FindMolChiralCenters(
        mol, force=True, includeUnassigned=True, useLegacyImplementation=False)
    assigned = [c for c in centers if c[1] != "?"]
    unassigned = [c for c in centers if c[1] == "?"]
    n_bond_unassigned = unassigned_stereo_bonds(Chem.MolFromSmiles(smi))

    opts = StereoEnumerationOptions(onlyUnassigned=True, unique=True,
                                    maxIsomers=max_isomers)
    isomers = list(EnumerateStereoisomers(mol, options=opts))
    iso_records = []
    for iso in isomers:
        Chem.AssignStereochemistry(iso, cleanIt=True, force=True)
        iso_smiles = Chem.MolToSmiles(iso)
        labels = Chem.FindMolChiralCenters(
            iso, force=True, includeUnassigned=False, useLegacyImplementation=False)
        iso_records.append({"smiles": iso_smiles,
                            "cip": {f"a{a}": lab for a, lab in labels}})
    # Canonical ordering so isoNN indices are stable across runs.
    iso_records.sort(key=lambda r: r["smiles"])

    stereo_explicit = (len(unassigned) == 0 and n_bond_unassigned == 0)
    return {
        "name": name,
        "smiles": smi,
        "n_stereocenters": len(centers),
        "n_assigned": len(assigned),
        "n_unassigned": len(unassigned),
        "n_unassigned_bonds": n_bond_unassigned,
        "n_isomers": len(iso_records),
        "stereo_explicit": stereo_explicit,
        "isomers": iso_records,
    }


def embed_sdf(smiles, out_path):
    """Embed one isomeric SMILES to a single 3D conformer SDF (deterministic seed)."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D  # fixed so the prepared ligand is reproducible
    if AllChem.EmbedMolecule(mol, params) != 0:
        # fall back to a plain distance-geometry embed if ETKDG fails
        if AllChem.EmbedMolecule(mol, randomSeed=0xF00D) != 0:
            raise RuntimeError(f"embed failed for {smiles}")
    AllChem.MMFFOptimizeMolecule(mol)
    w = Chem.SDWriter(str(out_path))
    w.write(mol)
    w.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smi", default="work/candidates/shortlist_focus.smi")
    ap.add_argument("--ligand-dir", default="work/candidates/stereo_ligands")
    ap.add_argument("--audit-out", default="work/candidates/stereo_audit.json")
    ap.add_argument("--max-isomers", type=int, default=16)
    ap.add_argument("--no-embed", action="store_true",
                    help="audit only; do not write isomer SDFs")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    os.chdir(ROOT)
    mols = read_smi(args.smi)
    records = [audit_one(smi, name, args.max_isomers) for smi, name in mols]

    lig_dir = pathlib.Path(args.ligand_dir)
    stereo_map = {}
    if not args.no_embed:
        lig_dir.mkdir(parents=True, exist_ok=True)
        for rec in records:
            if rec["stereo_explicit"]:
                continue  # single defined isomer — nothing to enumerate/dock
            for k, iso in enumerate(rec["isomers"], 1):
                tag = f"{rec['name']}__iso{k:02d}"
                embed_sdf(iso["smiles"], lig_dir / f"{tag}.sdf")
                stereo_map[tag] = {"candidate": rec["name"],
                                   "isomer_index": k,
                                   "smiles": iso["smiles"],
                                   "cip": iso["cip"]}
        (lig_dir / "stereo_map.json").write_text(json.dumps(stereo_map, indent=2))

    audit = {"n_molecules": len(records),
             "max_isomers": args.max_isomers,
             "candidates": records}
    pathlib.Path(args.audit_out).write_text(json.dumps(audit, indent=2))

    if args.json:
        print(json.dumps(audit, indent=2))
        return 0

    print("=" * 88)
    print("STEREOCHEMISTRY ENUMERATION AUDIT (#85)")
    print(f"smiles: {args.smi}   (source = supplied shortlist SMILES)")
    print("=" * 88)
    print("%-16s %7s %8s %10s %6s %8s  %s" % (
        "candidate", "centers", "assigned", "unassigned", "bonds", "isomers", "verdict"))
    print("-" * 88)
    for r in records:
        verdict = "EXPLICIT" if r["stereo_explicit"] else f"AMBIGUOUS ({r['n_isomers']} iso)"
        print("%-16s %7d %8d %10d %6d %8d  %s" % (
            r["name"], r["n_stereocenters"], r["n_assigned"], r["n_unassigned"],
            r["n_unassigned_bonds"], r["n_isomers"], verdict))
    ambiguous = [r["name"] for r in records if not r["stereo_explicit"]]
    print("-" * 88)
    print(f"stereo-explicit: {len(records) - len(ambiguous)}/{len(records)}   "
          f"ambiguous (isomer-docked): {ambiguous or 'none'}")
    if not args.no_embed and stereo_map:
        print(f"wrote {len(stereo_map)} isomer SDFs -> {lig_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
