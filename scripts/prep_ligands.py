"""SMILES -> 3D structure -> PDBQT, ready for docking.

Pipeline stage: ligands/
Each molecule: parse -> add hydrogens -> embed a 3D conformer -> optimize -> write SDF.
Unparseable molecules are skipped, logged, and counted (never silently dropped).
"""
import sys, os
from rdkit import Chem
from rdkit.Chem import AllChem

src, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)

ok, failed = [], []
for line in open(src):
    parts = line.rstrip("\n").split("\t")
    if len(parts) != 2:
        continue
    smi, name = parts
    m = Chem.MolFromSmiles(smi)
    if m is None:
        failed.append((name, "unparseable SMILES")); continue
    m = Chem.AddHs(m)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42                     # deterministic: same input -> same 3D
    if AllChem.EmbedMolecule(m, params) != 0:
        failed.append((name, "3D embedding failed")); continue
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    except Exception as e:
        failed.append((name, f"optimization failed: {e}")); continue
    path = os.path.join(outdir, f"{name}.sdf")
    Chem.SDWriter(path).write(m)
    rb = AllChem.CalcNumRotatableBonds(m)
    ok.append((name, m.GetNumHeavyAtoms(), rb))
    print(f"  {name:15} heavy atoms {m.GetNumHeavyAtoms():3d}   rotatable bonds {rb:2d}")

print()
print(f"prepared: {len(ok)}   failed: {len(failed)}")
for n, why in failed:
    print(f"  FAILED {n}: {why}")
