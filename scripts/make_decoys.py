"""Generate property-matched decoys for the known azoles.

Pipeline stage: ligands/ (decoy half)

Why property matching is mandatory
----------------------------------
Measured on this system: raw docking score correlates with heavy-atom count at
r = -0.91 (see work/RESULTS_all_azoles.md). If decoys are smaller than the actives,
the pipeline "enriches" by preferring heavy molecules and the validation gate passes
for entirely the wrong reason. Matching on size/lipophilicity removes that shortcut.

Why decoys must also contain an aromatic nitrogen
-------------------------------------------------
The validation gate requires an azole nitrogen within 3 A of the heme iron. Every real
azole has aromatic nitrogens. If decoys lacked them, the gate would separate actives
from decoys perfectly while really only asking "does this molecule contain a nitrogen?"
That is circular. Requiring decoys to be plausible iron-coordinators makes the test
measure binding quality instead.

Selection criteria per decoy
----------------------------
  matched   : MW +/- 25, cLogP +/- 1.5, rotatable bonds +/- 2, HBD +/- 1, HBA +/- 2
  plausible : >= 1 aromatic nitrogen (could coordinate the heme iron)
  different : Morgan(r=2) Tanimoto < 0.35 against EVERY active (not just its pair)

The last criterion is what makes a decoy a presumed non-binder: same physical profile,
different chemical skeleton.
"""
import json
import random
import sys
import time
import urllib.parse
import urllib.request

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/molecule"
MW_TOL, LOGP_TOL, ROTB_TOL, HBD_TOL, HBA_TOL = 25.0, 1.5, 2, 1, 2
MAX_TANIMOTO = 0.35
DECOYS_PER_ACTIVE = 50
AROMATIC_N = Chem.MolFromSmarts("[n]")

_fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def props(mol):
    return {
        "mw": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "rotb": Descriptors.NumRotatableBonds(mol),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "heavy": mol.GetNumHeavyAtoms(),
    }


def fetch_pool(mw_lo, mw_hi, limit=1000, offset=0):
    q = urllib.parse.urlencode({
        "molecule_properties__mw_freebase__gte": f"{mw_lo:.1f}",
        "molecule_properties__mw_freebase__lte": f"{mw_hi:.1f}",
        "limit": limit,
        "offset": offset,
        "format": "json",
    })
    try:
        with urllib.request.urlopen(f"{CHEMBL}?{q}", timeout=60) as r:
            data = json.load(r)
    except Exception as e:
        print(f"    (fetch failed: {e})")
        return []
    out = []
    for m in data.get("molecules", []):
        smi = (m.get("molecule_structures") or {}).get("canonical_smiles")
        if smi:
            out.append((m["molecule_chembl_id"], smi))
    return out


def main():
    actives = []
    for line in open("data/ligands/actives.smi"):
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 2:
            smi, name = parts
            m = Chem.MolFromSmiles(smi)
            actives.append({"name": name, "mol": m, "fp": _fpgen.GetFingerprint(m),
                            **props(m)})

    active_fps = [a["fp"] for a in actives]
    chosen, seen_smiles = [], set()
    random.seed(42)

    for a in actives:
        print(f"{a['name']}  (MW {a['mw']:.0f}, logP {a['logp']:.1f}, "
              f"rotB {a['rotb']}, heavy {a['heavy']})")
        picked, tried = [], 0
        for offset in (0, 1000, 2000, 4000, 8000):
            if len(picked) >= DECOYS_PER_ACTIVE:
                break
            pool = fetch_pool(a["mw"] - MW_TOL, a["mw"] + MW_TOL, offset=offset)
            time.sleep(0.2)
            for cid, smi in pool:
                tried += 1
                if smi in seen_smiles:
                    continue
                m = Chem.MolFromSmiles(smi)
                if m is None or not m.GetSubstructMatches(AROMATIC_N):
                    continue
                p = props(m)
                if (abs(p["logp"] - a["logp"]) > LOGP_TOL
                        or abs(p["rotb"] - a["rotb"]) > ROTB_TOL
                        or abs(p["hbd"] - a["hbd"]) > HBD_TOL
                        or abs(p["hba"] - a["hba"]) > HBA_TOL):
                    continue
                fp = _fpgen.GetFingerprint(m)
                if max(DataStructs.BulkTanimotoSimilarity(fp, active_fps)) >= MAX_TANIMOTO:
                    continue
                picked.append((cid, smi, a["name"], p))
                seen_smiles.add(smi)
                if len(picked) >= DECOYS_PER_ACTIVE:
                    break
        print(f"    screened {tried} candidates -> kept {len(picked)} decoys")
        chosen.extend(picked)

    out = sys.argv[1] if len(sys.argv) > 1 else "data/ligands/decoys.smi"
    with open(out, "w") as fh:
        for cid, smi, paired, _ in chosen:
            fh.write(f"{smi}\tdecoy_{cid}_for_{paired}\n")
    print(f"\nwrote {len(chosen)} decoys to {out}")

    if chosen:
        import statistics as st
        print("\nproperty match check (decoys vs actives):")
        for key, label in (("mw", "MW"), ("logp", "cLogP"), ("heavy", "heavy atoms"),
                           ("rotb", "rot. bonds")):
            av = st.mean(a[key] for a in actives)
            dv = st.mean(c[3][key] for c in chosen)
            print(f"  {label:12} actives {av:8.1f}   decoys {dv:8.1f}   diff {abs(av-dv):6.1f}")


if __name__ == "__main__":
    main()
