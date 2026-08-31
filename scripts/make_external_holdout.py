"""Build the issue-#82 external validation holdout — a fresh active/decoy set no gate has
ever scored, so the frozen geometric criterion can be tested on data it has never touched.

Pipeline stage: ligands/ (external holdout)

The overfitting objection issue #82 closes: every AUC/EF the project reports was measured on
sets the criterion (the N-Fe < 3.0 A iron-approach filter, its cutoff, and the frozen run-2
protocol) was developed and re-examined against. This builds ONE clean holdout, disjoint from
all of it, and grades the FROZEN protocol on it with no re-tuning (scripts/validate_gate_external.py).

Discipline (mirrors data/ligands/PROVENANCE.md — only the provenance of "held out" changes):

  actives  a deterministic seed-82 random sample of the measured-active C. albicans azoles in
           verified_actives.smi that are NOT in verified_actives_sample.smi — i.e. the ~3,400
           actives no gate ever docked (the ensemble-confirm look #9 scored only the seed-42
           n=200 sample). Same curated population, but every molecule here is genuinely
           never-scored. Domain re-asserted: <= 45 heavy atoms, >= 1 aromatic nitrogen.

  decoys   property-matched presumed-inactives, built by the SAME frozen rule as decoys.smi /
           verified_inactives.smi (openafr.inactives.select_inactives: MW +/-25, cLogP +/-1.5,
           rotB +/-2, HBD +/-1, HBA +/-2; >= 1 aromatic N; Morgan r=2 Tanimoto < 0.35 vs every
           external active; <= 45 heavy), fetched fresh from ChEMBL under seed 82, and disjoint
           by ChEMBL id from EVERY set the project has used (all actives, decoys, verified
           inactives, potency sets) so no decoy can secretly be a known active or a reused decoy.

Deterministic: same verified_actives population + same ChEMBL cache -> byte-identical .smi.

Usage:
    python scripts/make_external_holdout.py actives            # offline, no network
    python scripts/make_external_holdout.py decoys             # fetches ChEMBL
"""
import argparse
import hashlib
import pathlib
import random
import sys
import time
import urllib.parse
import urllib.request

from rdkit import Chem, RDLogger

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import inactives

RDLogger.DisableLog("rdApp.*")

VERIFIED = "data/ligands/verified_actives.smi"
SCORED = "data/ligands/verified_actives_sample.smi"       # the seed-42 n=200 already docked
ACTIVES_OUT = "data/ligands/external_actives.smi"
DECOYS_OUT = "data/ligands/external_decoys.smi"
SEED = 82                                                 # the issue number; NOT 42 (the tuning seed)
N_ACTIVES = 50
PER_ACTIVE = 8

# every set the criterion has ever seen — decoys must be disjoint from all of them by id
USED_ID_FILES = [
    "data/ligands/verified_actives.smi", "data/ligands/verified_inactives.smi",
    "data/ligands/decoys.smi", "data/ligands/decoys_holdout.smi",
    "data/ligands/potency_enzyme.smi", "data/ligands/potency_wholecell.smi",
]
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/molecule"


def _read_smi(path):
    """[(smiles, name)] from a tab-separated .smi."""
    rows = []
    for line in open(path):
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 2:
            rows.append((parts[0], parts[1]))
    return rows


def _id_of(name):
    """Pull the CHEMBL id out of a set name like active_CHEMBL1327 / decoy_CHEMBL6362_for_x."""
    for tok in name.replace("_", " ").split():
        if tok.startswith("CHEMBL"):
            return tok
    return None


def _sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def cmd_actives(args):
    pool = _read_smi(VERIFIED)
    scored_ids = {_id_of(n) for _s, n in _read_smi(SCORED)}
    disjoint = [(smi, name) for smi, name in pool if _id_of(name) not in scored_ids]
    print(f"verified actives         : {len(pool)}")
    print(f"already scored (seed-42) : {len(scored_ids)}")
    print(f"never-scored disjoint    : {len(disjoint)}")

    # deterministic: sort by id, then a fixed-seed sample (independent of file order)
    disjoint.sort(key=lambda r: int(_id_of(r[1])[6:]) if _id_of(r[1])[6:].isdigit() else 10**12)
    rng = random.Random(SEED)
    idx = sorted(rng.sample(range(len(disjoint)), N_ACTIVES))
    chosen = [disjoint[i] for i in idx]

    # re-assert the frozen applicability domain (verified_actives already enforces it; verify)
    for smi, name in chosen:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            sys.exit(f"unparseable active {name}")
        if mol.GetNumHeavyAtoms() > inactives.MAX_HEAVY:
            sys.exit(f"{name} over domain ({mol.GetNumHeavyAtoms()} heavy) — should be impossible")
        if not inactives.has_aromatic_nitrogen(mol):
            sys.exit(f"{name} has no aromatic nitrogen — should be impossible")

    with open(ACTIVES_OUT, "w") as fh:
        for smi, name in chosen:
            fh.write(f"{smi}\t{name}\n")
    heavy = [Chem.MolFromSmiles(s).GetNumHeavyAtoms() for s, _ in chosen]
    print(f"\nwrote {len(chosen)} external actives (seed {SEED}) to {ACTIVES_OUT}")
    print(f"heavy atoms: min {min(heavy)}  mean {sum(heavy)/len(heavy):.1f}  max {max(heavy)}")
    print(f"sha256 {_sha(ACTIVES_OUT)}")
    return 0


def _fetch_pool(mw_lo, mw_hi, offset, limit=1000):
    q = urllib.parse.urlencode({
        "molecule_properties__mw_freebase__gte": f"{mw_lo:.1f}",
        "molecule_properties__mw_freebase__lte": f"{mw_hi:.1f}",
        "limit": limit, "offset": offset, "format": "json",
    })
    try:
        with urllib.request.urlopen(f"{CHEMBL}?{q}", timeout=90) as r:
            import json
            data = json.load(r)
    except Exception as e:
        print(f"    (fetch failed at offset {offset}: {e})")
        return []
    out = []
    for m in data.get("molecules", []):
        smi = (m.get("molecule_structures") or {}).get("canonical_smiles")
        if smi:
            out.append((m["molecule_chembl_id"], smi))
    return out


def cmd_decoys(args):
    actives = inactives.load_actives(ACTIVES_OUT)
    active_fps = [a["fp"] for a in actives]
    mws = [a["mw"] for a in actives]
    print(f"external actives         : {len(actives)}  (MW {min(mws):.0f}-{max(mws):.0f})")

    used = set()
    for f in USED_ID_FILES:
        used |= {_id_of(n) for _s, n in _read_smi(f)}
    used |= {_id_of(a["name"]) for a in actives}
    print(f"excluded ids (all sets)  : {len(used)}")

    # Fetch fresh ChEMBL molecules across the actives' MW band. Deterministic: fixed offsets,
    # candidates de-duplicated and sorted by id before matching, so the pool order never
    # depends on network timing. select_inactives then applies the frozen property/novelty rule.
    band_lo, band_hi = min(mws) - inactives.MW_TOL, max(mws) + inactives.MW_TOL
    cand = {}
    for offset in range(0, args.pages * 1000, 1000):
        pool = _fetch_pool(band_lo, band_hi, offset)
        time.sleep(0.2)
        for cid, smi in pool:
            if cid not in used and cid not in cand:
                cand[cid] = smi
        print(f"    offset {offset:5d}: pool {len(pool):4d}  distinct usable so far {len(cand)}")
        if not pool:
            break
    candidates = sorted(cand.items(), key=lambda kv: int(kv[0][6:]) if kv[0][6:].isdigit() else 10**12)
    print(f"candidate pool           : {len(candidates)} distinct, disjoint by id")

    chosen, rejects = inactives.select_inactives(candidates, actives, active_fps, PER_ACTIVE)
    print("\nfilter funnel:")
    for k, v in sorted(rejects.items(), key=lambda kv: -kv[1]):
        print(f"  rejected {k:20} {v:6d}")
    print(f"  -> KEPT {len(chosen)} property-matched decoys")

    chosen.sort(key=lambda c: int(c["chembl_id"][6:]) if c["chembl_id"][6:].isdigit() else 10**12)
    with open(DECOYS_OUT, "w") as fh:
        for c in chosen:
            fh.write(f"{c['smiles']}\tdecoy_{c['chembl_id']}_for_{c['paired']}\n")
    print(f"\nwrote {len(chosen)} external decoys to {DECOYS_OUT}")
    print(f"sha256 {_sha(DECOYS_OUT)}")

    if chosen:
        print("\nproperty match (external decoys vs external actives):")
        for key, label in (("mw", "MW"), ("logp", "cLogP"), ("heavy", "heavy atoms"),
                           ("rotb", "rot. bonds")):
            av = sum(a[key] for a in actives) / len(actives)
            dv = sum(c[key] for c in chosen) / len(chosen)
            print(f"  {label:12} actives {av:8.1f}   decoys {dv:8.1f}   diff {abs(av-dv):6.1f}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("actives", help="offline: sample the never-scored disjoint actives")
    a.set_defaults(func=cmd_actives)
    d = sub.add_parser("decoys", help="fetch ChEMBL + build property-matched disjoint decoys")
    d.add_argument("--pages", type=int, default=12, help="1000-molecule ChEMBL pages to scan")
    d.set_defaults(func=cmd_decoys)
    args = ap.parse_args(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    sys.exit(main())
