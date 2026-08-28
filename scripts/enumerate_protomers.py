"""Protonation/tautomer enumeration + ligand prep for the per-candidate protomer axis (#84).

Deterministic, network-free, docking-free. For every shortlist molecule
(`work/candidates/shortlist_focus.smi`) it reports, straight from the supplied SMILES, the set of
protonation/tautomer microstates plausibly populated at pH 7.4:

  - PROTONATION microstates: the distinct microspecies Open Babel returns at pH 6.4, 7.4 and 8.4
    (`obabel -p`), i.e. the states obabel's own pKa model flips across the physiological band
    (pH 7.4 ± 1). The pH-7.4 microspecies is the reference — the exact form the scorecard dock
    used. No external pKa predictor is added, so the pinned environment is unchanged.
  - TAUTOMERS: RDKit `TautomerEnumerator` (default transforms, MaxTautomers=16) applied to EACH
    protonation microstate; every distinct enumerated tautomer is kept (no energy pruning — the
    conservative choice, see the pre-registration).
  - protomer_explicit = the (protonation × tautomer) product collapses to a SINGLE distinct
    canonical SMILES (== the reference): protonation is stable across 6.4–8.4 AND there is one
    tautomer. When true, the shortlist rank already refers to one relevant microstate and needs
    no extra docking.

For molecules that are NOT protomer-explicit it writes one embedded 3D SDF per microstate into
--ligand-dir (named <candidate>__stNN.sdf, plus a protomer_map.json recording each microstate's
canonical SMILES + how it arose), so the consensus screen can dock every microstate under the
frozen protocol. The grader scripts/candidate_protomer.py reads those per-microstate dock trees.

This file produces NO nitrogen-to-iron geometry; it is safe to run before the pre-registration is
frozen (the enumeration is a deterministic function of the input SMILES + obabel/RDKit, not a
gradeable result). Mirrors scripts/enumerate_stereo.py (#85), swapping stereoisomers for
protonation/tautomer microstates.

Usage:
    python scripts/enumerate_protomers.py \
        [--smi work/candidates/shortlist_focus.smi] \
        [--ligand-dir work/candidates/protomer_ligands] \
        [--audit-out work/candidates/protomer_audit.json] \
        [--max-states 16] [--no-embed] [--json]
"""
import argparse
import json
import os
import pathlib
import subprocess
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = pathlib.Path(__file__).resolve().parent.parent

# Physiological band sampled by Open Babel's pKa model. 7.4 is the reference (the pH the
# scorecard dock protonates at); ±1.0 brackets the range a group can titrate across in vivo.
PH_BAND = (6.4, 7.4, 8.4)
PH_REF = 7.4
MAX_TAUTOMERS = 16


def read_smi(path):
    """[(smiles, name)] from a whitespace/tab SMILES file (name = last field).

    Strips any CXSMILES extension (' |...|') so RDKit/obabel get a plain SMILES."""
    out = []
    for line in open(path):
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        parts = line.rsplit("\t", 1) if "\t" in line else line.rsplit(None, 1)
        if len(parts) != 2:
            continue
        smi, name = parts
        smi = smi.split(" |")[0].strip()  # drop CXSMILES stereo annotation if present
        out.append((smi, name.strip()))
    return out


def obabel_protonate(smi, ph):
    """Canonical SMILES of the Open Babel microspecies of `smi` at pH `ph` (None on failure)."""
    r = subprocess.run(["obabel", "-:" + smi, "-osmi", "-p", str(ph)],
                       capture_output=True, text=True)
    out = r.stdout.strip().split("\t")[0].strip()
    if not out:
        return None
    m = Chem.MolFromSmiles(out)
    return None if m is None else Chem.MolToSmiles(m)


def enumerate_microstates(smi, max_states):
    """Distinct canonical microstate SMILES over (obabel protonation band) x (RDKit tautomers).

    Returns (states, meta) where states is a list of dicts sorted by canonical SMILES:
      {smiles, source} with source describing how the microstate arose (protonation pH and
      whether it is the reference tautomer or an enumerated alternate). The reference microspecies
      (obabel -p 7.4, canonical tautomer) is guaranteed present and tagged ref=True.
    """
    ref = obabel_protonate(smi, PH_REF)
    if ref is None:
        raise ValueError(f"obabel could not protonate SMILES: {smi}")

    # Distinct protonation microspecies across the physiological band.
    prot = {}  # canonical smiles -> the pH(s) it appeared at
    for ph in PH_BAND:
        s = obabel_protonate(smi, ph)
        if s is not None:
            prot.setdefault(s, []).append(ph)

    te = rdMolStandardize.TautomerEnumerator()
    te.SetMaxTautomers(MAX_TAUTOMERS)

    states = {}  # canonical smiles -> source dict (first-seen wins for provenance)
    for pstate, phs in prot.items():
        m = Chem.MolFromSmiles(pstate)
        if m is None:
            continue
        tauts = list(te.Enumerate(m))
        if not tauts:
            tauts = [m]
        for k, t in enumerate(tauts):
            cs = Chem.MolToSmiles(t)
            if cs in states:
                continue
            states[cs] = {
                "smiles": cs,
                "protonation_pH": sorted(phs),
                "is_reference_protonation": abs(min(phs) - PH_REF) < 1e-6
                or PH_REF in phs,
                "tautomer_index": k,
            }

    # The reference (obabel -p 7.4, unchanged tautomer) MUST be a member — assert it.
    if ref not in states:
        # Reference tautomer may not be RDKit's canonical enumeration output; add it verbatim.
        states[ref] = {"smiles": ref, "protonation_pH": [PH_REF],
                       "is_reference_protonation": True, "tautomer_index": -1}
    for cs, src in states.items():
        src["is_reference"] = (cs == ref)

    ordered = sorted(states.values(), key=lambda s: s["smiles"])[:max_states]
    # Guarantee the reference survives the cap.
    if not any(s["is_reference"] for s in ordered):
        ordered = ordered[:max_states - 1] + [states[ref]]
        ordered.sort(key=lambda s: s["smiles"])
    return ordered, {"reference_smiles": ref, "n_protonation_states": len(prot)}


def audit_one(smi, name, max_states):
    """Deterministic protonation/tautomer audit + enumerated microstates for one molecule."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f"unparseable SMILES for {name}: {smi}")
    states, meta = enumerate_microstates(smi, max_states)
    protomer_explicit = (len(states) == 1)
    return {
        "name": name,
        "smiles": smi,
        "reference_smiles": meta["reference_smiles"],
        "n_protonation_states": meta["n_protonation_states"],
        "n_states": len(states),
        "protomer_explicit": protomer_explicit,
        "states": states,
    }


def embed_sdf(smiles, out_path):
    """Embed one microstate SMILES to a single 3D conformer SDF (deterministic seed)."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D  # fixed so the prepared ligand is reproducible
    if AllChem.EmbedMolecule(mol, params) != 0:
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
    ap.add_argument("--ligand-dir", default="work/candidates/protomer_ligands")
    ap.add_argument("--audit-out", default="work/candidates/protomer_audit.json")
    ap.add_argument("--max-states", type=int, default=16)
    ap.add_argument("--no-embed", action="store_true",
                    help="audit only; do not write microstate SDFs")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    os.chdir(ROOT)
    mols = read_smi(args.smi)
    records = [audit_one(smi, name, args.max_states) for smi, name in mols]

    lig_dir = pathlib.Path(args.ligand_dir)
    protomer_map = {}
    if not args.no_embed:
        lig_dir.mkdir(parents=True, exist_ok=True)
        for rec in records:
            if rec["protomer_explicit"]:
                continue  # single relevant microstate — nothing to enumerate/dock
            for k, st in enumerate(rec["states"], 1):
                tag = f"{rec['name']}__st{k:02d}"
                embed_sdf(st["smiles"], lig_dir / f"{tag}.sdf")
                protomer_map[tag] = {"candidate": rec["name"], "state_index": k, **st}
        (lig_dir / "protomer_map.json").write_text(json.dumps(protomer_map, indent=2))

    audit = {"n_molecules": len(records),
             "ph_band": PH_BAND, "ph_reference": PH_REF,
             "max_states": args.max_states, "candidates": records}
    pathlib.Path(args.audit_out).write_text(json.dumps(audit, indent=2))

    if args.json:
        print(json.dumps(audit, indent=2))
        return 0

    print("=" * 92)
    print("PROTONATION / TAUTOMER ENUMERATION AUDIT (#84)")
    print(f"smiles: {args.smi}   protonation band pH {PH_BAND} (ref {PH_REF})")
    print("=" * 92)
    print("%-16s %10s %10s %8s  %s" %
          ("candidate", "protStates", "states", "explicit", "verdict"))
    print("-" * 92)
    for r in records:
        verdict = ("EXPLICIT" if r["protomer_explicit"]
                   else f"MULTI-STATE ({r['n_states']} states)")
        print("%-16s %10d %10d %8s  %s" %
              (r["name"], r["n_protonation_states"], r["n_states"],
               str(r["protomer_explicit"]), verdict))
    multi = [r["name"] for r in records if not r["protomer_explicit"]]
    print("-" * 92)
    print(f"protomer-explicit: {len(records) - len(multi)}/{len(records)}   "
          f"multi-state (microstate-docked): {multi or 'none'}")
    if not args.no_embed and protomer_map:
        print(f"wrote {len(protomer_map)} microstate SDFs -> {lig_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
