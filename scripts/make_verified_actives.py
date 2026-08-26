"""Build an INDEPENDENT measured-active azole set — the mirror of verified_inactives.

Pipeline stage: ligands/ (active half, independent confirmation)

The ensemble look (work/RESULTS_ensemble.md) passed at AUC 0.750 on the 7 held-out azoles,
but with a wide CI (n=7) — so the ensemble pre-registration requires a "separately
pre-registered confirmation on an INDEPENDENT active set" before any within-class claim. This
builds that set with the SAME provenance discipline the verified-inactive decoys were built
under (openafr/inactives.py): from the ChEMBL C. albicans bioactivity cache, by a frozen rule,
with ChEMBL's own canonical SMILES (never AI-generated), so no molecule is hand-picked to make
the gate pass.

The rule (single-variable mirror of the inactive construction — only the sign of the evidence
flips, everything else is held fixed):

    keep a compound iff  ALL of:
      * verdict == has-active-evidence   (>= 1 active record; openafr.inactives.compound_verdict,
        whose asymmetry rule — one active outranks any inactives — is CONSERVATIVE for an active
        benchmark too: a false active ranks low and only DEPRESSES measured separation)
      * azole-bearing                    (same AZOLE_SMARTS the gate's subgroup split uses, so
        this is a within-class test: measured-active azoles vs measured-inactive azoles)
      * in the applicability domain      (<= 45 heavy atoms, the frozen run-2 domain)
      * >= 1 aromatic nitrogen           (the mechanism requirement, as for the decoys)
      * NOVEL vs every molecule already used  (max Morgan r=2 Tanimoto < 0.70 against the 8
        training + 7 held-out actives AND oteseconazole/VT-1161, the 5TZ1 co-crystal ligand —
        so the set is genuinely distinct chemistry, exact-matches of the 15 are excluded, and
        the 5TZ1 co-crystal cannot leak in as a self-docking active)

Deterministic: same cache -> same .smi (compounds emitted in ascending ChEMBL-id order).

Usage:
    python scripts/make_verified_actives.py build [--cache DIR] [--out FILE]
"""
import argparse
import collections
import hashlib
import pathlib
import sys

from rdkit import Chem, DataStructs

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import inactives
from scripts.make_verified_inactives import ALL_ACTIVES, DEFAULT_CACHE, load_cache
from scripts.validate_gate_verified import AZOLE_SMARTS

DEFAULT_OUT = "data/ligands/verified_actives.smi"
NOVELTY_MAX = 0.70            # discipline constraint #2 (data/ligands/PROVENANCE.md)
MAX_HEAVY = inactives.MAX_HEAVY
# The 5TZ1 co-crystal azole (oteseconazole / VT-1161), CID 77050711, SMILES from PubChem
# (data/ligands/PROVENANCE.md). Added to the novelty reference so the reference ligand — and
# anything ~identical to it — cannot enter the active set and leak through 5TZ1 self-docking.
OTESECONAZOLE = ("C1=CC(=CC=C1C2=CN=C(C=C2)C([C@](CN3C=NN=N3)"
                 "(C4=C(C=C(C=C4)F)F)O)(F)F)OCC(F)(F)F")


def is_azole(mol):
    return any(mol.HasSubstructMatch(q) for q in AZOLE_SMARTS)


def cmd_build(args):
    by_cmp, smiles, counts = load_cache(args.cache)
    print("cache: " + "  ".join(f"{k}={v} records" for k, v in counts.items())
          + f"  ({len(by_cmp)} distinct compounds)")

    # Novelty reference: every active the project has ever used + the 5TZ1 co-crystal ligand.
    ref_fps = []
    for path in ALL_ACTIVES:
        ref_fps += [a["fp"] for a in inactives.load_actives(path)]
    ref_fps.append(inactives.fingerprint(Chem.MolFromSmiles(OTESECONAZOLE)))

    def sort_key(cid):
        return int(cid.replace("CHEMBL", "")) if cid[6:].isdigit() else 10**12

    verdicts = collections.Counter()
    rejects = collections.Counter()
    chosen, seen = [], set()
    for cid in sorted(by_cmp, key=sort_key):
        verdict, _tally = inactives.compound_verdict(by_cmp[cid])
        verdicts[verdict] += 1
        if verdict != "has-active-evidence" or cid not in smiles:
            continue
        smi = smiles[cid]
        if smi in seen:
            rejects["duplicate"] += 1
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            rejects["unparseable"] += 1
            continue
        if not is_azole(mol):
            rejects["not-azole"] += 1
            continue
        if not inactives.has_aromatic_nitrogen(mol):
            rejects["no-aromatic-n"] += 1
            continue
        if mol.GetNumHeavyAtoms() > MAX_HEAVY:
            rejects["over-domain"] += 1
            continue
        fp = inactives.fingerprint(mol)
        if max(DataStructs.BulkTanimotoSimilarity(fp, ref_fps)) >= NOVELTY_MAX:
            rejects["too-similar-to-existing"] += 1
            continue
        seen.add(smi)
        chosen.append((cid, smi, mol.GetNumHeavyAtoms()))

    print("\nevidence verdicts over every compound ever tested on C. albicans:")
    for k, v in verdicts.most_common():
        print(f"  {k:22} {v:6d}")
    print(f"\nfilter funnel (has-active-evidence compounds with a structure):")
    for k, v in rejects.most_common():
        print(f"  rejected {k:24} {v:6d}")
    print(f"  -> KEPT {len(chosen)} independent measured-active azoles")

    if chosen:
        import statistics as st
        print(f"\n  heavy atoms: min {min(c[2] for c in chosen)}  "
              f"mean {st.mean(c[2] for c in chosen):.1f}  max {max(c[2] for c in chosen)}")

    with open(args.out, "w") as fh:
        for cid, smi, _h in chosen:
            fh.write(f"{smi}\tactive_{cid}\n")
    digest = hashlib.sha256(open(args.out, "rb").read()).hexdigest()
    print(f"\nwrote {len(chosen)} independent actives to {args.out}")
    print(f"sha256 {digest}")
    return 0


def cmd_sample(args):
    """Deterministic random sample of the built population -> the frozen confirmation test set.

    The full population (thousands of measured-active azoles) is too large to dock; a
    seeded random sample is unbiased (unlike a head-of-list slice) and reproducible. The
    population is sorted first so the shuffle is independent of file order, then a fixed-seed
    shuffle picks N. Same population + same seed + same N -> byte-identical sample.
    """
    import random
    pop = sorted(l for l in open(args.infile) if "\t" in l)
    if args.n > len(pop):
        sys.exit(f"--n {args.n} exceeds population {len(pop)}")
    rng = random.Random(args.seed)
    idx = sorted(rng.sample(range(len(pop)), args.n))
    with open(args.out, "w") as fh:
        fh.writelines(pop[i] for i in idx)
    digest = hashlib.sha256(open(args.out, "rb").read()).hexdigest()
    print(f"sampled {args.n} of {len(pop)} (seed {args.seed}) -> {args.out}")
    print(f"sha256 {digest}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="offline: classify + filter + write the population .smi")
    b.add_argument("--cache", default=DEFAULT_CACHE)
    b.add_argument("--out", default=DEFAULT_OUT)
    b.set_defaults(func=cmd_build)
    s = sub.add_parser("sample", help="deterministic random sample -> frozen test set")
    s.add_argument("--infile", default=DEFAULT_OUT)
    s.add_argument("--out", default="data/ligands/verified_actives_sample.smi")
    s.add_argument("--n", type=int, default=200)
    s.add_argument("--seed", type=int, default=42)
    s.set_defaults(func=cmd_sample)
    args = ap.parse_args(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    sys.exit(main())
