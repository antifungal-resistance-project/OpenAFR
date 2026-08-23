"""Stage 0: how big is the resistance-active validation pool if the domain ceiling is lifted?

Pipeline stage: pre-receptor audit (the cheap, local, potentially-fatal gate on the
flexible/ensemble-receptor work — docs/designs/flexible-receptor-domain-ceiling.md, Stage 0).

Why this exists
---------------
T1 concluded that a resistance-*breaker* validation is not achievable with the RIGID receptor,
because "resistance-retention" is confounded with "over the receptor's <=45-heavy-atom domain"
(data/ligands/PROVENANCE.md). A flexible/ensemble receptor could lift that size cap. But the
receptor is only worth building if lifting the cap actually rebuilds the validation pool to a
size that supports a differential-activity permutation test (the design's working bar: n>=4).

This script answers that WITHOUT any docking. It reuses the exact novelty definition the
project already uses (Morgan r=2, 2048-bit, Tanimoto < 0.70 vs BOTH the training and holdout
actives — PROVENANCE constraint #2) and the exact applicability-domain cap (<=45 heavy atoms,
scripts/filter_library.py). For each candidate resistance-retaining compound it reports:

  heavy atoms | max Tanimoto vs the 15-active panel | has C. auris evidence

and classifies it into the pool under two regimes:

  RIGID         : evidence AND novelty(<0.70) AND heavy<=45   (today's tool)
  DOMAIN-LIFTED : evidence AND novelty(<0.70)                 (if a flexible receptor is built)

The delta between the two counts is exactly what the receptor work would buy. If the
DOMAIN-LIFTED pool is still < 4, the compound universe -- not the receptor -- is the binding
constraint, and the flexible-receptor build is not justified on its own.

Deterministic and rdkit-only: runs anywhere (no vina/obabel/network), same as the gate math.

Usage:
    conda run -n openafr python scripts/domain_lifted_pool_audit.py
    conda run -n openafr python scripts/domain_lifted_pool_audit.py \\
        --candidates data/ligands/resistance_active_candidates.tsv --min-pool 4
"""
import argparse
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs

RDLogger.DisableLog("rdApp.*")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIG = os.path.join(ROOT, "data", "ligands")
TRAIN = os.path.join(LIG, "actives.smi")
HOLDOUT = os.path.join(LIG, "actives_holdout_final.smi")
CANDIDATES = os.path.join(LIG, "resistance_active_candidates.tsv")

NOVELTY_MAX = 0.70   # PROVENANCE constraint #2: Morgan r=2 Tanimoto < 0.70 vs training AND holdout
MAX_HEAVY = 45       # applicability domain, scripts/filter_library.py / PREREGISTRATION_run2.md

_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    return _GEN.GetFingerprint(mol), mol.GetNumHeavyAtoms()


def read_smi(path):
    """`SMILES<TAB>name` -> [(name, smiles)], the project's active-set format."""
    rows = []
    for line in open(path):
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        smi, _, name = line.partition("\t")
        rows.append((name.strip(), smi.strip()))
    return rows


def read_candidates(path):
    """TSV with header `name cid cauris_evidence smiles` (# comments skipped)."""
    rows = []
    for line in open(path):
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if parts[0] == "name":            # header
            continue
        if len(parts) != 4:
            sys.exit(f"ERROR: candidate row is not 4 tab-separated fields: {line!r}")
        name, cid, evidence, smiles = (p.strip() for p in parts)
        rows.append((name, cid, evidence, smiles))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates", default=CANDIDATES)
    ap.add_argument("--min-pool", type=int, default=4,
                    help="pool size the design deems sufficient for a permutation test (default 4)")
    args = ap.parse_args()

    panel = read_smi(TRAIN) + read_smi(HOLDOUT)
    panel_fps = []
    for name, smi in panel:
        f, _ = fp(smi)
        if f is None:
            sys.exit(f"ERROR: unparseable panel SMILES for {name}")
        panel_fps.append((name, f))
    print(f"reference panel : {len(panel_fps)} actives "
          f"({len(read_smi(TRAIN))} training + {len(read_smi(HOLDOUT))} holdout)")
    print(f"novelty rule    : Morgan r=2, 2048-bit, max Tanimoto < {NOVELTY_MAX}")
    print(f"domain cap      : <= {MAX_HEAVY} heavy atoms (rigid regime only)\n")

    cands = read_candidates(args.candidates)
    hdr = f"{'compound':<16} {'heavy':>5} {'maxTc':>6} {'nearest':<14} {'evid':>5} " \
          f"{'novel':>6} {'rigid':>6} {'lifted':>7}"
    print(hdr)
    print("-" * len(hdr))

    rigid_pool, lifted_pool = [], []
    for name, cid, evidence, smiles in cands:
        f, heavy = fp(smiles)
        if f is None:
            sys.exit(f"ERROR: unparseable candidate SMILES for {name} (CID {cid})")
        sims = [(DataStructs.TanimotoSimilarity(f, pf), pn) for pn, pf in panel_fps]
        max_tc, nearest = max(sims)
        has_evidence = evidence.lower() != "none"
        is_novel = max_tc < NOVELTY_MAX
        in_rigid = has_evidence and is_novel and heavy <= MAX_HEAVY
        in_lifted = has_evidence and is_novel
        if in_rigid:
            rigid_pool.append(name)
        if in_lifted:
            lifted_pool.append(name)
        print(f"{name:<16} {heavy:>5} {max_tc:>6.3f} {nearest:<14} "
              f"{'yes' if has_evidence else 'NO':>5} {'yes' if is_novel else 'NO':>6} "
              f"{'IN' if in_rigid else '-':>6} {'IN' if in_lifted else '-':>7}")

    print()
    print(f"RIGID pool         (evidence + novel + <=45 heavy) : n={len(rigid_pool)}  {rigid_pool}")
    print(f"DOMAIN-LIFTED pool (evidence + novel, no size cap) : n={len(lifted_pool)}  {lifted_pool}")
    unlocked = [c for c in lifted_pool if c not in rigid_pool]
    print(f"unlocked by lifting the domain cap                 : n={len(unlocked)}  {unlocked}")
    print()

    if len(lifted_pool) >= args.min_pool:
        print(f"VERDICT: domain-lifted pool n={len(lifted_pool)} >= {args.min_pool} — a "
              "flexible/ensemble receptor could support a differential-activity gate.")
    else:
        print(f"VERDICT: domain-lifted pool n={len(lifted_pool)} < {args.min_pool} — lifting the "
              "domain cap alone does NOT reach a testable pool.")
        print("         The binding constraint is the COMPOUND UNIVERSE (verified, novel, "
              "C.auris-active CYP51 inhibitors), not only the rigid receptor.")


if __name__ == "__main__":
    main()
