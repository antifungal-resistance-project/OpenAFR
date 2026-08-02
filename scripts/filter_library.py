"""Filter a novel screening library down to the validated applicability domain.

Pipeline stage: filter/  (runs on a novel SMILES library, before prep_ligands.py)

Why this stage exists
---------------------
The gate validates ONE ranking criterion — closest approach of a ligand nitrogen
to the heme iron — under a rigid receptor and a frozen protocol. That validation
does not extend to every molecule; the pre-registration declares, in advance, the
envelope the result is claimed to hold over. Ranking molecules outside that
envelope would present numbers the pipeline has no validated basis for. This stage
enforces the envelope in code, before any docking compute is spent, so the ranked
list only ever contains molecules the gate's evidence actually covers.

The two exclusion rules, both declared — not tuned
--------------------------------------------------
  applicability domain : <= 45 heavy atoms
      Pre-registered in work/PREREGISTRATION_run2.md, section "Applicability
      domain". Larger ligands (posaconazole/itraconazole class) failed reproducibly
      in run 1, consistent with induced fit a rigid receptor cannot model. This is
      the project's declared blind spot, not a knob turned to improve a score.

  mechanistic scope    : >= 1 aromatic nitrogen
      The validated criterion is a nitrogen-to-iron distance. A molecule with no
      aromatic nitrogen cannot coordinate the heme iron and is undefined under that
      criterion — it is outside the mechanism the gate tested, not merely a poor
      binder. (rank_candidates.py would already sink it to last for want of a pose;
      excluding it here just avoids docking a molecule the criterion cannot score.)

Deliberately NOT applied: drug-likeness filters (Lipinski, PAINS, etc.). None were
pre-registered, and adding un-validated criteria to flatter the candidate list is
exactly the anti-tuning trap the rest of this pipeline is built to avoid.

Nothing is silently dropped. Every excluded or unparseable molecule is written to
the excluded log with a reason and counted, in keeping with prep_ligands.py.

Usage: python scripts/filter_library.py LIBRARY.smi [-o KEEP.smi] [-x EXCLUDED.tsv]
  LIBRARY.smi   tab-separated `SMILES<TAB>name`, one per line (prep_ligands format)
  -o / --keep   in-domain molecules, ready for prep_ligands.py (default: stdout)
  -x / --excluded  excluded molecules with reasons        (default: stderr summary only)
"""
import argparse
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

MAX_HEAVY_ATOMS = 45                     # applicability domain (PREREGISTRATION_run2.md)
AROMATIC_N = Chem.MolFromSmarts("[n]")   # same iron-coordinator probe as make_decoys.py


def classify(smiles):
    """Return (keep: bool, reason: str) for one SMILES string.

    reason is "" when kept, else a short human-readable exclusion cause. The order
    of checks is fixed so a molecule failing several rules reports one stable cause.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "unparseable SMILES"
    heavy = mol.GetNumHeavyAtoms()
    if heavy > MAX_HEAVY_ATOMS:
        return False, f"{heavy} heavy atoms > {MAX_HEAVY_ATOMS} (outside applicability domain)"
    if not mol.HasSubstructMatch(AROMATIC_N):
        return False, "no aromatic nitrogen (cannot coordinate the heme iron)"
    return True, ""


def filter_library(lines):
    """Partition `SMILES<TAB>name` lines into (kept, excluded).

    kept     : [(smiles, name), ...] preserving input order
    excluded : [(name, reason), ...] preserving input order
    Malformed lines (not exactly two tab-separated fields) are skipped, matching
    prep_ligands.py, which reads the same format.
    """
    kept, excluded = [], []
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 2:
            continue
        smi, name = parts
        keep, reason = classify(smi)
        if keep:
            kept.append((smi, name))
        else:
            excluded.append((name, reason))
    return kept, excluded


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("library", help="tab-separated SMILES<TAB>name library")
    ap.add_argument("-o", "--keep", help="write in-domain SMILES here (default: stdout)")
    ap.add_argument("-x", "--excluded", help="write excluded molecules + reasons here (TSV)")
    args = ap.parse_args()

    with open(args.library) as fh:
        kept, excluded = filter_library(fh)

    keep_text = "".join(f"{smi}\t{name}\n" for smi, name in kept)
    if args.keep:
        with open(args.keep, "w") as fh:
            fh.write(keep_text)
    else:
        sys.stdout.write(keep_text)

    if args.excluded:
        with open(args.excluded, "w") as fh:
            fh.write("name\treason\n")
            for name, reason in excluded:
                fh.write(f"{name}\t{reason}\n")

    total = len(kept) + len(excluded)
    print(f"\nfiltered: {total} in -> {len(kept)} in-domain, {len(excluded)} excluded",
          file=sys.stderr)
    for name, reason in excluded:
        print(f"  EXCLUDED {name}: {reason}", file=sys.stderr)
    if args.keep:
        print(f"# in-domain library -> {args.keep} (feed to prep_ligands.py)", file=sys.stderr)


if __name__ == "__main__":
    main()
