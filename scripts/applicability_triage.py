"""Per-molecule applicability triage: is the validated CYP51 method even claimed to cover THIS molecule?

Pipeline stage: the outward-facing front door. Given a candidate SMILES (a collaborator's
compound, or a whole library), answer — before any docking compute is spent, and without a
cloud VM — the one question the preprint's honest-envelope message implies:

    Does the geometry-over-docking criterion, as validated, apply to this molecule, and is it
    genuinely new chemistry vs the compounds the validation was built on?

This is NOT an activity prediction and NOT a docking score. The N-to-Fe geometric score itself
needs a docked pose (vina, gated on the Linux/x86 VM). This tool is the cheap triage that runs
first: it reports whether a molecule sits inside the pre-registered envelope the result is
claimed over, and how novel it is relative to the 15 validated actives — so nobody wastes a
docking run (or a claim) on a molecule the method has no validated basis for.

The three checks, all reused from already-validated code / pre-registrations
--------------------------------------------------------------------------
  mechanistic scope   >= 1 aromatic nitrogen  (scripts.filter_library.AROMATIC_N)
      The criterion is a nitrogen-to-iron distance; a molecule with no aromatic N cannot
      coordinate the heme iron and is undefined under it — out of scope, not merely a poor hit.
  applicability domain  <= 45 heavy atoms      (scripts.filter_library.MAX_HEAVY_ATOMS)
      Pre-registered blind spot (PREREGISTRATION_run2.md): larger ligands failed reproducibly
      in run 1 from induced fit a rigid receptor cannot model. NOTE this is honest to a fault —
      posaconazole/itraconazole (the defining azoles) are themselves out-of-domain here.
  novelty vs the validated panel  Morgan r=2 / 2048-bit, max Tanimoto < 0.70
      PROVENANCE.md constraint #2. Distinguishes genuinely new chemistry from a near-duplicate
      re-discovery of a compound the criterion was already built on.

Verdict (one stable cause, fixed priority order — same discipline as filter_library):
    unparseable       -> bad SMILES
    out-of-scope      -> no aromatic nitrogen; criterion undefined
    out-of-domain     -> > 45 heavy atoms; outside the validated envelope
    near-known        -> in envelope, but Tanimoto >= 0.70 to a validated active (rediscovery)
    in-envelope-novel -> in envelope AND novel; the method's evidence covers it

Nothing is silently dropped: every molecule gets a verdict line.

Deterministic and rdkit-only (no vina/obabel/network) — runs anywhere.

Usage:
    conda run -n openafr python scripts/applicability_triage.py --smiles "<SMILES>" [--name X]
    conda run -n openafr python scripts/applicability_triage.py LIBRARY.smi   # SMILES<TAB>name
"""
import argparse
import os
import sys

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)  # make the repo root importable when run as `python scripts/...`

# Single source of truth for the two envelope rules — do not re-hardcode them here.
from scripts.filter_library import AROMATIC_N, MAX_HEAVY_ATOMS

RDLogger.DisableLog("rdApp.*")
LIG = os.path.join(ROOT, "data", "ligands")
PANEL_FILES = [os.path.join(LIG, "actives.smi"),
               os.path.join(LIG, "actives_holdout_final.smi")]

NOVELTY_MAX = 0.70   # PROVENANCE.md constraint #2 (Morgan r=2, 2048-bit)
_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

# Human-readable, decision-first verdict text.
VERDICTS = {
    "unparseable":       "unparseable SMILES",
    "out-of-scope":      "no aromatic nitrogen — cannot coordinate the heme iron; criterion undefined",
    "out-of-domain":     f"> {MAX_HEAVY_ATOMS} heavy atoms — outside the validated applicability domain",
    "near-known":        "in envelope, but a near-duplicate of a validated active (Tanimoto >= 0.70)",
    "in-envelope-novel": "in the validated envelope AND novel chemistry — the method's evidence covers it",
}


def load_panel():
    """The 15 validated actives (training + holdout) as [(name, fingerprint)]."""
    panel = []
    for path in PANEL_FILES:
        for line in open(path):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            smi, _, name = line.partition("\t")
            mol = Chem.MolFromSmiles(smi.strip())
            if mol is None:
                sys.exit(f"ERROR: unparseable panel SMILES for {name.strip()} in {path}")
            panel.append((name.strip(), _GEN.GetFingerprint(mol)))
    return panel


def classify(smiles, panel):
    """Return a dict verdict for one SMILES. Fixed-priority, one stable cause."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"verdict": "unparseable", "heavy": None, "max_tc": None, "nearest": None}
    heavy = mol.GetNumHeavyAtoms()
    max_tc, nearest = None, None
    if panel:
        fp = _GEN.GetFingerprint(mol)
        sims = [(DataStructs.TanimotoSimilarity(fp, pf), pn) for pn, pf in panel]
        max_tc, nearest = max(sims)

    if not mol.HasSubstructMatch(AROMATIC_N):
        verdict = "out-of-scope"
    elif heavy > MAX_HEAVY_ATOMS:
        verdict = "out-of-domain"
    elif max_tc is not None and max_tc >= NOVELTY_MAX:
        verdict = "near-known"
    else:
        verdict = "in-envelope-novel"
    return {"verdict": verdict, "heavy": heavy, "max_tc": max_tc, "nearest": nearest}


def read_library(path):
    """`SMILES<TAB>name` lines -> [(name, smiles)] (# comments + blanks skipped)."""
    rows = []
    for line in open(path):
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        smi, _, name = line.partition("\t")
        rows.append((name.strip() or "(unnamed)", smi.strip()))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("library", nargs="?", help="SMILES<TAB>name library file")
    ap.add_argument("--smiles", help="a single SMILES to triage")
    ap.add_argument("--name", default="query", help="name for --smiles (default: query)")
    args = ap.parse_args()

    if bool(args.library) == bool(args.smiles):
        ap.error("give exactly one of: a LIBRARY file, or --smiles")

    rows = [(args.name, args.smiles)] if args.smiles else read_library(args.library)
    panel = load_panel()
    print(f"validated panel : {len(panel)} actives  |  domain <= {MAX_HEAVY_ATOMS} heavy  |  "
          f"novelty Tanimoto < {NOVELTY_MAX}")
    print("NOTE: applicability + novelty triage only — NOT an activity or docking prediction.\n")

    hdr = f"{'molecule':<20} {'heavy':>5} {'maxTc':>6} {'nearest':<14} {'verdict':<18} reason"
    print(hdr)
    print("-" * (len(hdr) + 30))
    counts = {}
    for name, smi in rows:
        r = classify(smi, panel)
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
        heavy = "-" if r["heavy"] is None else str(r["heavy"])
        tc = "-" if r["max_tc"] is None else f"{r['max_tc']:.3f}"
        near = r["nearest"] or "-"
        print(f"{name:<20} {heavy:>5} {tc:>6} {near:<14} {r['verdict']:<18} {VERDICTS[r['verdict']]}")

    print(f"\nsummary ({len(rows)} molecule{'s' if len(rows) != 1 else ''}):")
    for v in ("in-envelope-novel", "near-known", "out-of-domain", "out-of-scope", "unparseable"):
        if counts.get(v):
            print(f"  {counts[v]:>4}  {v}")


if __name__ == "__main__":
    main()
