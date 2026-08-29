"""Flag candidates that ChEMBL or PubChem has already measured against the target.

Pipeline stage: report/  (an annotation layer — run it on the same library you rank)

Why this exists
---------------
rank_candidates.py orders molecules by geometry; it has no idea whether a molecule was
already assayed against fungal CYP51 and found inactive. Because failing measurements are
under-published, proposing such a molecule for the bench risks re-running a lost experiment.
This stage asks ChEMBL and PubChem BioAssay, per candidate, "is there a deposited assay
record, and what did it say?" and writes the merged answer alongside the rank so a human can
weigh it. PubChem is the larger store of HTS *inactives*; ChEMBL carries the curated,
quantitative literature. A hit in either is caught; one active in either flips the verdict.

Read the verdict honestly (openafr/precedent.py says the same, louder):
  tested-active       measured active against the CYP51 enzyme — a KNOWN inhibitor, not novel.
  tested-inactive     measured INACTIVE against the enzyme — the file-drawer catch; think hard
                      before spending bench time re-testing it.
  tested-ambiguous    measured, but weak/mixed — inconclusive precedent.
  record-inconclusive found in ChEMBL against the enzyme, but no usable quantitative verdict.
  no-precedent        no ChEMBL record found. This is NOT "untested" — a failing measurement
                      may simply never have been deposited. Absence is not evidence.
  unparseable         RDKit could not read the SMILES; cannot look it up.

The 'phenotypic' column is a whole-cell MIC against a FUNGUS and is reported separately and more
weakly than an enzyme verdict: a cell assay folds in permeability and efflux, so an inactive MIC
does not prove the enzyme was untouched — treat it as a hint to read the underlying assays, not a
fungal verdict. It is now organism-filtered (#79): a µg/mL record against a bacterium or a
mammalian cell line is NOT antifungal evidence and is counted as off-target instead, so a
`tested-active` here means a fungal cell result, not (e.g.) an antibacterial one. The
'phenotypic_organism' column names the fungal organisms behind that verdict (and, when the
fungal bucket is empty, the non-fungal/unspecified organisms whose µg/mL MICs were reclassified),
so the verdict carries the organism it rests on. Off-target ChEMBL records are shown only as a
count — a molecule already characterised against other targets is a known compound, a novelty
caveat, not a hit here.

This makes network calls to ChEMBL (www.ebi.ac.uk) and PubChem (pubchem.ncbi.nlm.nih.gov).
It fails per-molecule: a lookup error is recorded as 'lookup-failed' for that row and the run
continues.

A candidate whose exact molecule has no record can still be the close analogue of one that was
measured inactive — a soft negative. Pass --near TESTED.smi to add the nearest previously-tested
neighbour (Tanimoto ≥ 0.85) as three columns. The verified-inactives this repo builds
(scripts/make_verified_inactives.py) are the intended input: "nearest tested neighbour 0.91,
tested-inactive" tells a reader an all-but-identical compound already failed. It never changes
the verdict — it annotates it.

Usage: python scripts/check_precedent.py LIBRARY.smi [-o PRECEDENT.tsv] [--source both|chembl|pubchem] [--near TESTED.smi]
  LIBRARY.smi   tab-separated `SMILES<TAB>name`, one per line (filter_library.py format)
  -o / --out    write the annotated TSV here (default: stdout)
  --source      which databases to query (default: both)
  --near        a previously-tested set (`smiles<TAB>id[<TAB>verdict]`) for the neighbour check
"""
import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import precedent
from openafr import filedrawer

_SOURCES = {"both": ("chembl", "pubchem"), "chembl": ("chembl",), "pubchem": ("pubchem",)}


def check(lines, sources=("chembl", "pubchem")):
    """Yield (name, smiles, summary) for each `SMILES<TAB>name` line.

    summary is precedent.lookup()'s dict, or a {'verdict': 'lookup-failed', ...} stand-in when
    the network call raised — one bad lookup never aborts the batch. Malformed lines (not two
    tab-separated fields) are skipped, matching filter_library.py.
    """
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 2:
            continue
        smi, name = parts
        try:
            summary = precedent.lookup(smi, sources=sources)
        except Exception as e:                       # network/parse failure for this molecule
            summary = {"verdict": "lookup-failed", "chembl_id": None, "pubchem_cid": None,
                       "inchikey": None, "phenotypic": None, "n_off_target": 0,
                       "n_records": 0, "error": str(e)}
        yield name, smi, summary


def _pheno(summary):
    p = summary.get("phenotypic")
    return p["verdict"] if p else "—"


_ORG_CAP = 6  # keep the TSV cell legible; the full list lives in the summary dict / JSON


def _cap(names):
    """A short, legible join of an organism list — first _ORG_CAP names, then '+N more'."""
    if len(names) <= _ORG_CAP:
        return "; ".join(names)
    return "; ".join(names[:_ORG_CAP]) + f"; +{len(names) - _ORG_CAP} more"


def _pheno_organism(summary):
    """The organisms behind the phenotypic verdict (#79), for the TSV.

    Reports the fungal organisms the whole-cell verdict rests on (capped for legibility — the
    full set is in the summary dict). When no fungal MIC exists but non-fungal µg/mL records were
    reclassified out, it shows those (prefixed 'reclassified:') so a reader can see why a
    previously-flagged molecule (e.g. an antibacterial) now reads '—'.
    """
    orgs = summary.get("phenotypic_organisms") or {}
    fungal = orgs.get("fungal") or []
    if fungal:
        return f"{len(fungal)} fungal: " + _cap(fungal)
    nonfungal = orgs.get("nonfungal") or []
    if nonfungal:
        return "reclassified: " + _cap(nonfungal)
    # A phenotypic verdict with no ChEMBL organism breakdown came from PubChem, whose phenotypic
    # bucket is matched by fungal assay *name* (no structured organism) and so is already
    # fungal-gated — say so rather than print a bare '—' that reads as 'no evidence'.
    if summary.get("phenotypic"):
        return "PubChem (fungal name-matched)"
    return "—"


def _near_cols(smiles, neighbours):
    """(near_id, near_sim, near_verdict) for one candidate, or ('—','—','—') when none close."""
    if neighbours is None:
        return None
    fp = precedent.fingerprint(smiles)
    hit = precedent.nearest_tested_neighbour(fp, neighbours) if fp is not None else None
    if hit is None:
        return ("—", "—", "—")
    return (hit["id"], f"{hit['similarity']:.2f}", hit["verdict"])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("library", help="tab-separated SMILES<TAB>name library")
    ap.add_argument("-o", "--out", help="write the annotated TSV here (default: stdout)")
    ap.add_argument("--source", choices=sorted(_SOURCES), default="both",
                    help="which databases to query (default: both)")
    ap.add_argument("--near", help="previously-tested set for the nearest-neighbour column")
    args = ap.parse_args()

    neighbours = precedent.load_neighbours(args.near) if args.near else None

    header = ("name\tverdict\tchembl_id\tpubchem_cid\tphenotypic\tphenotypic_organism\t"
              "n_off_target\tn_records")
    if neighbours is not None:
        header += "\tnear_id\tnear_sim\tnear_verdict"
    rows = [header]
    counts = {}
    for name, smi, s in check(open(args.library), _SOURCES[args.source]):
        counts[s["verdict"]] = counts.get(s["verdict"], 0) + 1
        row = (f"{name}\t{s['verdict']}\t{s.get('chembl_id') or '—'}\t"
               f"{s.get('pubchem_cid') or '—'}\t{_pheno(s)}\t{_pheno_organism(s)}\t"
               f"{s.get('n_off_target', 0)}\t{s.get('n_records', 0)}")
        near = _near_cols(smi, neighbours)
        if near is not None:
            row += "\t" + "\t".join(near)
        rows.append(row)
    text = "\n".join(rows) + "\n"

    if args.out:
        pathlib.Path(args.out).write_text(text)
        print(f"# wrote {args.out}", file=sys.stderr)
    else:
        sys.stdout.write(text)

    print("\n# precedent tally: "
          + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())), file=sys.stderr)
    print("# 'no-precedent' means no enzyme record was found in the queried source(s) — NOT that "
          "the molecule was never tested; a failing measurement may simply never have been "
          "deposited.", file=sys.stderr)
    print("# " + filedrawer.caveat_line(), file=sys.stderr)
    print("# 'tested-inactive' is the one to heed: a wet lab already measured it inactive against "
          "the CYP51 enzyme.", file=sys.stderr)


if __name__ == "__main__":
    main()
