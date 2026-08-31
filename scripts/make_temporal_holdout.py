"""Build the issue-#82 TEMPORAL holdout — a publication-date split that holds out the
NEWEST-literature azoles and tests whether the frozen geometric criterion, defined on the
classic approved azoles, ranks recent research chemistry it never saw.

Pipeline stage: ligands/ (temporal holdout)

The stronger follow-up #82 asks for. The disjoint-population holdout already landed (#108,
scripts/make_external_holdout.py + validate_gate_external.py): a never-*scored* subset of the
curated actives. Its one honest residual — "a clean never-touched holdout, but not a later
time slice" — is exactly what this closes. Here the split is TEMPORAL:

  the criterion (the N-Fe < 3.0 A iron filter, its cutoff, the frozen run-2 protocol) was
  DEFINED on the classic approved azoles (fluconazole, itraconazole, ... — antifungal
  literature reaching back to the 1980s-2000s). This holdout is the azoles whose EARLIEST
  antifungal / CYP51 literature appearance is recent (document year > CUTOFF_YEAR = 2023, i.e.
  >= 2024) and that no gate has ever docked. So it is both never-touched AND a genuine
  later-time slice of the chemistry — new research compounds, not old drugs.

Scope, stated honestly (see work/PREREGISTRATION_temporal.md): this is a publication-year split
of the CURRENT ChEMBL release (ChEMBL_37, 2026-05-01 — the release the corpus was pulled from).
It is NOT a *future*-prospective set post-dating the release; ChEMBL_37 is the latest release,
so that stronger variant is calendar-gated for everyone until ChEMBL_38 and re-runs this exact
protocol when it drops. What this DOES deliver now: a real temporal generalization test the
criterion never saw, executed end-to-end.

Discipline (mirrors data/ligands/PROVENANCE.md — only the provenance of "held out" changes, to
"earliest antifungal literature is recent"):

  actives  a deterministic seed-82 sample of the verified_actives azoles whose earliest
           antifungal / CYP51 document year > 2023, that are NOT in the seed-42 docked sample
           and NOT in the #108 external holdout (so never docked by any gate), in the frozen
           domain (<= 45 heavy, >= 1 aromatic N; verified_actives already enforces both).

  decoys   property-matched presumed-inactives, built by the SAME frozen rule as every other
           decoy set (openafr.inactives.select_inactives), fetched fresh from ChEMBL and
           disjoint by ChEMBL id from every used set + the temporal actives. Decoys are NOT
           date-constrained: the temporal variable is the ACTIVES' publication era; the single
           changed variable vs #108 is which actives are held out.

The per-molecule earliest antifungal year comes from data/ligands/chembl_document_years.json
(built once by `years`, pinned for offline reproducibility). Deterministic: same year map +
same verified_actives -> byte-identical temporal_actives.smi.

Usage:
    python scripts/make_temporal_holdout.py years              # (re)fetch the year map (network)
    python scripts/make_temporal_holdout.py actives            # offline, deterministic
    python scripts/make_temporal_holdout.py decoys             # fetches ChEMBL
"""
import argparse
import json
import pathlib
import random
import sys
import time
import urllib.parse
import urllib.request

from rdkit import Chem, RDLogger

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import inactives
from scripts.make_external_holdout import (USED_ID_FILES, _fetch_pool, _id_of, _read_smi, _sha)

RDLogger.DisableLog("rdApp.*")

# The temporal boundary (FROZEN — see work/PREREGISTRATION_temporal.md). The criterion was
# defined on the classic approved azoles; the holdout is azoles whose earliest antifungal
# literature is recent. > 2023 (i.e. >= 2024) gives a well-powered, never-docked pool (320
# candidates) that is a clear later-time slice from the old approved drugs.
CUTOFF_YEAR = 2023

VERIFIED = "data/ligands/verified_actives.smi"
SCORED = "data/ligands/verified_actives_sample.smi"           # the seed-42 n=200 already docked
EXTERNAL = "data/ligands/external_actives.smi"                 # the #108 holdout (also never re-use)
YEARS = "data/ligands/chembl_document_years.json"             # {molecule_id: earliest antifungal year}
ACTIVES_OUT = "data/ligands/temporal_actives.smi"
DECOYS_OUT = "data/ligands/temporal_decoys.smi"
SEED = 82                                                      # #82's provenance seed; NOT 42 (tuning)
N_ACTIVES = 50                                                 # match #108 for a like-for-like contrast
PER_ACTIVE = 8

# every set the criterion has ever seen, INCLUDING the #108 external sets — decoys must be
# disjoint from all of them by id.
USED_FILES = list(USED_ID_FILES) + [EXTERNAL, "data/ligands/external_decoys.smi"]

# the two antifungal targets the verified sets were built from (whole-cell C. albicans + the
# C. albicans CYP51 enzyme), for the year fetch.
TARGET_QUERIES = ({"target_organism": "Candida albicans"}, {"target_chembl_id": "CHEMBL1780"})
ACTIVITY = "https://www.ebi.ac.uk/chembl/api/data/activity.json"


# --- pure, offline-testable temporal logic -------------------------------------------------

def is_temporal(earliest_year, cutoff=CUTOFF_YEAR):
    """A molecule is in the temporal holdout iff its earliest antifungal literature appearance
    is strictly AFTER the cutoff era. An unknown year (None) is NOT temporal — absence of proof
    that it is recent is not proof."""
    return earliest_year is not None and earliest_year > cutoff


def load_years(path=YEARS):
    """{molecule_id: earliest antifungal / CYP51 document year} from the pinned map."""
    return json.load(open(path))


def select_temporal_actives(pool, years, scored_ids, external_ids, cutoff=CUTOFF_YEAR,
                            n=N_ACTIVES, seed=SEED):
    """The deterministic core (unit-tested offline). From verified_actives, keep the never-docked
    molecules whose earliest antifungal year > cutoff, then a fixed-seed sample of n. Sorted by
    id first so the sample is independent of file order; same inputs -> same output."""
    cand = [(smi, name) for smi, name in pool
            if _id_of(name) not in scored_ids and _id_of(name) not in external_ids
            and is_temporal(years.get(_id_of(name)), cutoff)]
    cand.sort(key=lambda r: int(_id_of(r[1])[6:]) if _id_of(r[1])[6:].isdigit() else 10**12)
    if n >= len(cand):
        return cand
    rng = random.Random(seed)
    return [cand[i] for i in sorted(rng.sample(range(len(cand)), n))]


def used_ids(files=USED_FILES):
    used = set()
    for f in files:
        if pathlib.Path(f).exists():
            used |= {_id_of(n) for _s, n in _read_smi(f)}
    return used


# --- network shell -------------------------------------------------------------------------

def _get(url, retries=4):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except Exception as e:
            if attempt == retries - 1:
                raise SystemExit(f"ChEMBL fetch failed after {retries} tries: {e}\n  {url}")
            time.sleep(2 ** attempt)


def cmd_years(args):
    """Fetch the earliest antifungal / CYP51 document year for every C. albicans compound and
    pin it. Deterministic given a ChEMBL release; re-run only to refresh against a new release."""
    minyear = {}
    for sel in TARGET_QUERIES:
        params = dict(sel, only="molecule_chembl_id,document_year", limit=1000, offset=0, format="json")
        total = _get(f"{ACTIVITY}?{urllib.parse.urlencode(dict(params, limit=1))}")["page_meta"]["total_count"]
        n = 0
        while n < total:
            page = _get(f"{ACTIVITY}?{urllib.parse.urlencode(dict(params, offset=n))}")
            rows = page.get("activities", [])
            if not rows:
                break
            for rec in rows:
                cid, yr = rec.get("molecule_chembl_id"), rec.get("document_year")
                if cid and isinstance(yr, int) and (cid not in minyear or yr < minyear[cid]):
                    minyear[cid] = yr
            n += len(rows)
            print(f"  {list(sel.values())[0]}: {n}/{total}", end="\r", flush=True)
        print()
    json.dump(minyear, open(YEARS, "w"), sort_keys=True)
    print(f"wrote earliest-year map for {len(minyear)} compounds to {YEARS}")
    print(f"sha256 {_sha(YEARS)}")
    return 0


def cmd_actives(args):
    pool = _read_smi(VERIFIED)
    years = load_years()
    scored = {_id_of(n) for _s, n in _read_smi(SCORED)}
    external = {_id_of(n) for _s, n in _read_smi(EXTERNAL)}
    n_recent = sum(1 for _s, name in pool
                   if is_temporal(years.get(_id_of(name))) and _id_of(name) not in scored
                   and _id_of(name) not in external)
    print(f"verified actives              : {len(pool)}")
    print(f"never-docked, earliest yr > {CUTOFF_YEAR} : {n_recent}")

    chosen = select_temporal_actives(pool, years, scored, external)
    for smi, name in chosen:                                   # re-assert the frozen domain
        mol = Chem.MolFromSmiles(smi)
        if mol is None or mol.GetNumHeavyAtoms() > inactives.MAX_HEAVY \
                or not inactives.has_aromatic_nitrogen(mol):
            sys.exit(f"{name} out of domain — should be impossible for a verified active")
    with open(ACTIVES_OUT, "w") as fh:
        for smi, name in chosen:
            fh.write(f"{smi}\t{name}\n")
    yrs = [years[_id_of(n)] for _s, n in chosen]
    heavy = [Chem.MolFromSmiles(s).GetNumHeavyAtoms() for s, _ in chosen]
    print(f"\nwrote {len(chosen)} temporal actives (seed {SEED}) to {ACTIVES_OUT}")
    print(f"earliest antifungal year: min {min(yrs)}  max {max(yrs)}   heavy: {min(heavy)}-{max(heavy)}")
    print(f"sha256 {_sha(ACTIVES_OUT)}")
    return 0


def cmd_decoys(args):
    actives = inactives.load_actives(ACTIVES_OUT)
    active_fps = [a["fp"] for a in actives]
    mws = [a["mw"] for a in actives]
    print(f"temporal actives              : {len(actives)}  (MW {min(mws):.0f}-{max(mws):.0f})")
    used = used_ids() | {_id_of(a["name"]) for a in actives}
    print(f"excluded ids (all sets)       : {len(used)}")

    band_lo, band_hi = min(mws) - inactives.MW_TOL, max(mws) + inactives.MW_TOL
    cand = {}
    for offset in range(0, args.pages * 1000, 1000):
        pool = _fetch_pool(band_lo, band_hi, offset)
        time.sleep(0.2)
        for cid, smi in pool:
            if cid not in used and cid not in cand:
                cand[cid] = smi
        print(f"    offset {offset:5d}: pool {len(pool):4d}  distinct usable {len(cand)}")
        if not pool:
            break
    candidates = sorted(cand.items(), key=lambda kv: int(kv[0][6:]) if kv[0][6:].isdigit() else 10**12)
    print(f"candidate pool                : {len(candidates)} distinct, disjoint by id")

    chosen, rejects = inactives.select_inactives(candidates, actives, active_fps, PER_ACTIVE)
    print("\nfilter funnel:")
    for k, v in sorted(rejects.items(), key=lambda kv: -kv[1]):
        print(f"  rejected {k:20} {v:6d}")
    print(f"  -> KEPT {len(chosen)} property-matched decoys")
    chosen.sort(key=lambda c: int(c["chembl_id"][6:]) if c["chembl_id"][6:].isdigit() else 10**12)
    with open(DECOYS_OUT, "w") as fh:
        for c in chosen:
            fh.write(f"{c['smiles']}\tdecoy_{c['chembl_id']}_for_{c['paired']}\n")
    print(f"\nwrote {len(chosen)} temporal decoys to {DECOYS_OUT}")
    print(f"sha256 {_sha(DECOYS_OUT)}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("years", help="fetch + pin the earliest antifungal document year map (network)"
                   ).set_defaults(func=cmd_years)
    sub.add_parser("actives", help="offline: deterministic newest-literature never-docked sample"
                   ).set_defaults(func=cmd_actives)
    d = sub.add_parser("decoys", help="fetch ChEMBL + build property-matched disjoint decoys")
    d.add_argument("--pages", type=int, default=12, help="1000-molecule ChEMBL pages to scan")
    d.set_defaults(func=cmd_decoys)
    args = ap.parse_args(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    sys.exit(main())
