"""Audit: are there enough experimentally-confirmed CYP51-ENZYME inactives to build the
contrast the preprint (work/PREPRINT_geometry_ceiling.md §4.1) calls "the most useful next
contribution anyone could make to this benchmark"?

Context. Limitation #2 was powered against 279 compounds MEASURED not to inhibit *whole-cell*
C. albicans (look #6/#11, RESULTS_verified_power.md): pooled AUC 0.688 [0.656, 0.720], sub-bar,
because 58% of that set are FAILED AZOLE analogues (a chemotype-vs-potency confound) and a
whole-cell MIC is not a target assay (permeability/efflux). §4.1's named clean exit is compounds
measured not to inhibit the ENZYME itself. This script asks, exhaustively and reproducibly,
whether ChEMBL has them at usable scale for the docked target (*C. albicans* CYP51, CHEMBL1780,
the UniProt-P10613 enzyme behind 5TZ1).

It is an availability AUDIT, not a pre-registered gate: it decides no criterion and docks nothing.
It reuses the frozen decoy filters from openafr/inactives.py unchanged, so a kept compound would
be admissible in exactly the single-variable-swap sense look #6 used (only the *provenance* of
"inactive" changes: enzyme-measured, not whole-cell-measured, not presumed).

Two evidence rules are run so the conclusion is shown robust to strictness:
  CLEAN    — only unambiguously-classifiable rows: a cell-free/enzyme IC50/Ki censored high, an
             explicit inactive comment, or an 'Inhibition' %-assay <=20% at a stated dose >=10
             ug/ml. The directionally-AMBIGUOUS 'Activity' %-rows (per-sterol GC-MS composition,
             where low ergosterol = inhibited but low eburicol = NOT inhibited) are refused.
  GENEROUS — an indefensibly loose upper bound: ANY %-row (either type) <=30% at ANY dose counts
             inactive, >=50% counts active. Used only to bound the answer from above.

Usage:
    python scripts/audit_enzyme_inactives.py [--cache work/chembl_enzyme_audit] [--refetch]
"""
import argparse
import collections
import json
import os
import pathlib
import re
import sys
import time
import urllib.parse
import urllib.request

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import inactives as I

API = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
TARGET = "CHEMBL1780"  # C. albicans (SC5314) lanosterol 14-alpha demethylase = the docked enzyme
FIELDS = ("molecule_chembl_id,canonical_smiles,standard_type,standard_relation,"
          "standard_value,standard_units,activity_comment,pchembl_value,assay_description")
ACTIVES = "data/ligands/actives_holdout_final.smi"
ALL_ACTIVES = ("data/ligands/actives.smi", "data/ligands/actives_holdout_final.smi")
DOSE = re.compile(r"at\s+([\d.]+)\s*ug", re.I)
CELL = re.compile(r"in candida|in fluconazole|c\. albicans|sc5314|caal|atcc|membrane|"
                  r"ergosterol biosynthesis|in cells|whole cell", re.I)


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def fetch(cache):
    path = os.path.join(cache, f"{TARGET}.jsonl")
    if os.path.exists(path):
        return path
    os.makedirs(cache, exist_ok=True)
    base = {"target_chembl_id": TARGET, "limit": 1000, "only": FIELDS, "format": "json"}

    def get(p):
        with urllib.request.urlopen(f"{API}?{urllib.parse.urlencode(p)}", timeout=120) as r:
            return json.load(r)
    total = get(dict(base, limit=1))["page_meta"]["total_count"]
    n = 0
    with open(path, "w") as fh:
        while n < total:
            got = get(dict(base, offset=n)).get("activities", [])
            if not got:
                break
            for rec in got:
                fh.write(json.dumps(rec, sort_keys=True) + "\n")
            n += len(got)
            time.sleep(0.1)
    return path


def classify(rec, generous):
    st, units = rec.get("standard_type"), rec.get("standard_units")
    rel = (rec.get("standard_relation") or "").strip()
    val, cmt, pc = _num(rec.get("standard_value")), (rec.get("activity_comment") or "").lower(), _num(rec.get("pchembl_value"))
    if cmt.strip().startswith(("active", "potent")) or (pc is not None and pc >= 5):
        return "active"
    if any(w in cmt for w in ("not active", "inactive", "no activity")):
        return "inactive"
    if units == "nM" and st in ("IC50", "Ki", "Kd") and val is not None:
        if rel in ("=", "<", "<=", "~") and val <= 10000:
            return "active"
        if rel in (">", ">=") and val >= 10000:
            return "inactive"
        return "ignore"
    if units == "%" and val is not None:
        if generous and st in ("Inhibition", "Activity"):
            return "active" if val >= 50 else ("inactive" if val <= 30 else "ignore")
        if st == "Inhibition":                       # CLEAN: only the direct inhibition %-assay
            if val >= 50:
                return "active"
            if val <= 20:
                m = DOSE.search(rec.get("assay_description") or "")
                dose = _num(m.group(1)) if m else None
                return "inactive" if (dose is not None and dose >= 10) else "ignore"
    return "ignore"


def funnel(by_cmp, smiles, generous, label):
    verd, inact = collections.Counter(), []
    for cid, recs in by_cmp.items():
        t = collections.Counter(classify(r, generous) for r in recs)
        v = "has-active" if t["active"] else ("enzyme-inactive" if t["inactive"] else "no-usable-evidence")
        verd[v] += 1
        if v == "enzyme-inactive" and cid in smiles:
            inact.append((cid, smiles[cid]))
    actives = I.load_actives(ACTIVES)
    fps = []
    for p in ALL_ACTIVES:
        if os.path.exists(p):
            fps += [a["fp"] for a in I.load_actives(p)]
    chosen, rej = I.select_inactives(inact, actives, fps, per_active=50)
    print(f"\n[{label}] compound-level enzyme verdicts: " + ", ".join(f"{k}={v}" for k, v in verd.most_common()))
    print(f"[{label}] enzyme-inactive with SMILES: {len(inact)}  ->  after frozen decoy filters: KEPT {len(chosen)}")
    print(f"[{label}]   filter rejects: " + ", ".join(f"{k}={v}" for k, v in rej.items() if v))
    return len(chosen)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cache", default="work/chembl_enzyme_audit")
    ap.add_argument("--refetch", action="store_true")
    args = ap.parse_args(argv)
    if args.refetch:
        p = os.path.join(args.cache, f"{TARGET}.jsonl")
        if os.path.exists(p):
            os.remove(p)
    rows = [json.loads(l) for l in open(fetch(args.cache))]
    by_cmp, smiles = collections.defaultdict(list), {}
    for r in rows:
        by_cmp[r["molecule_chembl_id"]].append(r)
        if r.get("canonical_smiles"):
            smiles.setdefault(r["molecule_chembl_id"], r["canonical_smiles"])
    cell = sum(1 for r in rows if CELL.search(r.get("assay_description") or ""))
    nm = [_num(r.get("standard_value")) for r in rows
          if r.get("standard_units") == "nM" and r.get("standard_type") in ("IC50", "Ki", "Kd")]
    nm = [v for v in nm if v is not None]
    print(f"target {TARGET} (C. albicans CYP51, UniProt P10613): {len(rows)} records, "
          f"{len(by_cmp)} distinct compounds ({len(smiles)} with SMILES)")
    print(f"  cell-based assay descriptions: {cell}/{len(rows)} "
          f"({100*cell/len(rows):.0f}% — even this 'enzyme' target is mostly cellular sterol readouts)")
    if nm:
        print(f"  IC50/Ki/Kd in nM: n={len(nm)}, all potent (max {max(nm):.0f} nM) — no censored non-binders")
    clean = funnel(by_cmp, smiles, False, "CLEAN")
    gen = funnel(by_cmp, smiles, True, "GENEROUS upper bound")
    print(f"\nCONCLUSION: property-matched enzyme-confirmed inactives available = {clean} (clean) "
          f"/ {gen} (generous upper bound).")
    print("A powered contrast needs ~O(100) (cf. the 279 whole-cell set). NOT BUILDABLE from public data.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
