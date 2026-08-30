"""Build the continuous-potency validation set: CYP51 ligands with measured pChEMBL potency (#81).

Pipeline stage: ligands/ (the continuous-potency axis, criterion-validation)

The binary active/inactive gates (validate_gate2.py, validate_gate_verified.py) show the
N-Fe geometry *enriches actives*. Issue #81 asks the stronger question: does it track *how
potent* a molecule is? That needs a continuous axis, not a label. This script assembles one
from ChEMBL's `pchembl_value` (-log10 molar potency), applies the SAME applicability triage
the front door uses, and emits the frozen dataset the grader (validate_gate_potency.py) then
correlates against docked geometry.

Two datasets, both pre-registered (work/PREREGISTRATION_potency.md):

    enzyme      target CHEMBL1780 (C. albicans CYP51 / sterol 14a-demethylase) — the DIRECT
                on-target binding potency. Small n, wide CI, but no whole-cell confound. PRIMARY.
    wholecell   target_organism "Candida albicans" — whole-cell growth-inhibition potency.
                Well-powered but confounded by uptake / efflux / target abundance. SECONDARY.

Two commands, split so the network half runs once (mirrors make_verified_inactives.py):

    fetch   pull every ChEMBL bioactivity record CARRYING A pchembl_value for each target
            into a local JSONL cache. Network-bound, resumable. Nothing is decided here.
    build   offline + deterministic: aggregate each compound to one median pChEMBL
            (openafr.potency.median_pchembl), run the applicability triage
            (scripts.applicability_triage.classify), and write a per-set .tsv + an
            in-domain .smi. Same cache -> same output.

Usage:
    conda run -n openafr python scripts/make_potency_set.py fetch [--cache DIR]
    conda run -n openafr python scripts/make_potency_set.py build [--cache DIR]
"""
import argparse
import collections
import hashlib
import json
import os
import pathlib
import sys
import time
import urllib.parse
import urllib.request

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import potency
from scripts.applicability_triage import classify, load_panel

API = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
PAGE = 1000
# Only records that carry a pchembl_value — the continuous axis. Everything else (censored,
# raw MIC, textual) is irrelevant to a potency correlation and would only bloat the cache.
QUERIES = {
    "enzyme": {"target_chembl_id": "CHEMBL1780", "pchembl_value__isnull": "false"},
    "wholecell": {"target_organism": "Candida albicans", "pchembl_value__isnull": "false"},
}
FIELDS = ("molecule_chembl_id,canonical_smiles,standard_type,standard_relation,"
          "standard_value,standard_units,pchembl_value,assay_chembl_id,"
          "document_chembl_id,target_chembl_id")
DEFAULT_CACHE = "work/chembl_potency"
LIG = "data/ligands"


def _get(params, retries=4):
    url = f"{API}?{urllib.parse.urlencode(params)}"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            if attempt == retries - 1:
                raise SystemExit(f"ChEMBL fetch failed after {retries} tries: {e}\n  {url}")
            time.sleep(2 ** attempt)


def cmd_fetch(args):
    os.makedirs(args.cache, exist_ok=True)
    manifest = {}
    for tag, sel in QUERIES.items():
        path = os.path.join(args.cache, f"{tag}.jsonl")
        params = dict(sel, limit=PAGE, offset=0, only=FIELDS, format="json")
        total = _get(dict(params, limit=1))["page_meta"]["total_count"]
        print(f"{tag}: {total} pchembl records -> {path}")
        n = 0
        with open(path, "w") as fh:
            while n < total:
                page = _get(dict(params, offset=n))
                rows = page.get("activities", [])
                if not rows:
                    break
                for rec in rows:
                    fh.write(json.dumps(rec, sort_keys=True) + "\n")
                n += len(rows)
                print(f"  {n}/{total}", end="\r", flush=True)
                time.sleep(0.1)
        digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
        manifest[tag] = {"query": sel, "records": n, "sha256": digest}
        print(f"  {n}/{total}  sha256 {digest[:16]}...")
    manifest["fetched_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(os.path.join(args.cache, "MANIFEST.json"), "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
    print(f"\ncache manifest: {os.path.join(args.cache, 'MANIFEST.json')}")


def _load(cache, tag):
    """cache/tag.jsonl -> ({cid: [records]}, {cid: smiles})."""
    path = os.path.join(cache, f"{tag}.jsonl")
    if not os.path.exists(path):
        raise SystemExit(f"missing cache {path} — run: {sys.argv[0]} fetch")
    by_cmp, smiles = collections.defaultdict(list), {}
    for line in open(path):
        rec = json.loads(line)
        cid = rec.get("molecule_chembl_id")
        if not cid:
            continue
        by_cmp[cid].append(rec)
        if rec.get("canonical_smiles"):
            smiles.setdefault(cid, rec["canonical_smiles"])
    return by_cmp, smiles


def _sort_key(cid):
    return int(cid.replace("CHEMBL", "")) if cid[6:].isdigit() else 10 ** 12


def _build_one(cache, tag, panel):
    by_cmp, smiles = _load(cache, tag)
    rows = []
    for cid in sorted(by_cmp, key=_sort_key):
        med, n_used, types = potency.median_pchembl(by_cmp[cid])
        if med is None or cid not in smiles:
            continue
        c = classify(smiles[cid], panel)
        in_domain = c["verdict"] in ("in-envelope-novel", "near-known")
        rows.append({
            "cid": cid, "smiles": smiles[cid], "pchembl": round(med, 3),
            "n_records": n_used, "types": "|".join(types),
            "heavy": c["heavy"], "verdict": c["verdict"],
            "max_tc": None if c["max_tc"] is None else round(c["max_tc"], 3),
            "nearest_active": c["nearest"], "in_domain": in_domain,
        })
    return rows


def _write(tag, rows):
    tsv = os.path.join(LIG, f"potency_{tag}.tsv")
    cols = ["cid", "smiles", "pchembl", "n_records", "types", "heavy", "verdict",
            "max_tc", "nearest_active", "in_domain"]
    with open(tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")
    dockable = [r for r in rows if r["in_domain"]]
    smi = os.path.join(LIG, f"potency_{tag}.smi")
    with open(smi, "w") as fh:
        for r in dockable:
            fh.write(f"{r['smiles']}\t{r['cid']}\n")
    tsv_sha = hashlib.sha256(open(tsv, "rb").read()).hexdigest()
    smi_sha = hashlib.sha256(open(smi, "rb").read()).hexdigest()
    n_novel = sum(r["verdict"] == "in-envelope-novel" for r in dockable)
    print(f"\n[{tag}] {len(rows)} compounds with a median pChEMBL")
    print(f"   in-domain (dockable) : {len(dockable)}  ({n_novel} novel, "
          f"{len(dockable) - n_novel} near-known)")
    print(f"   pChEMBL range        : {min(r['pchembl'] for r in rows):.2f} .. "
          f"{max(r['pchembl'] for r in rows):.2f}")
    print(f"   {os.path.basename(tsv)}  sha256 {tsv_sha[:16]}...")
    print(f"   {os.path.basename(smi)}  sha256 {smi_sha[:16]}...")


def cmd_build(args):
    panel = load_panel()
    for tag in QUERIES:
        _write(tag, _build_one(args.cache, tag, panel))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["fetch", "build"])
    ap.add_argument("--cache", default=DEFAULT_CACHE)
    args = ap.parse_args()
    (cmd_fetch if args.cmd == "fetch" else cmd_build)(args)


if __name__ == "__main__":
    main()
