"""Build the verified-inactive decoy set: compounds MEASURED not to inhibit C. albicans.

Pipeline stage: ligands/ (decoy half, v2)

This is the independent dataset the preprint's stopping rule requires
(work/PREPRINT_geometry_ceiling.md §5: "a further attempt requires an independent dataset
— new actives or a new decoy construction"), and the exit §4.1 names for the
decoy-construction bind. The evidence rule and every filter live in openafr/inactives.py;
this script is the network + I/O shell around them.

Two commands, deliberately split so the expensive half runs once:

    fetch   pull every ChEMBL bioactivity record for C. albicans (whole-cell) and for the
            C. albicans CYP51 enzyme target into a local JSONL cache. Network-bound,
            resumable, ~88 pages. Nothing is decided here.
    build   offline: classify each compound from the cache, apply the frozen filters, and
            write the .smi + a provenance stanza. Deterministic — same cache, same output.

Usage:
    python scripts/make_verified_inactives.py fetch [--cache DIR]
    python scripts/make_verified_inactives.py build [--cache DIR] [--out FILE]
                                                    [--actives FILE] [--per-active N]
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
from openafr import inactives

API = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
PAGE = 1000
# Whole-cell C. albicans (target_organism string) + the C. albicans CYP51 enzyme target,
# which carries the strain-qualified organism name and so is NOT covered by the first pull.
QUERIES = {
    "calbicans_cell": {"target_organism": "Candida albicans"},
    "calbicans_cyp51": {"target_chembl_id": "CHEMBL1780"},
}
# Only the fields the evidence rule reads, so the cache stays small and auditable.
FIELDS = ("molecule_chembl_id,canonical_smiles,standard_type,standard_relation,"
          "standard_value,standard_units,activity_comment,pchembl_value,"
          "assay_chembl_id,document_chembl_id,target_chembl_id")
DEFAULT_CACHE = "work/chembl"
DEFAULT_OUT = "data/ligands/verified_inactives.smi"
DEFAULT_ACTIVES = "data/ligands/actives_holdout_final.smi"
ALL_ACTIVES = ("data/ligands/actives.smi", "data/ligands/actives_holdout_final.smi")
PER_ACTIVE = 50


def _get(params, retries=4):
    url = f"{API}?{urllib.parse.urlencode(params)}"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except Exception as e:
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
        print(f"{tag}: {total} records -> {path}")
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


def load_cache(cache):
    """cache dir -> ({compound_id: [records]}, {compound_id: smiles}, source counts)."""
    by_cmp, smiles, counts = collections.defaultdict(list), {}, {}
    for tag in QUERIES:
        path = os.path.join(cache, f"{tag}.jsonl")
        if not os.path.exists(path):
            raise SystemExit(f"missing cache {path} — run: {sys.argv[0]} fetch")
        n = 0
        for line in open(path):
            rec = json.loads(line)
            cid = rec.get("molecule_chembl_id")
            if not cid:
                continue
            by_cmp[cid].append(rec)
            if rec.get("canonical_smiles"):
                smiles.setdefault(cid, rec["canonical_smiles"])
            n += 1
        counts[tag] = n
    return by_cmp, smiles, counts


def _chembl_sort_key(cid):
    """Deterministic, human-meaningful ordering: numeric ChEMBL id ascending."""
    return int(cid.replace("CHEMBL", "")) if cid[6:].isdigit() else 10**12


def cmd_build(args):
    by_cmp, smiles, counts = load_cache(args.cache)
    print("cache: " + "  ".join(f"{k}={v} records" for k, v in counts.items())
          + f"  ({len(by_cmp)} distinct compounds)")

    verdicts = collections.Counter()
    verified = []
    for cid in sorted(by_cmp, key=_chembl_sort_key):
        verdict, _tally = inactives.compound_verdict(by_cmp[cid])
        verdicts[verdict] += 1
        if verdict == "verified-inactive" and cid in smiles:
            verified.append((cid, smiles[cid]))
    print("\nevidence verdicts over every compound ever tested on C. albicans:")
    for k, v in verdicts.most_common():
        print(f"  {k:22} {v:6d}")
    print(f"  -> {len(verified)} verified inactives carrying a structure")

    actives = inactives.load_actives(args.actives)
    # Novelty is measured against EVERY active the project uses (training + held-out), so a
    # decoy can never be a near-neighbour of a molecule the criterion was defined on.
    all_fps = []
    for path in ALL_ACTIVES:
        all_fps += [a["fp"] for a in inactives.load_actives(path)]

    chosen, rejects = inactives.select_inactives(verified, actives, all_fps, args.per_active)
    print(f"\nfilter funnel ({len(verified)} verified inactives in):")
    for k, v in rejects.items():
        print(f"  rejected {k:18} {v:6d}")
    print(f"  -> KEPT {len(chosen)}")

    per = collections.Counter(c["paired"] for c in chosen)
    print("\npaired per active:")
    for a in actives:
        print(f"  {a['name']:16} {per.get(a['name'], 0):3d}   "
              f"(MW {a['mw']:.0f}, logP {a['logp']:.1f}, heavy {a['heavy']})")

    if chosen:
        import statistics as st
        print("\nproperty match check (verified inactives vs held-out actives):")
        for key, label in (("mw", "MW"), ("logp", "cLogP"),
                           ("heavy", "heavy atoms"), ("rotb", "rot. bonds")):
            av = st.mean(a[key] for a in actives)
            dv = st.mean(c[key] for c in chosen)
            print(f"  {label:12} actives {av:8.1f}   inactives {dv:8.1f}   diff {abs(av - dv):6.1f}")

    with open(args.out, "w") as fh:
        for c in chosen:
            fh.write(f"{c['smiles']}\tinactive_{c['chembl_id']}_for_{c['paired']}\n")
    digest = hashlib.sha256(open(args.out, "rb").read()).hexdigest()
    print(f"\nwrote {len(chosen)} verified inactives to {args.out}")
    print(f"sha256 {digest}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fetch", help="pull ChEMBL bioactivity into a local cache")
    f.add_argument("--cache", default=DEFAULT_CACHE)
    f.set_defaults(func=cmd_fetch)
    b = sub.add_parser("build", help="offline: classify + filter + write the .smi")
    b.add_argument("--cache", default=DEFAULT_CACHE)
    b.add_argument("--out", default=DEFAULT_OUT)
    b.add_argument("--actives", default=DEFAULT_ACTIVES)
    b.add_argument("--per-active", type=int, default=PER_ACTIVE)
    b.set_defaults(func=cmd_build)
    args = ap.parse_args(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    sys.exit(main())
