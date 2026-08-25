"""Merge ranked repurposing hits with the precedent verdicts into one readable table.

Usage: merge_view.py [PRECEDENT.tsv]   (default work/candidates/top40_precedent.tsv)
Shows the near-neighbour columns when the precedent file has them.
"""
import sys

PREC = sys.argv[1] if len(sys.argv) > 1 else "work/candidates/top40_precedent.tsv"

KNOWN_AZOLES = {  # controls, not discoveries (from RESULTS_repurposing.md)
    "flutrimazole", "ravuconazole", "fosfluconazole", "fluconazole",
    "voriconazole", "luliconazole", "croconazole", "letrozole", "clotrimazole",
}

rank = {}
with open("work/repurposing/ranked.tsv") as f:
    next(f)
    for line in f:
        r, lig, nfe, vina = line.rstrip("\n").split("\t")
        rank[lig] = (int(r), nfe, vina)

prec = {}
with open(PREC) as f:
    hdr = next(f).rstrip("\n").split("\t")
    for line in f:
        p = line.rstrip("\n").split("\t")
        d = dict(zip(hdr, p))
        prec[d["name"]] = d

has_near = "near_verdict" in hdr
rows = sorted(prec.values(), key=lambda d: rank[d["name"]][0])
cols = ["rank", "name", "N-Fe", "enzyme_verdict", "phenotypic", "off"]
if has_near:
    cols += ["near_sim", "near_verdict"]
cols += ["ctrl?"]
w = " ".join("%-{}s".format(n) for n in (4, 16, 7, 19, 16, 5, 8, 16, 5)[:len(cols)])
print(w % tuple(cols))
print("-" * 100)
for d in rows:
    n = d["name"]
    r, nfe, vina = rank[n]
    ctrl = "CTRL" if n in KNOWN_AZOLES else ""
    vals = [r, n[:16], nfe, d["verdict"], d.get("phenotypic", ""), d.get("n_off_target", "")]
    if has_near:
        vals += [d.get("near_sim", ""), d.get("near_verdict", "")]
    vals += [ctrl]
    print(w % tuple(vals))
