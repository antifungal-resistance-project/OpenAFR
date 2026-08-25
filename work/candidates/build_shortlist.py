"""Join top-N ranked repurposing hits with their SMILES -> a check_precedent input."""
import sys

TOPN = int(sys.argv[1]) if len(sys.argv) > 1 else 40

smi = {}
with open("work/repurposing/library_in_domain.smi") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 2:
            smi[parts[1]] = parts[0]

rows = []
with open("work/repurposing/ranked.tsv") as f:
    next(f)
    for line in f:
        rank, ligand, nfe, vina = line.rstrip("\n").split("\t")
        rows.append((int(rank), ligand, nfe, vina))
        if len(rows) >= TOPN:
            break

miss = [l for _, l, _, _ in rows if l not in smi]
with open("work/candidates/top%d.smi" % TOPN, "w") as out:
    for rank, ligand, nfe, vina in rows:
        if ligand in smi:
            out.write("%s\t%s\n" % (smi[ligand], ligand))
print("wrote %d of %d (missing SMILES: %s)" % (
    sum(1 for _, l, _, _ in rows if l in smi), len(rows), miss))
