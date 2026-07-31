import math

def read_atoms(path, model=None):
    out, in_model = [], (model is None)
    for line in open(path):
        if line.startswith("MODEL"):
            in_model = (model is not None and int(line.split()[1]) == model); continue
        if line.startswith("ENDMDL"):
            if model is not None and in_model: break
            continue
        if not in_model: continue
        if line.startswith(("ATOM","HETATM")):
            name = line[12:16].strip()
            elem = (line[76:78].strip() or name[0]).upper()
            if elem == "H": continue
            out.append((name, float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return out

ref = {a[0]: a[1:] for a in read_atoms("work/ref_VT1.pdb")}
fe = [(float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in open("work/receptor_A.pdb")
      if l.startswith("HETATM") and l[76:78].strip().upper()=="FE"][0]

scores = [float(l.split()[3]) for l in open("work/redock_VT1.pdbqt") if "REMARK VINA RESULT" in l]

print(f"{'pose':>5}{'score':>9}{'RMSD(name-matched)':>20}{'N-Fe':>9}   verdict")
print("-"*64)
results=[]
for i, sc in enumerate(scores, 1):
    d = read_atoms("work/redock_VT1.pdbqt", model=i)
    matched = [(n,xyz) for n,*rest in [(a[0],a[1:]) for a in d] for xyz in [rest[0]] if n in ref]
    pairs = [(ref[n], xyz) for n, xyz in matched]
    rmsd = math.sqrt(sum(math.dist(a,b)**2 for a,b in pairs)/len(pairs))
    ns = [a for a in d if a[0].startswith("N")]
    dfe = min(math.dist(a[1:], fe) for a in ns)
    v = "*** CORRECT ***" if (rmsd < 2.0 and dfe < 3.0) else ("good pose" if rmsd < 2.0 else ("iron-bound" if dfe<3.0 else ""))
    print(f"{i:>5}{sc:>9.2f}{rmsd:>20.2f}{dfe:>9.2f}   {v}")
    results.append((i,sc,rmsd,dfe,len(pairs)))
print(f"\n(matched {results[0][4]}/{len(ref)} heavy atoms by name)")
b = min(results, key=lambda r: r[2])
print(f"Best pose by position: #{b[0]} — RMSD {b[2]:.2f} A, N-Fe {b[3]:.2f} A, score {b[1]:.2f}")
top = results[0]
print(f"Vina's TOP-RANKED pose:  #1 — RMSD {top[2]:.2f} A, N-Fe {top[3]:.2f} A")
print()
print("VERDICT: " + ("PASS — top pose reproduces the crystal structure (<2 A)" if top[2] < 2.0
      else "PARTIAL — a correct pose exists but is not ranked first" if b[2] < 2.0
      else "FAIL — no pose reproduces the crystal structure"))
