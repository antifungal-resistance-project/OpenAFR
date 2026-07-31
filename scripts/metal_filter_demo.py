import math
def read_atoms(path, model=None):
    out, in_model = [], (model is None)
    for line in open(path):
        if line.startswith("MODEL"):
            in_model = (model is not None and int(line.split()[1])==model); continue
        if line.startswith("ENDMDL"):
            if model is not None and in_model: break
            continue
        if not in_model: continue
        if line.startswith(("ATOM","HETATM")):
            nm=line[12:16].strip(); el=(line[76:78].strip() or nm[0]).upper()
            if el=="H": continue
            out.append((nm,float(line[30:38]),float(line[38:46]),float(line[46:54])))
    return out

ref={a[0]:a[1:] for a in read_atoms("work/ref_VT1.pdb")}
fe=[(float(l[30:38]),float(l[38:46]),float(l[46:54])) for l in open("work/receptor_A.pdb")
    if l.startswith("HETATM") and l[76:78].strip().upper()=="FE"][0]
scores=[float(l.split()[3]) for l in open("work/redock_VT1.pdbqt") if "REMARK VINA RESULT" in l]

poses=[]
for i,sc in enumerate(scores,1):
    d=read_atoms("work/redock_VT1.pdbqt",model=i)
    pairs=[(ref[n],xyz) for n,*r in [(a[0],a[1:]) for a in d] for xyz in [r[0]] if n in ref]
    rmsd=math.sqrt(sum(math.dist(a,b)**2 for a,b in pairs)/len(pairs))
    dfe=min(math.dist(a[1:],fe) for a in d if a[0].startswith("N"))
    poses.append((i,sc,rmsd,dfe))

print("WITHOUT the chemistry constraint — trust Vina's #1 ranking:")
p=poses[0]
print(f"   pick pose #{p[0]}  ->  RMSD {p[2]:.2f} A   {'CORRECT' if p[2]<2 else 'WRONG (misses the iron bond)'}")
print()
print("WITH the constraint — 'the azole nitrogen must bond the heme iron (<3 A)':")
kept=[p for p in poses if p[3]<3.0]
print(f"   {len(poses)} poses in  ->  {len(kept)} survive the filter")
for p in kept: print(f"      pose #{p[0]}: score {p[1]:.2f}, N-Fe {p[3]:.2f} A, RMSD {p[2]:.2f} A")
best=max(kept,key=lambda p:-p[1]) if kept else None
print()
if best:
    print(f"   best surviving pose: #{best[0]}  ->  RMSD {best[2]:.2f} A   {'*** CORRECT ***' if best[2]<2 else 'wrong'}")
