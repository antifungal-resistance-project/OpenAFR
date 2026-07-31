#!/bin/bash
# Pipeline stage: dock/  — dock every prepared ligand into the receptor.
set -u
CENTER_X=70.61; CENTER_Y=66.28; CENTER_Z=4.18
SIZE=26                       # 26 A box: fits the long azoles (posaconazole ~25 A)
EXHAUST=16; SEED=42           # fixed seed = reproducible runs
mkdir -p work/docked
for sdf in work/ligands/*.sdf; do
  name=$(basename "$sdf" .sdf)
  obabel "$sdf" -O "work/ligands/${name}.pdbqt" -p 7.4 >/dev/null 2>&1
  start=$(date +%s)
  vina --receptor work/receptor.pdbqt --ligand "work/ligands/${name}.pdbqt" \
       --center_x $CENTER_X --center_y $CENTER_Y --center_z $CENTER_Z \
       --size_x $SIZE --size_y $SIZE --size_z $SIZE \
       --exhaustiveness $EXHAUST --seed $SEED \
       --out "work/docked/${name}.pdbqt" >"work/docked/${name}.log" 2>&1
  rc=$?
  el=$(( $(date +%s) - start ))
  best=$(grep -m1 "REMARK VINA RESULT" "work/docked/${name}.pdbqt" 2>/dev/null | awk '{print $4}')
  echo "${name}: rc=${rc} best=${best:-none} ${el}s"
done
echo "BATCH_COMPLETE"
