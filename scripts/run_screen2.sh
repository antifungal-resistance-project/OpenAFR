#!/bin/bash
# Run 2 — protocol FIXED BY work/PREREGISTRATION_run2.md. Do not tune.
set -u
CX=70.61; CY=66.28; CZ=4.18; SIZE=26; EXHAUST=32; MODES=20; SEED=42
OUT=work/screen2; LIG=work/run2_pdbqt; mkdir -p "$OUT" "$LIG"
total=$(ls work/run2_ligands/*.sdf | wc -l | tr -d ' '); i=0; failed=0
for sdf in work/run2_ligands/*.sdf; do
  name=$(basename "$sdf" .sdf); i=$((i+1))
  [ -s "$OUT/${name}.pdbqt" ] && continue
  obabel "$sdf" -O "$LIG/${name}.pdbqt" -p 7.4 >/dev/null 2>&1
  [ -s "$LIG/${name}.pdbqt" ] || { echo "CONVERT_FAIL ${name}"; failed=$((failed+1)); continue; }
  vina --receptor work/receptor.pdbqt --ligand "$LIG/${name}.pdbqt" \
       --center_x $CX --center_y $CY --center_z $CZ \
       --size_x $SIZE --size_y $SIZE --size_z $SIZE \
       --exhaustiveness $EXHAUST --num_modes $MODES --seed $SEED \
       --out "$OUT/${name}.pdbqt" >/dev/null 2>&1
  [ -s "$OUT/${name}.pdbqt" ] || { echo "DOCK_FAIL ${name}"; failed=$((failed+1)); }
  [ $((i % 25)) -eq 0 ] && echo "progress: $i/$total failed=$failed"
done
echo "RUN2_COMPLETE docked=$(ls $OUT/*.pdbqt 2>/dev/null | wc -l | tr -d ' ') failed=$failed"
