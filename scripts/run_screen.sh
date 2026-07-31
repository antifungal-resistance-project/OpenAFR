#!/bin/bash
# Pipeline stage: screen/ (spike version)
# ONE protocol for every ligand. Actives and decoys are indistinguishable to this
# script — any per-ligand tuning would inflate enrichment and invalidate the test.
set -u
CX=67.4; CY=68.8; CZ=3.4     # midpoint of heme iron and the crystal ligand pocket
SIZE=30; EXHAUST=16; SEED=42
OUT=work/screen; LIG=work/screen_pdbqt
mkdir -p "$OUT" "$LIG"
total=$(ls work/screen_ligands/*.sdf | wc -l | tr -d ' ')
i=0; failed=0
for sdf in work/screen_ligands/*.sdf; do
  name=$(basename "$sdf" .sdf); i=$((i+1))
  [ -s "$OUT/${name}.pdbqt" ] && continue                      # resumable
  obabel "$sdf" -O "$LIG/${name}.pdbqt" -p 7.4 >/dev/null 2>&1
  if [ ! -s "$LIG/${name}.pdbqt" ]; then
    echo "CONVERT_FAIL ${name}"; failed=$((failed+1)); continue
  fi
  vina --receptor work/receptor.pdbqt --ligand "$LIG/${name}.pdbqt" \
       --center_x $CX --center_y $CY --center_z $CZ \
       --size_x $SIZE --size_y $SIZE --size_z $SIZE \
       --exhaustiveness $EXHAUST --seed $SEED \
       --out "$OUT/${name}.pdbqt" >/dev/null 2>&1
  [ -s "$OUT/${name}.pdbqt" ] || { echo "DOCK_FAIL ${name}"; failed=$((failed+1)); }
  if [ $((i % 25)) -eq 0 ]; then
    echo "progress: $i/$total  docked=$(ls $OUT/*.pdbqt 2>/dev/null | wc -l | tr -d ' ')  failed=$failed"
  fi
done
echo "SCREEN_COMPLETE docked=$(ls $OUT/*.pdbqt 2>/dev/null | wc -l | tr -d ' ') failed=$failed"
