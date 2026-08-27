#!/bin/bash
# Human CYP51 counter-screen (issue #73) — dock a ligand set into the HUMAN CYP51 receptor
# (scripts/prep_receptor_human.py, from 3LD6) under the SAME frozen protocol as the fungal
# screen, so the only difference between a fungal and a human score is the protein.
#
# It reuses the consensus engine (scripts/screen_consensus.sh -> screen.sh) unchanged; it only
# sets two things via the documented env overrides:
#   RECEPTOR            -> work/receptor_human.pdbqt      (the human protein)
#   CENTER_{X,Y,Z}_OVERRIDE -> the human active site, read from work/receptor_human_box.json
#                             (the crystallographic ketoconazole centroid, derived at prep time)
# Size, exhaustiveness, num_modes and the seed list stay from the frozen protocol.
#
# Usage: scripts/screen_human.sh LIGAND_DIR OUT_DIR [SEEDS]
set -eu
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

LIG_SDF="${1:?usage: screen_human.sh LIGAND_DIR OUT_DIR [SEEDS]}"
OUT="${2:?usage: screen_human.sh LIGAND_DIR OUT_DIR [SEEDS]}"
SEEDS="${3:-42 1 2 3 4}"
BOX="work/receptor_human_box.json"
[ -s "$BOX" ] || { echo "missing $BOX — run: python scripts/prep_receptor_human.py" >&2; exit 2; }
[ -s work/receptor_human.pdbqt ] || { echo "missing work/receptor_human.pdbqt — run: python scripts/prep_receptor_human.py" >&2; exit 2; }

# Pull the derived box center out of the JSON (python, no jq dependency).
read -r CX CY CZ < <(python3 -c "import json;d=json.load(open('$BOX'));print(d['center_x'],d['center_y'],d['center_z'])")
echo "human counter-screen: receptor=work/receptor_human.pdbqt box=($CX,$CY,$CZ) seeds=[$SEEDS]"

export RECEPTOR="work/receptor_human.pdbqt"
export CENTER_X_OVERRIDE="$CX" CENTER_Y_OVERRIDE="$CY" CENTER_Z_OVERRIDE="$CZ"
"$ROOT/scripts/screen_consensus.sh" "$LIG_SDF" "$OUT" "$SEEDS"
echo "HUMAN_SCREEN_COMPLETE out=$OUT"
