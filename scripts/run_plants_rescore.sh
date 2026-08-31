#!/bin/bash
# Pipeline stage: validate/ — PLANTS/ChemPLP metal-aware baseline rescore of the run-2 poses.
# Implements the PLANTS run step of work/PREREGISTRATION_baselines_multi.md (#83).
#
# Rescores (NOT re-docks) the EXACT held-out poses with ChemPLP, whose scoring function carries
# explicit metal-acceptor terms — the second metal-aware opponent alongside GNINA. For each
# Vina pose file <name>.pdbqt it writes PLANTS' per-pose features.csv to <OUT>/<name>.csv, which
# scripts/rescore_baselines.py parses (openafr.baselines.parse_plants_features) and grades.
#
# Split run (this script, PLANTS-capable Linux host) from grade (rescore_baselines.py, any host)
# for the same reason as GNINA: PLANTS + SPORES are Linux/x86 binaries — the osx-arm64 wall — so
# the numeric run happens on a cloud box while the hash-checked scoring stays reproducible.
#
# Usage: scripts/run_plants_rescore.sh [POSE_DIR] [OUT_DIR] [RECEPTOR_MOL2]
#   POSE_DIR       dir of Vina .pdbqt poses to rescore   (default work/screen2)
#   OUT_DIR        where per-ligand features.csv land     (default work/plants_scores)
#   RECEPTOR_MOL2  SPORES-prepared receptor .mol2         (default work/receptor_plants.mol2)
#
# Resumable: a ligand whose OUT_DIR/<name>.csv already exists is skipped.
#
# OPERATOR NOTE (cloud step, verify once): PLANTS rescoring needs mol2 inputs and a bindingsite
# centred on the frozen box (protocol.yaml: 70.61/66.28/4.18, radius derived from the 26 A box).
# Prepare the receptor with SPORES (`SPORES --mode complete receptor.pdb receptor_plants.mol2`)
# and confirm the `--mode rescore` config below matches your PLANTS build's flags before the run.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

POSE_DIR="${1:-work/screen2}"
OUT="${2:-work/plants_scores}"
RECEPTOR="${3:-work/receptor_plants.mol2}"
# Bindingsite: the frozen box centre; radius = half the 26 A box side (a sphere covering it).
CENTER="70.61 66.28 4.18"
RADIUS=13

command -v PLANTS >/dev/null 2>&1 || {
  echo "PLANTS not on PATH — this step needs a PLANTS-capable (Linux/x86) host." >&2
  echo "See work/PREREGISTRATION_baselines_multi.md; grading runs anywhere once *.csv exist." >&2
  exit 2; }
command -v SPORES >/dev/null 2>&1 || echo "warning: SPORES not on PATH — mol2 prep may be needed." >&2
[ -s "$RECEPTOR" ] || { echo "missing receptor mol2: $RECEPTOR (prepare with SPORES)" >&2; exit 2; }
[ -d "$POSE_DIR" ] || { echo "no pose dir: $POSE_DIR — run scripts/run_screen2.sh first" >&2; exit 2; }
total=$(ls "$POSE_DIR"/*.pdbqt 2>/dev/null | wc -l | tr -d ' ')
[ "$total" -gt 0 ] || { echo "no *.pdbqt poses in $POSE_DIR" >&2; exit 2; }
mkdir -p "$OUT"

PLANTS --version 2>&1 | head -3 | sed 's/^/plants: /' | tee "$OUT/PLANTS_VERSION.txt"
echo "plants rescore: $total ligands  receptor=$RECEPTOR  site=($CENTER) r=$RADIUS  out=$OUT"

scored=0
for pose in "$POSE_DIR"/*.pdbqt; do
  name=$(basename "$pose" .pdbqt)
  dest="$OUT/${name}.csv"
  [ -s "$dest" ] && { echo "SKIP  ${name}"; scored=$((scored + 1)); continue; }
  work=$(mktemp -d)
  # PLANTS reads mol2 ligands; convert the Vina poses (all modes) preserving order.
  obabel "$pose" -O "$work/lig.mol2" >/dev/null 2>&1
  cat > "$work/rescore.conf" <<CONF
scoring_function chemplp
protein_file $RECEPTOR
ligand_file $work/lig.mol2
bindingsite_center $CENTER
bindingsite_radius $RADIUS
output_dir $work/out
write_multi_mol2 0
CONF
  if PLANTS --mode rescore "$work/rescore.conf" >"$work/plants.log" 2>&1 \
        && [ -s "$work/out/features.csv" ]; then
    cp "$work/out/features.csv" "$dest"; echo "OK    ${name}"; scored=$((scored + 1))
  else
    echo "PLANTS_FAIL ${name} (see $work/plants.log)"
  fi
  rm -rf "$work"
done

echo "PLANTS_RESCORE_COMPLETE scored=${scored}/${total} out=${OUT}"
echo "next: python scripts/rescore_baselines.py $POSE_DIR work/receptor_A.pdb"
