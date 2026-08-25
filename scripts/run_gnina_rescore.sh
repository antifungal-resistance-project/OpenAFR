#!/bin/bash
# Pipeline stage: validate/ — GNINA metal-aware baseline rescore of the run-2 poses.
# Implements the run step of work/PREREGISTRATION_baselines.md.
#
# Runs GNINA in --score_only mode over the EXACT run-2 held-out poses (no re-docking):
# for each Vina pose file <name>.pdbqt it saves GNINA's stdout to <OUT>/<name>.gnina.
# scripts/rescore_gnina.py then parses those saved files + the same Vina poses and grades.
#
# Splitting run (this script, GNINA-capable host only) from grade (rescore_gnina.py, runs
# anywhere with rdkit) is deliberate: GNINA is a compiled CNN binary that is Linux/x86
# native — the same osx-arm64 wall the ERG11 re-caller hits — so the numeric run happens on
# a Linux box while the graded, hash-checked scoring stays reproducible on any host.
#
# Usage: scripts/run_gnina_rescore.sh [POSE_DIR] [OUT_DIR] [RECEPTOR]
#   POSE_DIR   dir of Vina .pdbqt poses to rescore   (default work/screen2)
#   OUT_DIR    where per-ligand .gnina outputs land   (default work/gnina_scores)
#   RECEPTOR   receptor GNINA scores against          (default work/receptor.pdbqt)
#
# Resumable: a ligand whose OUT_DIR/<name>.gnina already exists is skipped.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

POSE_DIR="${1:-work/screen2}"
OUT="${2:-work/gnina_scores}"
RECEPTOR="${3:-work/receptor.pdbqt}"

command -v gnina >/dev/null 2>&1 || {
  echo "gnina not on PATH — this step needs a GNINA-capable (Linux/x86) host." >&2
  echo "See work/PREREGISTRATION_baselines.md; grading runs anywhere once *.gnina exist." >&2
  exit 2; }
[ -s "$RECEPTOR" ] || { echo "missing receptor: $RECEPTOR" >&2; exit 2; }
[ -d "$POSE_DIR" ] || { echo "no pose dir: $POSE_DIR — run scripts/run_screen2.sh first" >&2; exit 2; }
total=$(ls "$POSE_DIR"/*.pdbqt 2>/dev/null | wc -l | tr -d ' ')
[ "$total" -gt 0 ] || { echo "no *.pdbqt poses in $POSE_DIR" >&2; exit 2; }
mkdir -p "$OUT"

# Record the exact GNINA build so the comparison is reproducible (prereg limitation 3).
gnina --version 2>&1 | head -3 | sed 's/^/gnina: /' | tee "$OUT/GNINA_VERSION.txt"

echo "gnina rescore: $total ligands  receptor=$RECEPTOR  out=$OUT"
scored=0
for pose in "$POSE_DIR"/*.pdbqt; do
  name=$(basename "$pose" .pdbqt)
  dest="$OUT/${name}.gnina"
  [ -s "$dest" ] && { echo "SKIP  ${name}"; scored=$((scored + 1)); continue; }
  # --score_only: score the poses as given, no minimization, no re-docking. Default CNN
  # model as shipped (prereg: opponent used as-is, recorded above). One block per pose,
  # in input-file order, so it aligns 1:1 with openafr.pdbqt.read_poses over the same file.
  if gnina --receptor "$RECEPTOR" --ligand "$pose" --score_only >"$dest" 2>"$dest.log"; then
    echo "OK    ${name}"; scored=$((scored + 1))
  else
    echo "GNINA_FAIL ${name} (see $dest.log)"; rm -f "$dest"
  fi
done

echo "GNINA_RESCORE_COMPLETE scored=${scored}/${total} out=${OUT}"
echo "next: python scripts/rescore_gnina.py $OUT $POSE_DIR <receptor.pdb with heme, e.g. work/receptor_A.pdb>"
