#!/bin/bash
# Pipeline stage: screen/ — the single place a ligand library gets docked.
#
# ONE frozen protocol for every molecule. Box + search parameters come from the
# hash-frozen protocol.yaml (the same source the gate reads), never hand-copied
# constants, so what is docked can never silently drift from what is graded.
# Every ligand is docked with the same fixed --seed and no per-ligand tuning —
# any tuning would inflate enrichment and invalidate the validation guarantee.
#
# Docking runs are independent, so ligands are dispatched in parallel, one Vina
# process per core. The fixed protocol seed is applied identically and
# independently to every ligand — there is no shared RNG across the pool — so a
# parallel run is bit-for-bit the same as a serial one, just faster.
#
# Usage: scripts/screen.sh [LIGAND_DIR] [OUT_DIR] [JOBS]
#   LIGAND_DIR  directory of prepared *.sdf ligands  (default work/screen_ligands)
#   OUT_DIR     where docked *.pdbqt poses land       (default work/screen)
#   JOBS        parallel Vina processes               (default: all cores)
#
# Resumable: a ligand whose OUT_DIR/<name>.pdbqt already exists is skipped, so an
# interrupted screen picks up where it stopped. Per-ligand output is captured to
# OUT_DIR/<name>.log; failures are counted and reported, never silently dropped.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

LIG_SDF="${1:-work/screen_ligands}"
OUT="${2:-work/screen}"
# core count: nproc (Linux) or sysctl (macOS / Apple Silicon), fall back to 4.
JOBS="${3:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)}"
PDBQT="${OUT}_pdbqt"
# Receptor to dock into. Defaults to the wild-type work/receptor.pdbqt so existing
# runs are unchanged; set RECEPTOR to dock a mutant (e.g. the C. auris resistance
# port, work/receptor_auris_Y132F.pdbqt). The box, seed and search effort still come
# from the frozen protocol — only the receptor structure changes.
RECEPTOR="${RECEPTOR:-work/receptor.pdbqt}"

# Box + search effort + seed, straight from the frozen protocol (single source of
# truth shared with the gate). Sets CX CY CZ SIZE EXHAUST MODES SEED.
eval "$(python openafr/protocol.py --shell docking)"
# Seed override for consensus/replicate runs (scripts/screen_consensus.sh). Unset =
# the frozen protocol seed, so a normal single-seed screen is byte-identical to before.
# Only the RNG seed changes; box, exhaustiveness and modes stay from the frozen protocol.
SEED="${SEED_OVERRIDE:-$SEED}"

# Box-center override for docking into a DIFFERENT receptor's active site (e.g. the human
# CYP51 counter-screen, scripts/prep_receptor_human.py). Unset = the frozen fungal box, so a
# normal screen is byte-identical to before. Only the box CENTER moves to the other pocket;
# size, exhaustiveness, num_modes and seed stay from the frozen protocol. Mirrors SEED_OVERRIDE:
# a different receptor needs its box centered on ITS ligand site, never the 5TZ1 centroid.
CX="${CENTER_X_OVERRIDE:-$CX}"
CY="${CENTER_Y_OVERRIDE:-$CY}"
CZ="${CENTER_Z_OVERRIDE:-$CZ}"

[ -s "$RECEPTOR" ] || { echo "missing receptor: $RECEPTOR — run: python scripts/prep_receptor.py" >&2; exit 2; }
[ -d "$LIG_SDF" ] || { echo "no ligand dir: $LIG_SDF" >&2; exit 2; }
total=$(ls "$LIG_SDF"/*.sdf 2>/dev/null | wc -l | tr -d ' ')
[ "$total" -gt 0 ] || { echo "no *.sdf ligands in $LIG_SDF" >&2; exit 2; }
mkdir -p "$OUT" "$PDBQT"

echo "screen: $total ligands  receptor=$RECEPTOR  jobs=$JOBS  box=($CX,$CY,$CZ) size=$SIZE exhaust=$EXHAUST modes=$MODES seed=$SEED"

# One self-contained docking job per ligand: convert to PDBQT at pH 7.4, then dock
# under the frozen protocol. Reads only its own files, so parallel copies never
# contend for shared state.
dock_one() {
  local sdf="$1" name log
  name=$(basename "$sdf" .sdf)
  log="$OUT/${name}.log"
  [ -s "$OUT/${name}.pdbqt" ] && { echo "SKIP  ${name}"; return 0; }
  obabel "$sdf" -O "$PDBQT/${name}.pdbqt" -p 7.4 >"$log" 2>&1
  [ -s "$PDBQT/${name}.pdbqt" ] || { echo "CONVERT_FAIL ${name}"; return 0; }
  vina --receptor "$RECEPTOR" --ligand "$PDBQT/${name}.pdbqt" \
       --center_x "$CX" --center_y "$CY" --center_z "$CZ" \
       --size_x "$SIZE" --size_y "$SIZE" --size_z "$SIZE" \
       --exhaustiveness "$EXHAUST" --num_modes "$MODES" --seed "$SEED" \
       --out "$OUT/${name}.pdbqt" >>"$log" 2>&1
  [ -s "$OUT/${name}.pdbqt" ] && echo "OK    ${name}" || echo "DOCK_FAIL ${name}"
}
export -f dock_one
export OUT PDBQT RECEPTOR CX CY CZ SIZE EXHAUST MODES SEED

ls "$LIG_SDF"/*.sdf | xargs -P "$JOBS" -I {} bash -c 'dock_one "$1"' _ {}

docked=$(ls "$OUT"/*.pdbqt 2>/dev/null | wc -l | tr -d ' ')
echo "SCREEN_COMPLETE docked=${docked}/${total} failed=$((total - docked)) out=${OUT}"
