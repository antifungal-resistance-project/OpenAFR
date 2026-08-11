#!/bin/bash
# Consensus screen — dock every ligand under R independent Vina searches (different
# seeds) so the min nitrogen-to-iron distance can be estimated with less search noise.
#
# Why: the single-search min N-Fe is a noisy estimator, especially for flexible actives —
# a real azole's best coordinating pose is hard to find in one search (measured spread up
# to ~1.1 A across seeds), while rigid decoys sample their floor every time. See
# .context/consensus_diagnostic.md. Docking each ligand R times and taking the BEST
# coordination is a less-biased estimate of the achievable minimum (which is what the
# criterion always meant to measure), and the run-to-run spread becomes a confidence signal.
#
# This does NOT redefine the protocol: box, exhaustiveness and num_modes stay from the
# hash-frozen protocol.yaml; ONLY the RNG seed varies across replicates. It reuses the one
# docking engine (scripts/screen.sh) unchanged — one normal single-seed screen per seed,
# each into its own subdirectory OUT/s<seed>/. The consensus grader
# (scripts/validate_gate_consensus.py) reads those subdirectories.
#
# Usage: scripts/screen_consensus.sh LIGAND_DIR OUT_DIR [SEEDS]
#   LIGAND_DIR  directory of prepared *.sdf ligands
#   OUT_DIR     parent dir; per-seed poses land in OUT_DIR/s<seed>/
#   SEEDS       space-separated seed list (default: "42 1 2 3 4"; 42 = the frozen seed,
#               so the seed-42 slice reproduces the single-seed baseline exactly)
#   RECEPTOR    (env) receptor pdbqt, forwarded to screen.sh (default work/receptor.pdbqt)
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

LIG_SDF="${1:?usage: screen_consensus.sh LIGAND_DIR OUT_DIR [SEEDS]}"
OUT="${2:?usage: screen_consensus.sh LIGAND_DIR OUT_DIR [SEEDS]}"
SEEDS="${3:-42 1 2 3 4}"

echo "consensus screen: seeds=[$SEEDS] ligands=$LIG_SDF out=$OUT receptor=${RECEPTOR:-work/receptor.pdbqt}"
for s in $SEEDS; do
  echo "=== seed $s -> $OUT/s$s ==="
  SEED_OVERRIDE="$s" "$ROOT/scripts/screen.sh" "$LIG_SDF" "$OUT/s$s"
done
echo "CONSENSUS_SCREEN_COMPLETE seeds=[$SEEDS] out=$OUT"
