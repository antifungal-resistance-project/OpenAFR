#!/bin/bash
# One-command runner for the multi-baseline comparison (#83).
# Implements work/PREREGISTRATION_baselines_multi.md end to end on a capable host.
#
# Does, in order:
#   0. PREFLIGHT   — scripts/preflight_baselines.py (fail fast on frozen-hash DRIFT;
#                    otherwise report what is still outstanding). Aborts on drift (exit 1).
#   1. GNINA       — scripts/run_gnina_rescore.sh   (skipped-with-note if gnina absent)
#   2. PLANTS      — scripts/run_plants_rescore.sh  (skipped-with-note if PLANTS absent)
#   3. GRADE       — scripts/rescore_baselines.py   (always; prints the table or the
#                    honest 'not-yet-run (cloud-gated)' state for any missing baseline)
#
# The two rescore scripts self-gate on their Linux/x86 binaries and exit 2 when absent;
# this wrapper treats that as "skip, not fatal" so the grader still runs and shows the
# current picture. The whole thing therefore also runs on a host with NO scorers — it
# just prints mode C plus the pending rows, exactly like calling the grader directly.
#
# Cloud recipe (pinned tool versions, receptor prep, pose regeneration): see
# work/RUNBOOK_baselines_multi.md. On a fresh host run pose regeneration + receptor prep
# from the RUNBOOK FIRST; this script assumes work/screen2 + the receptors exist.
#
# Usage: scripts/run_baselines_all.sh [POSE_DIR] [GRADE_RECEPTOR]
#   POSE_DIR        Vina poses to rescore              (default work/screen2)
#   GRADE_RECEPTOR  structure to read the heme Fe from (default work/receptor_A.pdb)
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

POSE_DIR="${1:-work/screen2}"
GRADE_RECEPTOR="${2:-work/receptor_A.pdb}"
PY="${PYTHON:-python3}"

echo "### [0/3] preflight"
if ! "$PY" scripts/preflight_baselines.py "$POSE_DIR"; then
  code=$?
  # 1 = frozen-hash DRIFT: refuse to run. 2 = run-steps-remain: expected, proceed.
  if [ "$code" -eq 1 ]; then
    echo "ABORT: preflight reported DRIFT — resolve the frozen-hash mismatch first." >&2
    exit 1
  fi
fi

echo; echo "### [1/3] GNINA rescore"
scripts/run_gnina_rescore.sh "$POSE_DIR" work/gnina_scores work/receptor.pdbqt \
  || echo "  (gnina step skipped — see message above; grader will mark it not-yet-run)"

echo; echo "### [2/3] PLANTS rescore"
scripts/run_plants_rescore.sh "$POSE_DIR" work/plants_scores work/receptor_plants.mol2 \
  || echo "  (plants step skipped — see message above; grader will mark it not-yet-run)"

echo; echo "### [3/3] grade"
"$PY" scripts/rescore_baselines.py "$POSE_DIR" "$GRADE_RECEPTOR"
grade=$?

echo
if [ "$grade" -eq 0 ]; then
  echo "COMPARISON COMPLETE — record the table in work/RESULTS_baselines_multi.md and"
  echo "interpret per work/PREREGISTRATION_baselines_multi.md (CI overlap = tie)."
else
  echo "COMPARISON INCOMPLETE — one or more baselines not yet run (Linux/x86 binaries)."
  echo "Fill them per work/RUNBOOK_baselines_multi.md, then re-run this script."
fi
exit "$grade"
