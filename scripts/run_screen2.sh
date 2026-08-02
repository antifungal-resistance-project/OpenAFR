#!/bin/bash
# Run 2 — the pre-registered held-out validation docking.
#
# Thin wrapper over the one docking engine (scripts/screen.sh): it docks the
# run-2 ligands (work/run2_ligands) into work/screen2 under the frozen protocol,
# exactly as before. Validation and any future novel-library screen therefore run
# through identical code — there is no separate, drift-prone docking loop here.
#
# Protocol FIXED BY work/PREREGISTRATION_run2.md. Do not tune. The box + search
# parameters are read from the hash-frozen protocol.yaml inside screen.sh.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
exec "$ROOT/scripts/screen.sh" work/run2_ligands work/screen2 "$@"
