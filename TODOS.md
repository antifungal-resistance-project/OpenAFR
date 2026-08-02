# TODOS

Deferred work captured during review. Each item has enough context to pick up cold.

## Parallelize the docking loop (screen stage) — DONE 2026-08-01

**Status:** Implemented. `scripts/screen.sh` is now the single docking engine and
dispatches one Vina process per core with `xargs -P` (core count via `nproc` on
Linux / `sysctl -n hw.ncpu` on Apple Silicon). Each ligand keeps the fixed
protocol `--seed`, applied independently — no shared RNG across the pool — so a
parallel run is bit-for-bit identical to a serial one. Box + search parameters are
still read from the hash-frozen `protocol.yaml` via
`eval "$(python openafr/protocol.py --shell docking)"`, so docking stays in sync
with grading. Progress is per-ligand (`OK`/`SKIP`/`*_FAIL` lines + per-ligand
`work/screen/<name>.log`); the run is resumable.

The three former docking scripts were consolidated: `run_screen.sh` and
`dock_batch.sh` were removed, and `run_screen2.sh` is now a thin wrapper over
`screen.sh` so validation and any novel-library screen share one code path.

**Remaining (still blocked):** running the screen at real scale needs a
wet-lab-collaborator-scoped compound library. Do not chase further throughput
tuning before there is something at scale to run. Screening is also only justified
after `validate_gate2.py` passes.

**Source:** /plan-eng-review finding P1, 2026-08-01. Task ID T6.
