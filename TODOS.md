# TODOS

Deferred work captured during review. Each item has enough context to pick up cold.

## Parallelize the docking loop (screen stage)

**What:** Replace the serial `for sdf in ...; do vina ...; done` loop with a
parallel dispatch (`xargs -P $(nproc)` or GNU `parallel`), one Vina process per
core, preserving the per-ligand fixed `--seed 42`.

**Why:** Each Vina run at `--exhaustiveness 32` takes minutes. At validation
scale (~355 molecules) serial is fine. The `screen` stage points the pipeline
at a novel library (potentially thousands of compounds) — serial that is days
of wall time on one machine. Docking runs are independent, so this is
embarrassingly parallel.

**Pros:** Turns a multi-day screen into an afternoon; keeps solo iteration fast.
**Cons:** None functional (runs are independent). Minor: interleaved progress
output needs per-ligand log files instead of a shared stdout counter.

**Context:** Current serial loop is `scripts/run_screen2.sh:7`. Determinism is
preserved as long as each ligand keeps its own fixed seed — do NOT share a seed
across the pool. Fold this into the single parameterized docking script from
eng-review task T1 rather than re-forking a script.

**Depends on / blocked by:** The `screen` stage existing, and ideally a
wet-lab-collaborator-scoped compound library to run it against (see the
mycology-foundation design doc). Do not build before there's something at scale
to run.

**Source:** /plan-eng-review finding P1, 2026-08-01. Task ID T6.
