# TODOS

Deferred work captured during review. Each item has enough context to pick up cold.

## Script the receptor preparation (5TZ1.pdb -> receptor.pdbqt)

**What:** Add a script that turns `data/structures/5TZ1.pdb` into the docking
receptor `work/receptor.pdbqt` (chain A + heme, correct protonation/charges),
the same one the run-2 result was produced against.

**Why:** This step was done by hand and never captured. `work/receptor.pdbqt` is
gitignored, so a fresh clone does not have it and `run_screen2.sh` cannot dock —
the reproduce recipe in the README dead-ends at step 1 for any new developer.
This is the single biggest onboarding blocker.

**Context:** The gate reads the heme-iron position from the tracked
`work/receptor_A.pdb` (a stripped chain-A + heme PDB, already in the repo); the
*docking* input is the separate `work/receptor.pdbqt`. The script must reproduce
that exact `.pdbqt` — same atoms, protonation, and Gasteiger charges — or docking
scores will not match the frozen result. Verify by re-docking a known active and
checking the top score against `work/RESULTS_redock_VT1.md`.

**Source:** README reproducibility gap, 2026-08-01.

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

**Context:** The serial loop is the `for sdf in ...; do ... vina ...; done` block
in `scripts/run_screen2.sh` (lines 10–22). Determinism is preserved as long as
each ligand keeps its own fixed seed — do NOT share a seed across the pool. The
box and search parameters already come from the hash-frozen `protocol.yaml` (via
`eval "$(python openafr/protocol.py --shell docking)"` at the top of the script);
any parallel version must keep reading them from there so what is docked stays in
sync with what is graded.

When the `screen` stage is built, put the parallel dispatch inside it rather than
forking yet another one-off docking script. Three docking scripts already exist
(`run_screen.sh`, `run_screen2.sh`, `dock_batch.sh`); do not add a fourth — the
screen stage should be the one place a novel library gets docked.

**Depends on / blocked by:** The `screen` stage existing, and ideally a
wet-lab-collaborator-scoped compound library to run it against (see the
mycology-foundation design doc). Do not build before there's something at scale
to run.

**Source:** /plan-eng-review finding P1, 2026-08-01. Task ID T6.
