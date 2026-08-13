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

## Widen early-warning detection to echinocandin/FKS1 resistance (v2) — DEFERRED 2026-08-12

**What:** Extend the resistance early-warning track (issues #20-27) to detect
emerging FKS1/echinocandin resistance in C. auris, in addition to azole/ERG11.

**Why:** C. auris is already largely azole-resistant at baseline, so novel
azole-mutation *emergence* may be a low-frequency signal. The genuinely dynamic
and clinically scary emergence axis for C. auris is echinocandin (FKS1) resistance
and pan-resistance. The azole MVP may build a beautiful pipeline that rarely fires.

**Trigger (data-gated, not a hunch):** the azole MVP's spike (#20) and backtest
(#27) are pre-registered to report azole-mutation *event frequency*. If that comes
back near-zero, this TODO is the redirection for v2.

**Context / where to start:** detection (extract FKS1 hot-spot regions, call the
S639/F-region mutations) is tractable and mirrors the ERG11 re-caller (T1). The
blocker is the STRUCTURAL layer: there is no FKS1 pocket tooling in the repo (the
whole moat is CYP51/ERG11 azole-docking). Widening detection without a structural
so-what turns the FKS1 half into a bare novelty feed — so v2 needs FKS1 structural
tooling built from scratch, or an honest "detection-only, no structural verdict"
label for that half.

**Depends on / blocked by:** azole MVP shipped (#20-27) + backtest event-frequency
result. Do NOT start before the MVP tells you whether azole emergence fires.

**Source:** /plan-eng-review D4 (kept azole-only MVP) + outside-voice finding 1,
2026-08-12.
