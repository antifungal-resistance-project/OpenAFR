# TODOS

Deferred work captured during review. Each item has enough context to pick up cold.

## Build the ERG11 re-caller (reads → azole-resistance call) — OPEN, TOP PRIORITY

**What:** the single missing piece that would make the genomic early-warning track (issues
#20–27) actually work. Given an NCBI *C. auris* isolate's linked **SRA raw reads**, call the
ERG11 (CYP51) azole-resistance substitutions — Y132F, K143R, F126L, the VF125AL/F126L
clade-III haplotype, TR34-style tandem duplications, plus TAC1B GOF / ERG3 where feasible —
and populate the `erg11_call` column that every downstream stage already consumes.

**Why it's the whole ballgame:** the spike (#20) proved NCBI runs *no* AMR pipeline on
*C. auris*, so there is no resistance call to read off the feed — it must be manufactured
from reads. Every other stage (snapshot, emergence, mapping, structural so-what, alert,
delivery) is built, tested, and — per the backtest's H3 mechanism proof — fires end to end
the moment `erg11_call` is populated. Until then the real-data backtest is an honest null
(H2: 0 flags over 8 years, because calls are absent, not because emergence is absent). This
re-caller is the entire gap between "scaffold" and "validated early warning."

**Where to start:** ~96.6% of isolates carry linked SRA `Run` accessions (assemblies only
3.4%, so work from reads, not assemblies). The snapshot schema already reserves
`erg11_call` + `resistance_source` and marks uncalled isolates `pending:sra-recaller`
(see `data/earlywarning/README.md`). The arrival budget is adequate (backtest H1: median
97-day deposit lag), so the effort is justified.

**Validation:** the backtest (#27) is pre-registered to use an *external* truth set
(published clade/mutation panels, CDC AR Isolate Bank strains with known ERG11 status), not
NCBI labels — only 23 isolates carry any AST phenotype. Wire the re-caller's calls through
`scripts/backtest_earlywarning.py replay` on a real snapshot to convert the H2 null into a
real-data result.

**Source:** spike #20 + backtest #27 (`work/RESULTS_backtest.md`), 2026-08-14. This is also
the shared blocker gating the FKS1/v2 decision below.

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

**Trigger status after #27 (2026-08-14): NOT fired — deferred behind the same
re-caller.** The backtest ([work/RESULTS_backtest.md](work/RESULTS_backtest.md))
found azole event frequency is currently **zero-*measurable***, not measured-near-zero:
`erg11_call` is empty on all real isolates because the ERG11 re-caller was never built,
so the real-data walk-forward is an honest null (H2). The redirection to FKS1/v2 is only
justified by a *measured* near-zero azole signal, which requires the re-caller first. So
the single blocker for both the azole MVP and this v2 decision is the same: build the
reads→ERG11 re-caller, then re-read the event frequency. The arrival budget is adequate
(H1: median deposit lag 97 d, PASS-arrival), so the re-caller is worth building.

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
