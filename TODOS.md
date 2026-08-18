# TODOS

Deferred work captured during review. Each item has enough context to pick up cold.

## Build the ERG11 re-caller (reads → azole-resistance call) — CORE DONE 2026-08-15, ORCHESTRATION + REAL-DATA VALIDATION OPEN

**What:** the single missing piece that would make the genomic early-warning track (issues
#20–27) actually work. Given an NCBI *C. auris* isolate's linked **SRA raw reads**, call the
ERG11 (CYP51) azole-resistance substitutions and populate the `erg11_call` column that every
downstream stage already consumes.

**What is now built (2026-08-15):**
- `data/earlywarning/erg11_reference/` — pinned B8441 reference CDS (accession
  XM_085597798.1, CDS sha256 764e6d1a…) + PROVENANCE.md. Verified: residues 125/126/132/143
  translate to V/F/Y/K, so an emitted `Y132F` is guaranteed to mean residue 132.
- `openafr/recaller.py` — the deterministic core: consensus ERG11 CDS → substitution tokens.
  Honesty enforced + tested: asserts the reference (hash + numbering), an ambiguous codon is
  *uncalled* not wild-type, a frameshift/length mismatch is *refused* not frame-guessed, the
  known panel (V125A/F126L/Y132F/K143R) is *tagged* but non-panel substitutions are emitted
  verbatim (novelty is emergence.py's job, not the re-caller's). 16 tests, all green.
- `scripts/recall_erg11.py` — CLI. `call` = deterministic (no network/tools). `recall`/`fill`
  = the reads→consensus orchestration (fasterq-dump + minimap2 + samtools consensus), gated on
  those binaries like the docking scripts gate on vina. `resistance_source` becomes
  `sra-recaller:{called,wild-type,partial(...),failed(...),refused(...)}`.

**What remains OPEN:**
1. **Run the orchestration on a real snapshot.** The reads→consensus half needs the external
   tools installed (bioconda) + network + multi-GB SRA downloads, and is meant for a Linux/x86
   host (bioconda has no osx-arm64 builds). What is now *de-risked offline* (2026-08-18): the
   orchestration **logic** is tested in `tests/test_recall_orchestration.py` — the exact
   fetch/align/consensus pipeline `reads_to_consensus` issues (short-read preset + min-depth
   threaded into `samtools consensus`), the `_read_consensus_fasta` bridge (joins a wrapped
   record; refuses empty/header-only/multi-record), and `cmd_fill`'s batch behavior (re-calls
   only pending rows carrying a run_acc, picks the first linked run, marks in place, leaves
   resolved rows untouched, honors --limit/--dry-run). So a `fill` bug can no longer silently
   corrupt a snapshot *after* the expensive downloads. What is still UNPROVEN and needs the
   real run: the external tool **binaries** themselves (are they invoked correctly end-to-end,
   do real reads assemble a callable consensus). `fill --limit N` is the entry point; start
   small, sanity-check via `scripts/recaller_sanity.py` against known-genotype strains, widen.
2. **Convert the backtest null into a real-data result.** The backtest (#27) is pre-registered
   to use an *external* truth set (published clade/mutation panels, CDC AR Isolate Bank strains
   with known ERG11 status), not NCBI labels — only 23 isolates carry any AST phenotype. Run a
   re-called snapshot through `scripts/backtest_earlywarning.py replay` to turn the H2 null
   (0 flags because calls were *absent*) into a measured azole event frequency.
3. **Out of scope for this CDS-level caller** (documented, not silently dropped): TR-type
   tandem-duplication/promoter calls and non-ERG11 loci (TAC1B GOF, ERG3) — none are point
   substitutions in the ERG11 CDS. Revisit if the measured azole signal warrants it.

**Source:** spike #20 + backtest #27 (`work/RESULTS_backtest.md`), core built 2026-08-15. Step
2's measured event frequency is also the shared blocker gating the FKS1/v2 decision below.

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

**Trigger status: NOT fired — still deferred behind a *measured* azole signal.** The
backtest ([work/RESULTS_backtest.md](work/RESULTS_backtest.md)) found azole event
frequency was **zero-*measurable***, not measured-near-zero, because `erg11_call` was
empty (no re-caller). As of 2026-08-15 the re-caller's **deterministic core is built and
tested**, but it has not yet been *run* on a real snapshot (the reads→consensus
orchestration needs the bioinformatics tools + SRA downloads; see the top TODO's OPEN
step 1). So the blocker for both the azole MVP and this v2 decision is no longer "build
the re-caller" but "**run** it and re-read the event frequency." The redirection to
FKS1/v2 stays justified only by a *measured* near-zero azole signal. The arrival budget is
adequate (H1: median deposit lag 97 d, PASS-arrival).

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
