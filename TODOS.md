# TODOS

Deferred work captured during review. Each item has enough context to pick up cold.

## Rank WITHIN the azole class — CLOSED (for good) 2026-08-26 by the independent confirmation (#9 FAIL)

> **UPDATE 2026-08-26 — look #9 (independent-active confirmation) FAILED. The within-azole line is
> now closed on BOTH the rigid and the exhaustive-ensemble receptor. Do NOT re-open without a
> genuinely new axis (not another conformer set).**
>
> Look #8's ensemble PASS (AUC 0.750) was on only 7 held-out azoles with a wide CI (0.683–0.825),
> and its own PASS branch demanded a pre-registered confirmation on an INDEPENDENT active set
> before any within-class claim. That confirmation ran: an independent N=200 (seed 42) sample of
> 3,609 ChEMBL measured-active azoles (novelty-filtered <0.70 vs all 15 actives + the 5TZ1
> co-crystal), same frozen ensemble/criterion/inactives — **only the actives changed**. Result:
> **ensemble-min mode C AUC 0.682 < 0.70** (tight CI 0.647–0.715 at n=200; permutation p<0.0001,
> so the signal is REAL but sub-threshold). The lift is genuine but small — 5TZ1 alone 0.638 →
> ensemble 0.682 — and which conformer helps is itself sample-dependent (5FSA carried look #8;
> 5V5Z is best here at 0.701). Novel-chemotype guard still holds at 0.823 (comparator 0.810), so
> the PRODUCT is intact — the failure is specifically within-class ranking.
> Pre-reg `work/PREREGISTRATION_ensemble_confirm.md` (frozen sha `7754050e…`); result
> `work/RESULTS_ensemble_confirm.md`; grader `scripts/validate_gate_ensemble_confirm.py`; active-set
> builder `scripts/make_verified_actives.py`. **Settled conclusion:** pose geometry over a fixed
> *C. albicans* CYP51 target (rigid OR full crystallographic ensemble) cannot rank working azoles
> above failed azole analogues to the 0.70 bar — a property of the QUESTION, not of any one
> feature or conformer. The defensible product remains novel-chemotype triage (AUC 0.810–0.823).

> **SUPERSEDED 2026-08-25 — look #8 (ensemble receptor) PASSED at AUC 0.750 and briefly re-opened
> this line, but look #9 above shows that PASS did not replicate on an independent set (n=7 → n=200).
> Look #8 stands as graded (`work/RESULTS_ensemble.md`) but is superseded in interpretation.**

> **SUPERSEDED 2026-08-25 — look #7 (mode F, channel engagement) FAILED, closing the line on the
> RIGID conformer. The ensemble looks #8/#9 extended the test to the conformer ensemble; #9 closes it.**
>
> Mode F AUC 0.590 on the matched 7-active + azole-bearing-inactive pool (mode C comparator:
> 0.650). An orthogonal, mechanistically-motivated tail-fit feature also cannot separate
> working from failed azoles on a single rigid conformer. Size-residualized AUC 0.612 —
> the failure is not explained by mode F re-measuring molecular weight. Result:
> `work/RESULTS_channel_engagement.md`.
>
> **Per pre-registration, this is the decisive evidence for Option B: re-scope the product
> to novel-chemotype triage (non-azole pool, AUC 0.810, EF@5% 4.98x).** The within-azole
> ceiling is a property of the rigid pose, not of any specific feature. Any further within-class
> attack requires ensemble/flexible-receptor poses (new dataset, new pre-registration).
> No further feature attempts on `work/screen_verified/` are permitted.

**History.** Look #6 (`work/RESULTS_verified_inactives.md`) passed overall (AUC 0.716) but
located the ceiling: AUC 0.810 vs non-azole chemotypes, AUC 0.650 vs azole-bearing inactives.
Mode C (iron-approach distance) detects chemotype, not potency. Look #7 tested mode F
(channel engagement) as the one permitted orthogonal feature — it also failed (0.590 < 0.70),
confirming the rigid-conformer limitation rather than any feature-specific weakness.

**Product re-scope (effective immediately):** the tool's defensible claim is
**novel-chemotype triage at AUC 0.810**. Shortlists should be filtered to exclude azole-class
molecules; within-class azole ranking is out of scope until ensemble/flexible poses exist.

**Source:** look #6, 2026-08-23. Depends on nothing; the dataset and poses reproduce from
`work/RESULTS_verified_inactives.md`.

## Build the ERG11 re-caller (reads → azole-resistance call) — DONE + RUN 2026-08-20 (representative fill measured)

> **UPDATE 2026-08-20 — the deliverable is measured.** The re-caller was run on a real,
> representative sample: `fill --random --seed 0 --limit 200` on snapshot PDG000000067.675
> (199 called / 1 partial / 0 failed; runlog `fill.jsonl` run_id `e4a0def6586a`; GCP VM
> since deleted). **Azole event frequency = 160/199 = 80.4%, Wilson 95% CI [74.3%, 85.3%]**
> ([work/RESULTS_prevalence.md](work/RESULTS_prevalence.md)). OPEN steps 1 and 2 below are
> now DONE — the orchestration binaries are proven end-to-end and the H2 null is converted
> to a measured frequency. The FKS1/v2 data gate below has **FIRED** (~80% baseline = azole
> resistance is near-saturated, so azole *emergence* is a weak early-warning signal).


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
   do real reads assemble a callable consensus). Run order for the real host: (0) pull a
   fresh snapshot (`snapshot_ncbi_auris.py pull` — verified live 2026-08-18, latest release
   PDG000000067.673, 29,153 isolates / 96.6% SRA-callable); (1) `recall_erg11.py plan
   <snapshot>` — offline preflight (no tools/network) that reports how many rows a `fill`
   would download/re-call and lists the accessions, so the multi-GB scope is known before
   committing; (2) `scripts/recaller_sanity.py --only B8441` — the wild-type positive
   control (a spurious panel hit there means the orchestration, not the science, is wrong),
   then widen the sanity set; (3) `fill --limit N`, start small and widen. Note: the tool
   gate now also **version-checks** samtools (needs `consensus`, >=1.13) before any download,
   so an old-samtools box fails instantly instead of after the reads are pulled.
2. **Convert the backtest null into a real-data result.** The backtest (#27) is pre-registered
   to use an *external* truth set (published clade/mutation panels, CDC AR Isolate Bank strains
   with known ERG11 status), not NCBI labels — only 23 isolates carry any AST phenotype. Run a
   re-called snapshot through `scripts/backtest_earlywarning.py replay` to turn the H2 null
   (0 flags because calls were *absent*) into a measured azole event frequency. The deliverable
   number itself now has a command: `backtest_earlywarning.py prevalence <snapshot>` reports
   the azole event frequency (panel-positive / *resolved* panel, with partial/failed/pending
   excluded not counted azole-negative) — undefined ("NOT MEASURED YET") until `fill` populates
   calls, then the measured value the FKS1/v2 decision below gates on.
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

## Widen early-warning detection to echinocandin/FKS1 resistance (v2) — BUILT (detection-only), 2026-09-01

> **UPDATE 2026-09-01 — the v2 detection-only pipeline is now built end-to-end.** The
> FKS1/echinocandin track runs the full loop *in parallel to* ERG11, honestly labelled
> detection-only (no structural verdict — echinocandins coordinate no metal, so the
> CYP51/heme-iron moat does not transfer). Shipped increments:
> - **T1 (#51)** — `openafr/fks1_caller.py`: windowed reads→consensus→S639-panel call,
>   pinned/numbering-verified GSC1 reference (`data/earlywarning/fks1_reference/`).
> - **T2 (#52)** — `scripts/recall_fks1.py` `call`/`recall`: the CLI, tool-gated like ERG11.
> - **T3 (#111)** — snapshot schema `fks1_call`/`fks1_resistance_source` + migration +
>   batch `plan`/`fill`/`prevalence`; `detect_emergence(call_field=...)` serves FKS1.
> - **Alert surface (#112)** — `alert.compose_fks1_alerts` + `render_fks1_markdown`;
>   `compose_alert.py --gene fks1` (WATCH/CONTEXT only, never ACT-NOW).
> - **Delivery (this change)** — `deliver_alert.py --gene fks1` with its OWN committed
>   digest/state (`DIGEST_FKS1.md`/`STATE_FKS1.json`); the scheduled workflow delivers
>   both genes independently. `openafr/delivery.py` routes off the result's `gene` marker.
>
> **What remains OPEN for v2:**
> 1. **A real FKS1 `fill` → measured prevalence/backtest.** `fks1_panel_prevalence` reports
>    "NOT MEASURED YET" until an actual re-caller `fill` runs — same wall as ERG11: bioconda
>    tools (fasterq-dump/minimap2/samtools≥1.13) + multi-GB SRA downloads on a Linux/x86
>    host (offline-tested orchestration; unproven binaries). Not doable on osx-arm64.
> 2. **T5 — the FKS1 structural so-what: still DEFERRED BY DESIGN**, not merely unbuilt.
>    Echinocandins coordinate no metal; the CYP51 geometry moat does not transfer to a
>    glucan synthase. Revisit only if a PI asks for it. The detection-only label is the
>    honest scope, and it is now enforced through composition, rendering, and delivery.

**What:** Extend the resistance early-warning track (issues #20-27) to detect
emerging FKS1/echinocandin resistance in C. auris, in addition to azole/ERG11.

**Why:** C. auris is already largely azole-resistant at baseline, so novel
azole-mutation *emergence* may be a low-frequency signal. The genuinely dynamic
and clinically scary emergence axis for C. auris is echinocandin (FKS1) resistance
and pan-resistance. The azole MVP may build a beautiful pipeline that rarely fires.

**Trigger (data-gated, not a hunch):** the azole MVP's spike (#20) and backtest
(#27) are pre-registered to report azole-mutation *event frequency*. If that comes
back near-zero, this TODO is the redirection for v2.

**Trigger status: FIRED (2026-08-20) — via baseline saturation, measured.** The re-caller
has now been run on a representative random sample (n=199 resolved): **azole event
frequency = 160/199 = 80.4%, Wilson 95% CI [74.3%, 85.3%]**
([work/RESULTS_prevalence.md](work/RESULTS_prevalence.md)). Read the nuance honestly: the
measured number is **high, not near-zero** — but the redirection was never really about a
near-zero *prevalence*; it was about a near-zero *emergence signal-to-noise*, which an ~80%
saturated baseline produces. With ~4 of 5 isolates already panel-positive there is little
headroom for azole *emergence* to be a useful early-warning trigger — the pipeline would
fire on a mutation that is already the norm. So the redirection premise ("C. auris is
already largely azole-resistant at baseline") is now **confirmed by data**, and v2 (FKS1)
is the justified next direction. The arrival budget is adequate (H1: median deposit lag
97 d, PASS-arrival).

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
