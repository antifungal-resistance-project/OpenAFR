# Pre-registration — validate the FKS1 re-caller against the full 100-isolate external benchmark

> **STATUS: FROZEN 2026-09-07, before any isolate beyond the 4-strain sanity control was
> scored.** This document fixes the metric, the pass bar, and the interpretation of *both*
> outcomes in advance. The grader (`scripts/validate_fks1_concordance.py`) re-checks this
> document's SHA-256, the pinned reference FKS1 CDS hash, and the truth-set fixture's hash,
> and refuses to certify a verdict if any moved. No concordance number over the full
> benchmark has been read at freeze time.

## Why this run is required — the caller's error rate is asserted, not measured

The FKS1 echinocandin re-caller (`openafr/fks1_caller.py`, orchestration
`scripts/recall_fks1.py`) has produced a **measured prevalence** on a real snapshot (2.3%
[1.2%, 4.1%], `work/RESULTS_fks1_prevalence.md`). But the caller's own **accuracy** has only
ever been checked as a 4-strain PASS/FAIL *control* (`scripts/recaller_sanity_fks1.py`,
`data/earlywarning/recaller_sanity/fks1_known_genotypes.tsv`): B11220 wild-type + one
positive each for S639F/S639P/S639Y. Four strains is a smoke test, not an error rate — the
exact gap the "widen n=7 → N=300" active-power run (`work/PREREGISTRATION_active_power.md`,
sha `9327563a…`) closed for the docking gate. This run does the analogous thing for the
caller: score it against a **fixed, published, purpose-built external truth set** and report
sensitivity, specificity and per-token identity with confidence intervals.

## The truth set (external, pinned)

**PMC12323592** — *"A benchmark dataset for validating FKS1 mutations in Candida auris"*
(Microbiol Spectr; doi 10.1128/spectrum.03147-24). It exists precisely to validate FKS1
variant callers. Its Table 1 lists **100 WGS isolates** with, per isolate, the SRA run
accession, the FKS1 genotype, and the echinocandin susceptibility interpretation. Reported
composition (from the paper, verified 2026-09-07): **47 resistant / 53 susceptible**; HS1
panel tokens **S639F ×15, S639Y ×13, S639P ×8**, plus non-HS1 changes the caller does not
tag (F635Y/F635C ×1 each, D642Y ×7, HS2 R1354S, HS3 W691L/W691C). The four strains already
in `fks1_known_genotypes.tsv` are a subset of these 100 and must keep the same labels.

**Fixture + hash discipline (load-bearing).** The 100 rows are transcribed from the paper
into `data/earlywarning/recaller_sanity/fks1_benchmark_100.tsv` (same column schema as the
sanity fixture: `strain clade run_acc expected_panel susceptibility note citation`). The
transcription is a **verified manual step** — the paper is not open-access, so it cannot be
scraped by a summarising model (an LLM scrape of the table produced a self-contradictory
duplicate row, exactly the failure mode this discipline exists to prevent). Every SRA run
accession is confirmed live against the ENA `read_run` API (Illumina/PAIRED/WGS,
`scientific_name` = *Candida auris*) before the fixture is frozen; the fixture's SHA-256 is
pinned in `work/PREREG_fks1_concordance.sha256` and re-checked by the grader. **This run does
not certify until that fixture exists and its hash is pinned.** A non-Illumina run (misaligns
under `-ax sr`) is excluded and noted, never scored.

## The metric (fixed now) — `openafr.concordance.evaluate`

The per-isolate evaluation runs the **identical** orchestration as the real `fill`
(`recall_fks1.reads_to_window_consensus` → `fks1_caller.call_windows`), so a pass here
vouches for the production path. Each isolate yields `(expected panel token set, called
panel token set, resolved?, resistant?)`; `openafr.concordance.evaluate` (pinned, unit-tested
in `tests/test_concordance.py`) computes:

1. **Panel DETECTION** — expected-panel-positive vs any panel call → **sensitivity**
   (TP/(TP+FN)) and **specificity** (TN/(TN+FP)), Wilson 95% CIs.
2. **Panel IDENTITY** — among detected positives, the fraction where the *exact* expected
   token was emitted (S639F ≠ S639P).
3. **Per-token sensitivity** — per expected-carrier, one record each for S639F/S639P/S639Y.
4. **Resistance DETECTION** — of resolved *phenotypically resistant* isolates, the fraction
   the HS1 caller flags. Reported as a **measured ceiling, not gated** (see below).

**Excluded, not guessed.** An isolate whose panel window is uncalled (coverage gap), refused
(in-window indel) or fails to orchestrate is `resolved=False`: excluded from every
denominator and reported as `n_unresolved`. This is the #22 undated discipline / the
`panel_prevalence` "excluded not counted negative" rule.

## Pre-committed pass bar (fixed now, both directions)

Computed over **resolved** isolates. The costly error for an early-warning system is a
**false emergence alert** (calling S639* on a wild-type isolate), so specificity is gated
hardest.

    PASS   if  specificity      point >= 0.95  AND  Wilson lower bound >= 0.90
           AND panel sensitivity point >= 0.90  AND  Wilson lower bound >= 0.75
           AND identity         point >= 0.95
           AND resolved fraction (n_resolved / 100) >= 0.80
    FAIL   if any gated metric's point estimate is below its bar
    UNDERPOWERED  if resolved fraction < 0.80 (coverage too thin to certify) — reported as
                  neither a pass nor a fail, with the achieved n and the fixes to try
                  (more reads / lower min-depth), mirroring the INCONCLUSIVE control outcome

The sensitivity lower bound is set at 0.75 (not 0.90) deliberately: with ~36 HS1-positive
isolates the interval is intrinsically wide, and the pre-commitment must be honest about the
power the benchmark actually provides rather than demand a bound the n cannot supply.

## Interpretation, committed now (both directions)

- **PASS.** The caller emits the right panel token on real reads at a measured, CI-bounded
  rate, and essentially never invents one on a wild-type isolate. The prevalence number and
  every emergence alert built on the caller inherit that measured error rate; the
  detection-only FKS1 track rests on a validated caller, not an asserted one. Folded into
  `work/RESULTS_fks1_prevalence.md` and the preprint/track write-up as the caller's accuracy.

- **FAIL.** The caller's real-reads accuracy does not meet the bar. This is a genuine,
  publishable negative that **blocks** promoting the FKS1 track past detection-of-record: a
  named, reproducible defect (which token, which direction) with the benchmark isolates that
  expose it. No prevalence or emergence claim may quote the caller as validated until fixed.

- **The resistance-detection number is NOT a pass/fail.** It is expected to sit **below**
  panel sensitivity because the benchmark contains resistant isolates whose mechanism is
  non-HS1 (F635/D642/R1354/W691) or non-*FKS1*, which an HS1 caller cannot and must not call.
  Reporting it quantifies, honestly, the ceiling on what a detection-only HS1 pipeline can
  catch — the same honesty as labelling the whole FKS1 half "detection-only, no structural
  verdict."

## What is NOT in scope (declared, not silently dropped)

- **TR/promoter and non-HS1 loci.** The caller is CDS point-substitution only; HS2 R1354 is
  reported-but-untagged, HS3 W691 and non-*FKS1* mechanisms are out of scope. They lower the
  resistance-detection ceiling (item 4) and that is reported, not hidden.
- **MIC regression.** The benchmark's MIC columns are not modelled; the phenotype is used
  only as the binary R/S anchor for item 4.
- **ERG11.** The identical harness applies to the azole caller against Lockhart 2017, but
  that truth set is 4 clade strains, not a 100-isolate benchmark; this pre-registration is
  the FKS1 run. An ERG11 concordance run is a separate (smaller-n) pre-registration.
