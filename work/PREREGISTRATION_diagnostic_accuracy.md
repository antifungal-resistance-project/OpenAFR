# Pre-registration — retrospective verdict accuracy vs known-MIC isolates (issue #137)

> **STATUS: FROZEN 2026-09-10, before any verdict has been scored against any phenotype.**
> This document fixes, in advance, the verdict→call mapping, the two clinically loaded error
> metrics, the pass bar, and the interpretation of *both* outcomes — *and* records that the
> run it governs is **BLOCKED on data**, not runnable today (see "The data this run needs").
> The grader (`scripts/validate_diagnostic_accuracy.py`) re-checks this document's SHA-256,
> the pinned caller reference hashes, and the truth-set fixture's hash, and refuses to certify
> if any moved or if the fixture is absent. Freezing the metric before the fixture is in hand
> removes the freedom to pick the error bar after seeing which isolates we could obtain. The
> SHA-256 is pinned in `work/PREREG_diagnostic_accuracy.sha256`.

## Why this run is required — the verdict-accuracy claim is asserted, not measured

The diagnostics MVP ships a per-drug-class **verdict** (`openafr/verdict.py`,
`docs/DIAGNOSTIC_VERDICT_CONTRACT.md`, #133): a four-value categorical call —
`RESISTANCE_MARKER_DETECTED` / `NO_KNOWN_MARKER` / `UNCHARACTERIZED_VARIANT` / `UNRESOLVED`.
Two credibility checks already exist and neither is this one:

- The **caller** concordance run (`work/PREREGISTRATION_fks1_concordance.md`, sha `da8e3827…`;
  `openafr/concordance.py`) measures TOKEN accuracy against a *genotype* truth set — "did the
  caller emit S639F when S639F was present?"
- The **calibration** run (`work/PREREGISTRATION_calibration.md`, #136) would measure whether
  a stamped probability is reliable — but the contract ships **no probability field**, so that
  track is deferred.

Neither answers the question a clinician actually asks of a resistance test: *when the engine
emits a verdict, how often does it agree with the isolate's phenotype, and how often does it
get the two dangerous directions wrong?* That is this run. It anchors on **phenotypic AST
truth** (the R/S interpretation), not genotype, and reports the CLSI-M23 error names.

## The metric (fixed now) — `openafr.accuracy.evaluate`

Per isolate, per drug class, the engine's verdict is mapped to an R/S CALL by the fixed rule
below, and scored against the phenotype. `openafr.accuracy.evaluate` (pinned, unit-tested in
`tests/test_accuracy.py`) computes, **per drug class** (azole VME and echinocandin VME are
different clinical facts and are never pooled):

**Verdict → call mapping (frozen):**

| verdict | call | in 2×2? |
|---|---|---|
| `RESISTANCE_MARKER_DETECTED` | positive ("R") | yes |
| `NO_KNOWN_MARKER` | negative ("not-R") | yes — the assay's negative prediction, **not** a susceptibility *claim*; the VME rate is exactly the honest price of reading it as one |
| `UNCHARACTERIZED_VARIANT` | — | **no** — abstention, reported as an abstention rate (contract rule 3, first-class "I don't know") |
| `UNRESOLVED` | — | **no** — excluded, reported as `n_unresolved`, never a negative (rule 4 / "excluded, not guessed") |

An isolate with no phenotype truth is `n_no_phenotype` and excluded. Abstentions and
unresolved verdicts are held OUT of the confusion matrix, so a low VME can never be bought by
quietly reclassifying hard isolates as negatives — the same discipline as
`concordance.evaluate` and `backtest.panel_prevalence`.

**Reported metrics (each a k/n with a Wilson 95% CI):**

1. **VERY MAJOR ERROR (VME)** — phenotype R, verdict `NO_KNOWN_MARKER`. Denominated over the
   **resistant** isolates (CLSI convention). The load-bearing error: the assay missed real
   resistance — precisely the failure the contract's rule 2 warns about.
2. **MAJOR ERROR (ME)** — phenotype S, verdict `RESISTANCE_MARKER_DETECTED`. Denominated over
   the **susceptible** isolates. A false resistance call.
3. **Sensitivity** (= 1 − VME) and **specificity** (= 1 − ME), same numbers, caller-facing framing.
4. **Categorical agreement** — concordant scored isolates / scored isolates.
5. **Abstention rate** — `UNCHARACTERIZED_VARIANT` over the characterised isolates (scored +
   abstained). Reported, never hidden: it bounds coverage honestly.

## The data this run needs does not exist yet (the honest block, committed now)

Scoring a verdict against phenotype needs, per isolate, **paired genotype + phenotypic AST**
for the target gene/drug. Two purpose-built external truth sets are named, and both are
currently blocked on the SAME upstream steps as their caller runs:

- **echinocandin / FKS1** — **PMC12323592**, the 100-isolate *C. auris* FKS1 benchmark
  (genotype + echinocandin R/S). This is the fixture the FKS1 concordance run already targets
  (`data/earlywarning/recaller_sanity/fks1_benchmark_100.tsv`), **blocked on a verified manual
  (non-scraped) transcription of the paywalled Table 1** plus a GCP run. This accuracy run
  reuses that exact fixture — the isolates carry the R/S column this run scores against.
- **azole / ERG11** — Lockhart 2017 provides only **4 clade reference strains** with paired
  data; that is a sanity control, not an error rate. A powered azole arm needs a larger paired
  *C. auris* ERG11 collection (or the pooled *C. albicans* option-C route named in #136), and
  is declared **separately underpowered** until such a fixture is transcribed and pinned.

Until a fixture exists and its hash is pinned in `work/PREREG_diagnostic_accuracy.sha256`, the
grader **refuses to certify** — exactly as `validate_fks1_concordance.py` and
`validate_calibration.py` refuse without their fixtures. The scoring math it would apply
(`openafr.accuracy`) IS unit-tested offline now.

**Expected fixture schema** (TSV, `#`-comments allowed) once a paired panel is frozen — one
row per resolvable isolate:

    isolate  species  run_acc  gene  drug_class  called_verdict  susceptibility  note  citation

where `called_verdict` is the enum the production path (`recall*` → caller → `build_verdict`)
emits for that isolate — filled by the orchestration, not hand-authored — and `susceptibility`
is the phenotype's `R`/`S` against a STATED CLSI/EUCAST breakpoint.

## Pre-committed pass bar (fixed now, both directions)

Computed per drug class over **scored** isolates. VME is gated hardest — a missed resistance
is the clinically costly error for a resistance test, and it is the direction the "concordance
≠ susceptibility" honesty of the contract exists to keep visible.

    PASS         if  VME point          <= 0.03   AND  VME Wilson UPPER bound <= 0.15
                 AND ME  point          <= 0.05
                 AND abstention         <= 0.30
    FAIL         if any gated metric's point estimate is above its bar
    UNDERPOWERED if the resistant arm < 15 scored OR the susceptible arm < 15 scored
                 (too thin to read an error rate) — reported as neither pass nor fail, with the
                 achieved n and the fix (more paired isolates), mirroring the calibration and
                 concordance INCONCLUSIVE outcomes.

The VME UPPER-bound gate (0.15, well above the 0.03 point bar) is deliberate: with the ~36–47
resistant isolates a benchmark of this size supplies, a 0/n VME still admits a double-digit
true rate, and the pre-commitment must be honest about the power the truth set actually
provides rather than demand a bound the n cannot support. The azole arm is expected to trip
UNDERPOWERED at MVP and that is a valid, reported outcome, not a failure.

## Interpretation, committed now (both directions)

- **PASS.** The engine's verdict agrees with phenotype at a measured, CI-bounded rate, and its
  very-major-error rate — the rate at which it calls `NO_KNOWN_MARKER` on a truly resistant
  isolate — is bounded low. The diagnostics MVP's day-one claim ("detects named resistance
  markers, with a measured miss rate") rests on a validated verdict, not an asserted one.
  Folded into the go/no-go doc (#139) and the track write-up as the engine's accuracy.

- **FAIL.** The verdict's real accuracy does not meet the bar — a genuine, publishable
  negative that **blocks** any accuracy claim for the MVP until fixed: a named, reproducible
  defect (which drug class, which direction, which isolates) with the truth-set isolates that
  expose it. An honest negative here is a valid outcome, per project ethos.

- **UNDERPOWERED (expected for azole at MVP).** The truth set is too thin to read an error
  rate. Reported with the achieved n and the data that would unblock it; never dressed up as a
  pass. This is the honest state of the azole arm until a larger paired ERG11 panel exists.

## What is NOT in scope (declared, not silently dropped)

- **The probability field.** This run scores the categorical verdict only; probability
  reliability is the separate, deferred #136 calibration track.
- **`NO_KNOWN_MARKER` as a positive susceptibility claim.** We score it as the assay's negative
  prediction to MEASURE the VME cost of doing so; the contract still forbids reporting it as a
  clinical "S". Measuring the error is not the same as endorsing the reading.
- **Non-target mechanisms.** Efflux (TAC1/MRR1/CDR1), ERG3 and non-*ERG11/FKS1* resistance are
  out of the callers' scope; resistant isolates carrying only those mechanisms drive part of
  the VME and that is reported as a mechanism-coverage ceiling, honestly, not hidden.
- **MIC regression.** The benchmark MIC columns are not modelled; the phenotype is used only as
  the binary R/S anchor.
