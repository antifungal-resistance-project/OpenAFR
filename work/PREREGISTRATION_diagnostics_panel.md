# Pre-registration — the calibration-grade panel gate for the diagnostics MVP (issue #132)

> **STATUS: FROZEN 2026-09-09, before any candidate panel's composition has been read.**
> This document fixes, in advance, what "calibration-grade" means numerically, the rule that
> picks the first target, and the interpretation of *both* outcomes (a panel is found / no
> panel qualifies). No candidate collection's size, drug coverage, or variant makeup has been
> looked up at freeze time. The point of freezing a *spike* is the same as freezing a metric
> run: it removes the freedom to set the bar after seeing what is available and then rationalize
> a marginal collection into "good enough." The SHA-256 is pinned in
> `work/PREREG_diagnostics_panel.sha256`.

## Why this gate is required — "calibrated probability" is the product, and it is unbacked

The diagnostics MVP (`.context/DIAGNOSTICS_CONCEPT.md`, issues #132–139) claims to emit a
**calibrated resistance probability** with a stated, testable reliability — not just a hard
S/R label. That claim has no data behind it yet. Calibrating a probability (fitting *and*
holding out a test split, then measuring a reliability curve, #136/#137) needs materially more
than the concordance check we already have: the FKS1 PMC12323592 100-isolate benchmark
(`work/PREREGISTRATION_fks1_concordance.md`, sha `da8e3827…`) is a **checks-the-caller-is-right**
set, not a **teaches-a-probability** set. If no obtainable collection is large or diverse enough
to calibrate, then the honest MVP is not "calibrated probability" — it is concordance-only, and
we must say so *before* building the contract in #133, not discover it at #136.

So this gate runs first and is a real kill gate.

## The target-choice rule (committed now, before looking)

First target defaults to **ERG11 / azole / *Candida*** — the re-callers already live there
([[openafr-resistance-early-warning]]), so #134 is a wrap, not a new caller.

- Search the ERG11/*Candida* panel first.
- **Flip to *Aspergillus* CYP51A / azole only if** ERG11/*Candida* FAILS the panel bar below
  **and** a CYP51A collection PASSES it. A flip pays a build cost (no CYP51A caller yet), so it
  is justified only by data availability, never by preference.
- If *both* qualify, stay on ERG11 (lower build cost). If *neither* qualifies, this is a FAIL
  (see fallbacks).

## What "calibration-grade" means — fixed acceptance criteria (all required)

A candidate collection qualifies **only if it meets every one** of these. Numbers are
pre-committed judgment calls; the rationale is stated so the bar is auditable, not so it can be
renegotiated after seeing a near-miss.

1. **Paired genotype + phenotype, per isolate.** Each isolate carries (a) the resolvable
   ERG11 (or CYP51A) sequence or variant call, and (b) a phenotype: quantitative **MIC**
   preferred, or a breakpoint-classified **S/I/R** against a *stated* CLSI or EUCAST breakpoint.
   The binary R/S is the calibration target; MIC, when present, is a bonus continuous target.
2. **Size — enough to split.** **N ≥ 150** resolvable isolates, with **both classes present at
   ≥ 40 each** (resistant and susceptible). Rationale: a ~70/30 fit/test split must leave a test
   set with ≥ ~15 events per class, the floor below which a reliability curve is noise.
3. **Variant diversity — enough to calibrate the two regimes separately.** **≥ 5 distinct
   resistance-associated variants**, each with **≥ 5 carriers.** #136 must calibrate the
   *known-variant* and *novel-variant* regimes separately (they will not share a curve); a panel
   dominated by one hotspot (e.g. Y132F only) cannot support that and fails here even if N is large.
4. **Drug coverage ≥ the claimed set.** At minimum the azoles the #133 contract will report
   (fluconazole + voriconazole for *Candida*; the relevant triazoles for CYP51A). A panel missing
   a drug we intend to call does not qualify us to calibrate that drug.
5. **Obtainable + licensed, with truth-set discipline.** Public accessions (SRA/ENA) or an
   open/licensed supplement we may redistribute as a pinned fixture. Genotype-from-genome must be
   derivable via the existing caller path. Any manual transcription follows the FKS1 rule —
   **verified manual step, hashed fixture, no LLM scrape** (the scrape failure mode in
   `PREREGISTRATION_fks1_concordance.md` is why).

## Pre-committed verdict (both directions, fixed now)

    PASS  — proceed to #133
        A collection meeting ALL five criteria is identified and pinned (source, N, class
        balance, variant table, drug coverage, license) into a frozen panel manifest before #133
        starts. The "calibrated probability" contract is backed and may be built.

    FAIL — do NOT build the calibrated-probability product
        No obtainable collection meets the bar (for either target). Pre-committed fallbacks, and
        which one is triggered is decided by WHY it failed, not by re-reading the bar:
          (a) A near-miss exists — real paired data but too small / too hotspot-skewed to
              calibrate (fails #2 or #3). Fallback: NARROW the day-one claim to
              concordance-only — drop "calibrated probability" from the #133 contract, keep the
              hard S/R call + the honest "uncharacterized" verdict, and record the panel size
              that WOULD unlock calibration.
          (b) No paired genotype+phenotype collection is obtainable at all (fails #1 or #5).
              Fallback: the interpretation-engine MVP is not yet buildable on public data;
              stop, and record what data would need to exist. Do not fabricate a panel.

**Anti-rationalization clause (load-bearing).** The five criteria do not move after a candidate
is seen. A collection that misses any one is a FAIL that triggers a fallback above — never a
prompt to redefine "calibration-grade" down to fit it. If the numbers in this doc are ever
judged wrong, that is a new, separately-frozen pre-registration with its own SHA and a written
reason, not an edit to this one.

## What is NOT in scope here (declared, not silently dropped)

- **Building anything.** This gate only decides *whether* the calibrated-probability MVP is
  buildable and *on which target*. The wrap (#134), calibration (#136) and validation (#137)
  are downstream and separately pre-registered where they set their own bars.
- **The clinical claim.** Everything here is RUO. A qualifying panel unlocks *calibration*, not
  clinical actionability; the go/no-go for that lives in #139.
- **MIC regression.** Even when quantitative MICs are present, the day-one calibration target is
  the binary R/S phenotype; modelling continuous MIC is a later, separate question.
