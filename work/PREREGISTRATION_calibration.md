# Pre-registration — calibrate the resistance probability, with a measured reliability (issue #136)

> **STATUS: FROZEN 2026-09-09, before any probability has been fit or any reliability number
> read.** This document fixes, in advance, the calibrator, the fit/test split, the reliability
> metrics, and the pass bar — *and* records that the run it governs is **BLOCKED on data**, not
> runnable today (see "The data this run needs does not exist yet"). The grader
> (`scripts/validate_calibration.py`) re-checks this document's SHA-256, the pinned caller
> reference hash, and the panel fixture's hash, and refuses to certify if any moved or if the
> fixture is absent. Freezing the metric before a panel exists is deliberate: it removes the
> freedom to choose the reliability bar after seeing which panel we could finally obtain. The
> SHA-256 is pinned in `work/PREREG_calibration.sha256`.

## Why this run is required — "calibrated probability" is a claim the contract defers to here

The diagnostics verdict contract (`docs/DIAGNOSTIC_VERDICT_CONTRACT.md`, issue #133) ships
**concordance-only**: a categorical S/R + "uncharacterized" verdict, with **no probability
field**, because the #132 spike (`work/RESULTS_diagnostics_panel_spike.md`) found no public
*Candida* collection that could calibrate one. The contract states the probability field
"re-enters only after that track produces a validated calibration — never before." This
pre-registration *is* that track's gate. Calibrating a probability means two things the
concordance check (`work/PREREGISTRATION_fks1_concordance.md`, sha `da8e3827…`) does not do:
**fit** a score→probability map, and **measure on a held-out split** whether a stamped 0.8
actually corresponds to ~80% observed resistance. Until that reliability is measured, any
emitted probability is an unbacked number, and this gate exists to keep it off the contract.

## The data this run needs does not exist yet (the honest block, committed now)

Per the frozen #132 gate (`work/PREREGISTRATION_diagnostics_panel.md`, sha `1211e973…`), **no
obtainable public *Candida* collection** clears the calibration bar (N≥150, ≥40/class, ≥5
variants with ≥5 carriers, paired genotype+MIC, single obtainable collection). The #132 verdict
was **FAIL → fallback (a), concordance-only**; the calibrated-probability product was deferred,
not killed, with its unlock route named: **option C — pool several public *C. albicans* ERG11
studies into one hashed fixture** (e.g. PRJNA592373 + Flowers/Morio/cohort studies). That route
has three pre-committed costs, none of which this document waives:

1. **A separate, separately-frozen prereg amendment** addressing cross-study MIC-method
   heterogeneity (pooling mixes CLSI/EUCAST breakpoints across labs — exactly what #132
   criterion #5's "single collection" rule guarded against). Pooling without that amendment is
   forbidden here.
2. **New *C. albicans* caller / structure work.** The re-callers and the structural port are
   *C. auris*-centric (5TZ1; [[openafr-auris-port]]); a *C. albicans* panel needs a caller path
   that resolves its ERG11, which does not exist yet.
3. **Verified manual transcription of the paired genotype+MIC tables**, hashed fixture, **no LLM
   scrape** — the same discipline (and the same paywall wall) as the FKS1 benchmark
   ([[fks1-caller-concordance]]).

So this run is **BLOCKED**. What is buildable and frozen *now* is the calibration **engine**
(`openafr/calibration.py`, unit-tested in `tests/test_calibration.py`) and this metric/bar,
mirroring how #129 froze the FKS1 concordance engine+prereg before its paywalled run. The
engine is ready; the run certifies nothing until the option-C fixture exists and is pinned.

## The score being calibrated (fixed now, per regime)

The two regimes the contract names **do not share a calibrator**, and one of them has **no
continuous score by design** — this is committed, not discovered later:

- **Known-variant regime.** Stratum = the specific panel token the caller emits (ERG11:
  V125A/F126L/Y132F/K143R; FKS1 HS1: S639F/S639P/S639Y). The calibrated probability is the
  **per-token empirical resistance rate** on the training split, with a Wilson 95% CI
  (`calibration.stratified_rates` → `apply_rates`). There is no sigmoid to fit: a named marker
  either is or is not present.
- **Novel-variant regime.** The structural module
  **`openafr/structural.py` deliberately emits a *direction*, never a magnitude** ("a direction,
  never a number") — so there is **no continuous score to feed a calibration curve**, and
  claiming one would fabricate a precision the module refuses to assert. The stratum here is the
  **ordinal structural evidence tier** (`direction` × `confidence`:
  weaker-fit/minimal-direct-effect/uncertain × measured/estimate/none), and its calibrated
  probability is the **per-tier empirical resistance rate** with a Wilson CI, again via
  `stratified_rates`. This is the honest calibration ceiling for the moat verdict.
- **Continuous calibration is pre-specified but not used day one.** If a genuinely continuous
  resistance score ever exists (a pooled-*C. albicans* **MIC regression**, #132 option C's bonus
  target), it is calibrated with **isotonic regression (PAVA)** as primary and **Platt scaling**
  as the parametric comparator (`calibration.fit_isotonic` / `fit_platt`), both already frozen
  and tested. The primary/day-one target remains the **binary R/S** phenotype, per #132.

## The split + procedure (fixed now)

- **Per drug, per regime, fit separately.** One calibrator is fit for each (drug ∈ reported
  azole/echinocandin set) × (regime ∈ known/novel) cell. They are never pooled across regimes
  (the contract's "will not share a curve"), and a cell with too few carriers to fit is reported
  as **UNDERPOWERED for that cell**, not force-fit.
- **Held-out test split, seeded.** A **70/30 stratified** fit/test split, stratified on the R/S
  label, with a **fixed seed recorded in this document: `seed = 20260909`**. The calibrator is
  fit on the 70% training split only; every reliability number below is computed on the disjoint
  30% test split. The engine enforces that fit and evaluate are separate calls
  (`tests/test_calibration.py::test_fit_train_then_evaluate_held_out_test_is_well_calibrated`).
- **Excluded, not guessed.** An UNRESOLVED isolate (no call) has no score and is excluded. A
  **test isolate whose stratum was never seen in training** has no honest probability and is
  reported as `n_unseen_stratum`, **never defaulted to the base rate** — the same "excluded, not
  counted negative" rule as `concordance.evaluate` / `backtest.panel_prevalence`.

## The reliability metrics (fixed now) — `openafr.calibration.evaluate_calibration`

Computed on the **held-out test split**, per cell:

1. **Brier score** — mean squared error of predicted probability vs the 0/1 outcome.
2. **Reliability curve** — 10 equal-width probability bins; per non-empty bin, mean predicted
   vs observed resistant frequency with a Wilson 95% CI. The diagonal is perfect calibration.
3. **Expected calibration error (ECE)** — count-weighted mean bin gap |predicted − observed|.
4. **Calibration-in-the-large** — mean predicted vs overall observed rate (catches a global
   bias a per-bin curve can mask).

## Pre-committed pass bar (fixed now, both directions)

Applied per cell, over the **held-out test split** of a panel that has cleared the frozen #132
bar (a PASS there is a precondition — this gate never runs on a sub-bar panel).

    PASS   if  ECE           <= 0.10
           AND Brier         <= 0.20
           AND calibration-in-the-large gap  |mean_pred - observed|  <= 0.10
           AND the test split has >= 15 resistant AND >= 15 susceptible isolates
               (the #132 criterion-2 floor, carried down to the TEST split)
    FAIL   if any gated metric is worse than its bar on a test split that IS powered
           (>=15/class). A named, reproducible miscalibration: which cell, which direction
           (over- or under-confident), with the test isolates that expose it. The probability
           field stays OFF the #133 contract until fixed.
    UNDERPOWERED  if the test split has < 15 of either class (coverage too thin to read a
           reliability curve) — reported as neither pass nor fail, with the achieved n and the
           panel size that would fix it, mirroring the FKS1 UNDERPOWERED outcome.

ECE ≤ 0.10 and Brier ≤ 0.20 are pre-committed judgment calls: a reliability curve within ~10
percentage points of the diagonal is the floor at which a probability is worth emitting over the
categorical verdict the contract already ships. They are not renegotiated after seeing a
near-miss (anti-rationalization clause, below).

## Interpretation, committed now (both directions)

- **PASS.** The stratified probability is reliable on held-out isolates to within the bar, per
  cell. The probability field **re-enters the #133 contract** for the passing cells only, each
  carrying its measured ECE/Brier and its Wilson-CI calibration map. Folded into the contract
  and the diagnostics write-up as the backing for "calibrated probability."
- **FAIL.** The probability is not reliable at the bar. A genuine, publishable negative that
  **keeps the contract concordance-only** for that cell: the honest S/R + "uncharacterized"
  verdict stands, the unbacked probability does not ship. The failing cell and direction are
  recorded.
- **The novel-regime result is a tier calibration, not a continuous curve.** By construction
  (structural.py emits no magnitude) it can only calibrate the ordinal evidence tiers; reporting
  it that way — not as a smooth score→probability curve — is the same honesty as the contract's
  "calibrated-LOW, never a measured call" on the structural best-guess (#135).

## Anti-rationalization clause (load-bearing)

The calibrator choice, the 70/30 split and seed, the four reliability metrics, and the ECE /
Brier / in-the-large / per-class-n bars **do not move after a panel is seen**. A cell that
misses any gated metric on a powered test split is a FAIL that keeps the contract
concordance-only for that cell — never a prompt to loosen the bar. Loosening the #132 panel bar
to manufacture a "calibratable" dataset is likewise forbidden here; that gate is frozen
separately. If any number in this document is later judged wrong, that is a new,
separately-frozen pre-registration with its own SHA and a written reason, not an edit to this
one.

## What is NOT in scope here (declared, not silently dropped)

- **Obtaining / pooling the panel.** This gate fixes *how* a probability is calibrated and
  measured; assembling the option-C *C. albicans* fixture (and its cross-study MIC-method
  amendment) is the separate, blocking prerequisite named above.
- **MIC regression as the day-one target.** Even when quantitative MICs are present, the
  calibration target is the binary R/S phenotype (#132 criterion 1); continuous-MIC modelling is
  the reserved isotonic/Platt path, not the day-one claim.
- **Clinical actionability.** Every calibrated probability is RUO; the go/no-go for a clinical
  claim is #139. A PASS here unlocks a *reliable research-use probability*, nothing more.
- **The caller's own error rate.** That is the concordance run (#129, FKS1) / the sanity
  harnesses; this gate assumes a resolved call and asks only whether the probability attached to
  it is reliable.
