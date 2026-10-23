# Pre-registration v2 — FKS1 (echinocandin) prevalence forecast, precision-weighted trend

**Written 2026-10-23, before the weighted fit is run or any v2 held-out error is read.** This
amends [PREREGISTRATION_fks1_forecast.md](PREREGISTRATION_fks1_forecast.md) (v1, sha frozen in
`work/PREREG_fks1_forecast.sha256`). v1 is **not** edited — it stays as the historical record of
the unweighted model and its honest NO-SKILL outcome ([RESULTS_fks1_forecast.md] run 2026-10-17).
This document's SHA-256 is pinned in `work/PREREG_fks1_forecast_v2.sha256`; the grader
(`scripts/backtest_forecast.py`) re-checks it and refuses to certify if it moved after freezing.

Everything below is fixed in advance. Deviations are recorded as deviations, never edited back.

## Why a v2 is required — the v1 skill comparison was not apples-to-apples

v1 fit the logit-prevalence trend by **ordinary (unweighted) least squares** (`forecast._ols`):
every deposit-month carried equal leverage. But a month's prevalence is a **binomial
proportion**, and its logit has sampling variance that *shrinks with the month's resolved count
N* — delta method, `Var(logit p̂) ≈ 1 / (N p̃(1-p̃))`. In release 710 the monthly denominators
span **201 → 651 resolved isolates**; treating October (N=201) as exactly as informative as
June (N=651) over-weights the noisy months.

The baseline v1 must beat, `naive_pooled`, **aggregates raw counts** and is therefore already
*implicitly N-weighted*. So v1 pitted an **equal-weighted** trend against an **N-weighted** flat
line — an asymmetry that handicaps the trend on a point of standard statistical method, not
taste. v2 removes that asymmetry by fitting the trend with the same N-weighting.

This correction is justified **a priori** (heteroskedasticity of binomial proportions; the
standard treatment for regressing rates), **independent of the v1 loss margin**. It is adopted
because it is the correct method, not because it is expected to win — the no-reverse-fit
discipline (`docs/DIAGNOSTIC_GO_NO_GO.md`) requires exactly this ordering.

## Method, frozen (one principled change; no new tunable)

The computation is the committed `openafr.forecast` module.

- **Series** — UNCHANGED from v1: `forecast.monthly_prevalence_series` — resolved-denominator
  (`emergence.resolved`), panel-gated (`fks1_caller.is_panel_token`; VUS never counted), binned
  by `target_creation_date[:7]`; undated/unparseable rows excluded, never guessed.

- **Model (the v2 change)** — `forecast.forecast_next(..., weighted=True)`: a **weighted**
  least-squares trend (`forecast._wls`) on the Anscombe-corrected logit
  `logit((n+0.5)/(N+1))` of monthly prevalence vs a 0-based month index, inverted to a
  probability. The per-month weight is the inverse delta-method logit variance,

  > **w_i = N_i · p̃_i · (1 − p̃_i),  p̃_i = (n_i + 0.5)/(N_i + 1)**  (`forecast._binomial_weight`)

  i.e. the binomial information the month actually contributes. p̃ is the *same* corrected
  prevalence v1 already fits — no new quantity, no free knob, no threshold. The weights differ
  across months whenever `N_i·p̃_i(1−p̃_i)` does (i.e. in general); the weighted fit coincides
  with v1's OLS only in the degenerate case of equal weights across all months.

- **Prediction interval** — Student-t (two-sided 95%, n−2 d.o.f.) on the logit scale, inverted,
  built from the WLS fit:

  > **se = √( s² · ( 1/w̄ + 1/Σw + (x\*−x̄_w)² / Sxx_w ) )**

  with weighted residual variance **s² = Σ wᵢ rᵢ² / (n−2)**, weighted mean
  **x̄_w = Σ wᵢxᵢ / Σw**, and **Sxx_w = Σ wᵢ (xᵢ−x̄_w)²**. The future month's own-scatter term
  uses the **training mean weight w̄ = Σw / n** (via `1/w̄`), because the next month's N\* is
  unknown at forecast time; this is the honest stand-in and is fixed here. (With equal weights
  `1/w̄ = 1`, recovering v1's `1 +` prediction-interval term.)

- **Naive baselines (the bar to beat)** — UNCHANGED: `naive_last` (persistence) and
  `naive_pooled` (count-pooled prevalence over the training months).

## What we report (interpretation of each outcome pre-committed)

### P1 — the monthly series (REAL, descriptive)
Identical series to v1 (same table, same gates). No bar; a measured trajectory.

### P2 — walk-forward skill test (the forecast claim)
`forecast.backtest(series, min_train=3, horizon=1)` over the **weighted** model: for each split
leaving ≥3 training months and one held-out month one step ahead, fit on the prefix, forecast
the held-out month, and record the absolute error of the model AND of both naive baselines.
No month is ever used to both fit and score itself.

**Pre-committed skill bar (`forecast.has_skill`, UNCHANGED from v1): the model's held-out MAE
must be strictly less than BOTH `naive_last_mae` AND `naive_pooled_mae`.** A tie, or beating
only one baseline, is **NOT** skill.

Interpretation, fixed now:
- **SKILL** (weighted model beats both baselines) → the page surfaces the one-step-ahead
  projection with its interval and a plain-language direction line, labelled a forecast, noting
  the fit is precision-weighted.
- **NO SKILL** (ties or loses to a naive baseline) → the page shows the monthly **trajectory
  only** and states plainly that no forecast skill over a flat baseline has been demonstrated at
  this series depth — a fully acceptable, likely outcome for a roughly-stationary series. We do
  **not** tune the model to manufacture skill; the weighting is the only change, and it is fixed
  above regardless of which verdict it yields.

### P3 — interval coverage (REAL, descriptive)
Fraction of held-out months whose actual prevalence fell inside the model's 95% prediction
interval. Descriptive calibration context; no bar (7 bins underpower a coverage *claim*).

## Pre-committed overall verdict

v2 delivers the **fair-method** version of the surveillance arm's first forward projection:
a precision-weighted trend, graded by the same walk-forward that can only promote a projection
to the page when it out-predicts naive persistence AND the (N-weighted) pooled mean. Whether the
fair comparison now shows skill or still shows a stationary series with none, the reported number
is the honest one. A deeper, multi-release series remains **blocked on repeated billable fills**,
the single named remaining blocker — v2 changes the estimator, not the data depth.
