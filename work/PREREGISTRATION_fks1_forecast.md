# Pre-registration — FKS1 (echinocandin) prevalence forecast

**Written 2026-10-17, as the forecasting companion to
[PREREGISTRATION_fks1_backtest.md](PREREGISTRATION_fks1_backtest.md).** Everything below is
fixed in advance. Deviations are recorded as deviations, never edited back into this
document. This extends the surveillance arm from a *hindcast* (what prevalence is / was) to a
*forecast* (where it is heading), under the same proof-before-claim discipline
(docs/MOONSHOT.md item 3).

## Why this prereg exists

The weather page reports the current and per-batch-historical echinocandin-resistance picture
but never projects forward — the one thing a "weather" service owes that a report does not.
This prereg fixes exactly what a forecast will claim, and the interpretation of **both**
outcomes (skill / no-skill), before reading any held-out error number.

## The data, and its honest limit (stated before we start)

The only committed series is the release-710 FKS1 call table
(`data/earlywarning/fks1_calls/PDG000000067.710.tsv`), one row per **resolved** isolate,
carrying `isolate_key, fks1_call, fks1_resistance_source, target_creation_date`. The only
dated field is `target_creation_date` (NCBI **deposit** date), spanning ~**2026-04 → 2026-10**
(7 month bins). It is read **offline**; no NCBI pull, no snapshot, no GCP. The exact table
sha256 is recorded in [RESULTS_fks1_forecast.md](RESULTS_fks1_forecast.md) so the run rebuilds.

Declared limits, true of every outcome:
1. **Deposit-time, not incidence.** The series is prevalence among *deposited* genomes per
   month. Collection lags deposit (~97-day median) and deposit volume tracks sequencing
   programmes, not transmission, so this is **not** epidemiological incidence. The page states
   this; the forecast is of the deposited series only.
2. **7 month bins.** A 7-point monthly series cannot support a strong parametric forecast. The
   model is deliberately the simplest honest one (below); the headline deliverable is the
   **certified attempt**, not a forced win.
3. Prevalence is **panel-gated** (known-resistance tokens via `fks1_caller.is_panel_token`;
   VUS never counted) over the **resolved** denominator (`emergence.resolved`) — the same
   gates the trend/geo sections use.

## Method, frozen (no new tunable)

The computation is the committed `openafr.forecast` module, not anything introduced to fit a
result:
- **Series**: `forecast.monthly_prevalence_series` — resolved-denominator, panel-gated,
  binned by `target_creation_date[:7]`; undated/unparseable rows excluded, never guessed.
- **Model**: `forecast.forecast_next` — least-squares trend on the Anscombe-corrected logit
  ((n+0.5)/(N+1)) of monthly prevalence vs a 0-based month index, inverted to a probability,
  with a Student-t 95% **prediction** interval (logit scale, inverted). No threshold, no knob
  that could be reverse-fit.
- **Naive baselines** (the bar to beat): `naive_last` (persistence — last observed month's
  prevalence) and `naive_pooled` (pooled prevalence over all training months).

## What we report (interpretation of each outcome pre-committed)

### P1 — the monthly series (REAL, descriptive)
The per-deposit-month panel-gated prevalence series over the committed table. No bar — a
measured descriptive trajectory, pre-committing only *how* it is computed.

### P2 — walk-forward skill test (the forecast claim)
`forecast.backtest(series, min_train=3, horizon=1)`: for each split leaving ≥3 training months
and one held-out month one step ahead, fit on the prefix, forecast the held-out month, and
record the absolute error of the model AND of both naive baselines against that month's actual
prevalence. No month is ever used to both fit and score itself.

**Pre-committed skill bar (`forecast.has_skill`): the model's held-out MAE must be strictly
less than BOTH `naive_last_mae` AND `naive_pooled_mae`.** A tie or beating only one baseline is
**NOT** skill.

Interpretation, fixed now:
- **SKILL** (model beats both baselines) → the page surfaces the one-step-ahead projection
  with its interval and a plain-language direction line, labelled a forecast.
- **NO SKILL** (model ties or loses to a naive baseline) → the page shows the monthly
  **trajectory only** and states plainly that no forecast skill over a flat baseline has been
  demonstrated at this series depth. This is a fully acceptable, likely outcome: if prevalence
  is roughly stationary, the pooled mean *should* be hard to beat, and saying so is the honest
  result. We do **not** tune the model to manufacture skill.

### P3 — interval coverage (REAL, descriptive)
The fraction of held-out months whose actual prevalence fell inside the model's 95% prediction
interval. Descriptive calibration context; no bar (7 bins make a coverage *claim* under­powered).

## Pre-committed overall verdict

What this delivers: the **first forward projection** of the surveillance arm, certified by a
walk-forward that can only promote a projection to the page when it out-predicts naive
persistence — and otherwise honestly reports a stationary series with no demonstrated forecast
skill. A deeper, multi-release series (which would let the trend model actually have something
to fit) remains **blocked on repeated billable fills**, the single named remaining blocker.
