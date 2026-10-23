# Results — FKS1 (echinocandin) prevalence forecast

**Run 2026-10-23.** Pre-registered in [PREREGISTRATION_fks1_forecast_v2.md](PREREGISTRATION_fks1_forecast_v2.md)
(hash re-checked by the grader before this file was written). Fully local — committed call
table read offline, no GCP.

- **Source table:** `data/earlywarning/fks1_calls/PDG000000067.710.tsv`
  sha256 `947ad9ab459c58be057d744be4bd54207a196762983fa96e512d3aa3c1d19ee7`
- **Series depth:** 7 deposit-month bins (2026-04 → 2026-10).

## P1 — monthly panel-gated prevalence (deposit date)

| Month | Panel-positive / resolved | Prevalence |
|---|---|---|
| 2026-04 | 19 / 369 | 5.15% |
| 2026-05 | 19 / 481 | 3.95% |
| 2026-06 | 18 / 651 | 2.76% |
| 2026-07 | 24 / 473 | 5.07% |
| 2026-08 | 9 / 641 | 1.40% |
| 2026-09 | 34 / 559 | 6.08% |
| 2026-10 | 5 / 201 | 2.49% |

Panel-gated (known-resistance tokens only) over the resolved denominator. Deposit-month axis,
**not** epidemiological incidence (see prereg limitation 1).

## P2 — walk-forward skill test

Pre-committed bar: model MAE must be strictly below **both** naive baselines.

| Metric | Value |
|---|---|
| Held-out folds | 4 (min train 3, horizon 1) |
| Model MAE | 2.93 pp |
| Naive-last MAE | 3.56 pp |
| Naive-pooled MAE | 2.01 pp |

| Fold (train → target) | Actual | Model | Naive-last | Naive-pooled | In 95% PI |
|---|---|---|---|---|---|
| 2026-06 → 2026-07 | 5.07% | 2.10% | 2.76% | 3.73% | ✓ |
| 2026-07 → 2026-08 | 1.40% | 4.05% | 5.07% | 4.05% | ✓ |
| 2026-08 → 2026-09 | 6.08% | 2.33% | 1.40% | 3.40% | ✓ |
| 2026-09 → 2026-10 | 2.49% | 4.84% | 6.08% | 3.88% | ✓ |

## P3 — interval coverage

Actual-in-95%-PI over held-out months: **100.00%** (descriptive; 7 bins
underpower a coverage *claim*).

## Verdict: NO SKILL — trajectory only

No forecast skill over a flat baseline at this series depth; the page shows the monthly trajectory only. Pooled prevalence ≈ 3.79%.

A deeper multi-release series (which would give the trend model something to fit) remains
blocked on repeated billable fills — the single named remaining blocker.
