# Results — FKS1 (echinocandin) real-data backtest

**Run 2026-10-16**, per the frozen [PREREGISTRATION_fks1_backtest.md](PREREGISTRATION_fks1_backtest.md).
First exercise of the OpenAFR surveillance detector on **real** *C. auris* resistance calls —
the azole arm's ["structurally impossible today"](PREREGISTRATION_backtest.md) blocker is
retired for FKS1 because the certified re-caller's billable output is committed.

## Provenance (rebuildable, offline)

- Call table: `data/earlywarning/fks1_calls/PDG000000067.710.tsv` (release **PDG000000067.710**,
  as-of 2026-04-14, `window_days=180`, 3375 resolved isolates).
- `INDEX.tsv` sha256: `c0456491ae3b30b8a04867ef7d8b9042c7a6747282734c4e58ca5204563a8103`.
- Read directly with `backtest.load_fks1_call_table` — **no NCBI pull, no snapshot blob, no GCP.**
- Detector/prevalence = the pre-existing committed `openafr.backtest` functions with the frozen
  track thresholds (`window_days=180, min_count=3, min_delta=0.05, backlog_frac=0.5`); no new tunable.

Reproduce:
```
python scripts/backtest_earlywarning.py fks1-prevalence data/earlywarning/fks1_calls/PDG000000067.710.tsv
python scripts/backtest_earlywarning.py fks1-replay     data/earlywarning/fks1_calls/PDG000000067.710.tsv
```

## P1 — echinocandin prevalence (REAL)

| quantity | value |
|---|---|
| resolved panel (denominator) | 3359 / 3375 (excluded: 15 partial, 1 failed, 0 pending) |
| panel-positive (known-resistance) | **126** |
| **echinocandin event frequency** | **3.8%** (126/3359) |
| 95% CI (Wilson) | **[3.2%, 4.4%]** (±0.6%) |
| resolved, non-panel variant only | 69 (variants of uncertain significance — **not** counted resistant) |

The resolved denominator excludes partial/failed isolates (never assumed echinocandin-negative),
the "uncalled is not wild type" discipline shared with the azole arm.

### ⚠ Reconciliation with the live page — a real discrepancy this backtest surfaces

The weather page's headline counts **197** "echinocandin-resistant FKS1 carriers"
(`earlywarning.count_fks1_resistant`, ≈5.9% of 3359), because it counts **any** non-empty
`fks1_call` token. The pre-registered, **panel-gated** event frequency is **126 (3.8%)**. The
gap is **71 isolates carrying only NON-panel variants** — hot-spot-region substitutions of
**uncertain significance**, not the caller's known-resistance panel:

`D642H, D642N, D642V, F635I, F635L, F635S, L638F, M690V, P643H, R641G, R1353G, R1354C, R1354H, W1351C`

These are *not* the known-resistance alleles (D642**Y**, F635**C/Y**, M690**I**, R1354**S**, …).
Labeling them "resistant" overstates confirmed echinocandin resistance by ~56% (197 vs 126). The
honest headline is **two numbers**: ~3.8% carry a *known* echinocandin-resistance substitution;
a further ~2% (69 isolates) carry a hot-spot variant of *uncertain significance* warranting review.

**Fixed in this change.** The weather page now headlines **known-resistance panel carriers**
(`earlywarning.count_fks1_panel_positive`) and surfaces the uncertain-significance carriers
separately, so it no longer calls a variant of uncertain significance "resistant." The page
count is **128/3375 ≈ 3.8%** (panel-token carriers over all committed resolved rows); this
pre-registered figure is the stricter **126/3359 = 3.8%** (panel-positive over *fully* panel-
resolved isolates — the 2-isolate / 16-isolate differences are partial-window rows the stricter
resolved gate excludes). The two agree to displayed precision; both are honest, panel-gated, and
far from the old 197. (cf. the earlier coverage-denominator fix.)

## P2 — allele census (REAL, resolved carriers)

| panel allele | carriers | share of resolved | first visible |
|---|---|---|---|
| D642Y | 48 | 1.4% | 2026-04 |
| M690I | 38 | 1.1% | 2026-04 |
| S639F | 20 | 0.6% | 2026-04 |
| S639Y | 9 | 0.3% | 2026-04 |
| F635Y | 5 | 0.1% | 2026-05 |
| S639P | 3 | 0.1% | 2026-05 |
| W691L | 2 | 0.1% | 2026-05 |
| R1354S | 1 | 0.0% | 2026-05 |

Every known allele is already present at **window-open (2026-04/05)**. There is no
first-appearance inside the committed window — a direct consequence of its ~6-month depth, and
the concrete reason no lead-time claim is made (see P3).

## P3 — real-data walk-forward (the azole H2 null, retired)

7 monthly as-of steps across the real creation-date range 2026-04-15 … 2026-10-09:

| as_of | callable | signals |
|---|---|---|
| 2026-04-01 | True | 12 (incl. S639F) |
| 2026-05-01 | True | 4 |
| 2026-06-01 | True | 2 |
| 2026-07-01 | True | 2 |
| 2026-08-01 | True | 1 |
| 2026-09-01 | True | 0 |
| 2026-10-01 | True | 0 |

**Every step is `callable=True`** — the detector runs on real calls at every timepoint. The
azole arm's "cannot detect yet" null is **retired for FKS1**: this is the first real-data
exercise of the surveillance detector.

**Interpretation (pre-committed).** The signal count *decays* 12 → 0 as the as-of date advances.
This is the pre-registered **baseline-starvation artifact**, not a real falling trend: with
`window_days=180` over a 180-day-deep table, the earliest steps compare the recent window against
a near-empty baseline, so everything present reads as "newly appearing" (12 signals, incl. the
S639F flag); as the baseline accrues, genuine newness drops to 0. The honest read is the *late*
steps: **no allele is newly rising once a baseline exists.** No emergence/rising or lead-time
claim is made.

## Overall verdict (pre-committed)

- **Delivered:** the first real-data run of the FKS1 surveillance detector (P3), plus committed,
  CI'd prevalence (P1) and allele census (P2) reproducible from committed data alone.
- **Surfaced:** a page-headline overcount (197 vs the panel-gated 126) — known resistance vs.
  variants of uncertain significance.
- **Still NOT-YET-VALIDATED:** historical emergence / lead-time, now blocked solely on a
  **wider-window fill (billable GCP)** — the single named remaining blocker, replacing the azole
  arm's "no re-caller."

## Limitations (per prereg)
1. **180-day-deep window** → no emergence/rising or lead-time claim; only the real detector
   replay + current-feed census.
2. **US-skewed feed** (~88% USA): "US + thin international," not global.
3. Prevalence is over the **resolved** denominator; partial/failed excluded and reported, never
   assumed echinocandin-negative.
