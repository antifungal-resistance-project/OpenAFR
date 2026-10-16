# Pre-registration — FKS1 (echinocandin) real-data backtest

**Written 2026-10-15, as the FKS1 companion to [PREREGISTRATION_backtest.md](PREREGISTRATION_backtest.md).**
Everything below is fixed in advance. Deviations are recorded as deviations, never edited
back into this document. This extends the azole-arm discipline to the first real-data run of
the surveillance detector.

## Why this prereg exists — the azole blocker is retired for FKS1

The azole-arm backtest ([PREREGISTRATION_backtest.md](PREREGISTRATION_backtest.md)) declared
its emergence claim **structurally impossible today** for one reason, on the record: the ERG11
re-caller *"was never built,"* so `erg11_call` is empty on every real isolate and the
walk-forward (H2) is a **null by construction**, the lead-time proof (H3) **synthetic**.

That blocker does **not** hold for FKS1. The FKS1 re-caller was built and certified 100%
sens/spec on the 98-genome PMC12323592 benchmark
([RESULTS_fks1_concordance.md](RESULTS_fks1_concordance.md)), and its billable GCP output is
committed as a durable call table (`data/earlywarning/fks1_calls/PDG000000067.710.tsv`). So for
the **first time** the detector can run on **real** *C. auris* resistance calls. This prereg
fixes exactly what we will report from that run, and the interpretation of both outcomes,
before reading any per-month number.

## The data, and its honest limit (stated before we start)

The committed call table is a **180-day window** (`INDEX.tsv`: `window_days=180`), one row per
RESOLVED isolate, carrying `isolate_key, fks1_call, fks1_resistance_source, target_creation_date`.
Visibility dates span ~**2026-04 → 2026-10** (~6 months deep). It is read **offline**; no NCBI
pull, no snapshot blob, no GCP. The exact release tag + `INDEX.tsv` sha256 are recorded in
[RESULTS_fks1_backtest.md](RESULTS_fks1_backtest.md) so the run is rebuildable.

A ~6-month-deep table **cannot** support a baseline-vs-recent emergence/rising claim, nor a
lead-time claim against a years-old literature recognition date: with `window_days=180` the
early as-of steps have a thin or empty baseline by construction. That is a property of the
committed window, not of the method. Widening it is a billable GCP fill, deliberately out of
scope here.

## Method, frozen (no new tunable)

The computation is the **already-committed** code, not anything introduced to fit a result:
- Prevalence: `openafr.backtest.fks1_panel_prevalence` (resolved denominator via
  `emergence.resolved(..., "fks1_call")`; panel positivity via `fks1_caller.is_panel_token`;
  Wilson 95% CI via `backtest.wilson_interval`).
- Walk-forward: `backtest.walk_forward(..., call_field="fks1_call")` →
  `emergence.detect_emergence` with the **frozen track defaults** shared with the azole arm
  (`window_days=180, min_count=3, min_delta=0.05, backlog_frac=0.5`). This change adds **no**
  new threshold that could be reverse-fit to a desired signal.

## What we report (three parts, interpretation of each outcome pre-committed)

### P1 — echinocandin prevalence (REAL, descriptive)
Panel-positive / resolved over the committed table, point estimate + Wilson 95% CI. **No
pass/fail bar** — this is a measured descriptive number, pre-committing only *how* it is
computed. Pre-committed cross-check: this panel-gated number is the fraction carrying a
**known-resistance** FKS1 substitution; it may differ from the page's headline carrier count
(`count_fks1_resistant`, which counts *any* FKS1 call token). If they differ we report the gap
and its cause honestly rather than reconciling by choosing the larger number.

### P2 — allele census (REAL, descriptive)
Per-panel-token resolved-carrier counts and the earliest `target_creation_date` month each
token is visible in the committed window. Descriptive; method pre-committed, no bar.

### P3 — real walk-forward (the azole H2 null, retired)
Step `as_of` monthly across the table's real `target_creation_date` range, running the detector
at each step with the frozen defaults. Pre-committed interpretation of both outcomes:
- The detector now runs on **real** calls — the azole "cannot detect yet" null is retired for
  FKS1. Every step reporting `callable=True` is the correct, expected output (calls exist), and
  whatever signals it emits are reported as-is.
- **Pre-committed caveat (fixed now):** because the table is 180-day-deep and `window_days=180`,
  early as-of steps have a thin/empty baseline, so any early "emergence" is **baseline-starved
  by construction**, not a validated rising trend; later steps, once a baseline accrues, are the
  honest read. We therefore make **no** emergence/rising or lead-time claim from this run. If the
  per-step signal count instead *rose* toward the present it would still not be claimed as real
  emergence without a deeper window.

## Pre-committed overall verdict

What this delivers: the **first real-data exercise** of the FKS1 surveillance detector (P3,
retiring the azole null for this arm), plus the committed, CI'd prevalence (P1) and allele
census (P2) that back the live page's numbers **reproducibly from committed data alone**. A
full historical emergence/lead-time backtest remains **NOT-YET-VALIDATED — blocked on a
wider-window fill (billable)**, which is now the single named remaining blocker, replacing the
azole arm's "no re-caller."

## Declared limitations (fixed in advance, true of every outcome)
1. **180-day-deep window** → no emergence/rising or lead-time claim; only the real detector
   replay + current-feed census.
2. **US-skewed feed** (~88% USA, inherited from the NCBI *C. auris* feed, spike #20): "US +
   thin international," not global.
3. Prevalence is over the **resolved** denominator; partial/failed isolates are excluded and
   reported as counts, never assumed echinocandin-negative — the "uncalled is not wild type"
   discipline shared with the azole arm.
