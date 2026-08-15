# Backtest / honest validation (issue #27)

**The credibility test for the whole resistance early-warning track (#20–26).**
Pre-registered in [PREREGISTRATION_backtest.md](PREREGISTRATION_backtest.md)
(sha256 `26e6d46f…`, locked **before** any number below was computed;
`work/PREREG_backtest.sha256`).

**Code:** `openafr/backtest.py` (library) · `scripts/backtest_earlywarning.py`
(CLI: `lag` / `replay` / `mechanism`) · `tests/test_backtest.py` (12 offline tests).
**Real snapshot:** NCBI release **PDG000000067.672**, 29,118 isolates, content sha256
`6e6d790fa9dd…` (`data/earlywarning/snapshots/INDEX.tsv`). Reproduce:

```bash
python scripts/snapshot_ncbi_auris.py pull
python scripts/backtest_earlywarning.py lag    data/earlywarning/snapshots/PDG000000067.672.tsv
python scripts/backtest_earlywarning.py replay data/earlywarning/snapshots/PDG000000067.672.tsv
python scripts/backtest_earlywarning.py mechanism
```

## Headline: NOT-YET-VALIDATED (blocked on the ERG11 re-caller) — honestly

The early-warning **claim** — "replayed over history, this would have flagged Y132F
before the published record" — **cannot be validated on real data today**, because the
pipeline's load-bearing input does not exist. NCBI runs no AMR pipeline on *C. auris*
(spike #20), so `erg11_call` is empty for all 29,118 isolates, and the **ERG11 re-caller
that would fill it was never built** (it is in no shipped milestone item; issue #23
became the pocket *mapping*). This is the pre-committed outcome, not a surprise, and per
the project ethos it is a valid result: *here is exactly what is missing, and here is the
arrival budget it would have to work within.*

Within that honest frame, the three pre-registered parts came back:

| Part | What it tests | Data | Pre-registered verdict |
|------|---------------|------|------------------------|
| **H1-arrival** | the lead-time *budget* (deposit lag) | **real** | ✅ **PASS-arrival** |
| **H2-null** | walk-forward on the real snapshot | **real** | ⬛ **HONEST NULL** (as committed) |
| **H3-mechanism** | lead time once calls exist | labelled/synthetic | ✅ **PASS-mechanism** |

## H1-arrival — the lead-time budget is real and adequate (PASS)

The one arrival-axis number spike #20 explicitly deferred to #27: the **deposit lag**
`collection_date → target_creation_date` (sample taken → visible in NCBI). It bounds any
achievable lead time — you cannot warn before the data is visible.

```
datable 23,376/29,118 (5,742 undatable, excluded not guessed)
lag days: median 97   p25 46   p75 527   p90 1107   [min 4, max 9211]
visible within 365d of collection: 16,154/23,376 (69%)
pre-registered verdict: PASS-arrival
  median lag 97d (<= 365d) and 69% visible within 365d (>= 60%)
```

**Median 97 days** from sample to NCBI visibility, **69%** visible within a year — both
clear the pre-committed bar (median ≤ 365 d, frac-within ≥ 60%). Against the *stated
assumption* that published surveillance of a novel resistance lineage lags collection by
1–3 years, a ~3-month median visibility lag leaves a genuine lead-time budget. The
arrival axis does **not** defeat the early-warning premise.

**The honest caveat is the tail, and we report it rather than smooth it:** p75 = 527 d
and p90 = 1,107 d — a long right tail of archival deposits (max 9,211 d ≈ 25 y). The
*typical* isolate is fast; a meaningful minority is a backlog, which is exactly why the
detector carries the backlog guard (#22). 5,742 isolates lack a usable date pair and are
excluded from the ratio, not imputed.

## H2-null — the honest real-data null (as committed)

Walk-forward over the snapshot's real creation-date range **2018-03-23 … 2026-08-13**
(102 monthly steps), tracking Y132F with the frozen detector defaults:

```
steps with >=1 called isolate in the recent window: 0/102
steps that flagged Y132F: 0
HONEST NULL: erg11_call is empty on real data -> every step returned
"cannot detect yet". The pipeline refusing to emit a false all-clear,
NOT a validated warning.
```

This is the pre-committed expectation, and it is the *correct* behaviour: across eight
years of replay the detector never once invents a warning from absent calls. A stable
"cannot detect yet" is the honesty contract holding, not a pass of the claim.

## H3-mechanism — the detector produces lead time once fed calls (PASS)

To locate the blocker — broken method, or starved of input? — the same walk-forward runs
over a **labelled (synthetic) replay** in which a Y132F-shaped mutation is first collected
mid-2023 and becomes NCBI-visible after a ~5-month deposit lag on three isolates:

```
walk-forward: 23 monthly steps, 4 flagged Y132F
first flag of Y132F: 2023-08-01
reference recognition: 2024-09-01
lead time: 397 days (13.2 months) BEFORE reference
pre-registered verdict: PASS-mechanism
```

The detector converts a first-appearance into a dated flag **397 days before** the
reference recognition date, using the frozen defaults and no per-run tuning. **This is a
mechanism proof, not real-data evidence** (its dates are illustrative, exactly as
`detect_emergence demo` is synthetic). It grants the track no real-data credibility — but
it does establish that H2's null is attributable to the **missing re-caller**, not to a
detector that cannot produce a warning.

## What this means for the track

1. **The single remaining blocker is named and isolated: the ERG11 re-caller.** Every
   other stage (#21–26) is built, tested, and — via H3 — shown to fire end-to-end once
   `erg11_call` is populated. The re-caller (reads → ERG11 substitution call, over the
   96.6% of isolates with linked SRA runs) is the whole gap between "scaffold" and
   "validated early warning."
2. **The arrival budget clears the bar (H1),** so building the re-caller is worth it: a
   ~3-month median visibility lag is compatible with beating the published record. This
   is the first evidence that the *early* in early-warning is defensible on real data.
3. **Azole-mutation event frequency on real data is currently zero-measurable** (H2), not
   because emergence is absent but because calls are absent. This is the input the FKS1/v2
   redirection TODO was data-gated on — see [../TODOS.md](../TODOS.md): the gate reads
   "cannot measure azole emergence until the re-caller exists," so v2 is **not** triggered
   by a measured near-zero; it is deferred behind the same re-caller.

## Issue #27 checklist

- [x] Replay historical NCBI snapshots and ask whether this would have flagged Y132F
      before surveillance — done as a real walk-forward (H2) + a mechanism replay (H3).
- [x] Report lead time (or lack of it) honestly — H1 gives the real budget (median 97 d);
      H2 gives the honest real-data null; H3 gives the mechanism lead time (397 d, synthetic).
- [x] Pre-register the success criterion before running — `PREREGISTRATION_backtest.md`,
      sha256-locked before any number was computed.
- [x] An honest negative is a valid outcome — the overall verdict is NOT-YET-VALIDATED
      (blocked on the re-caller), reported as such with the exact missing piece named.
```
