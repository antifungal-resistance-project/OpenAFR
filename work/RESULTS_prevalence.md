# Azole event frequency in *C. auris* — the measured deliverable (issue #27)

**The number the whole resistance early-warning track was gated on.** How common is a
known ERG11 (CYP51) azole-resistance substitution in the sequenced *C. auris* population?
Until now this was `NOT MEASURED YET` — `erg11_call` was empty for every isolate because
NCBI runs no AMR pipeline on *C. auris* and the ERG11 re-caller had not been *run* on a
real snapshot. This is the first representative measurement.

**Code:** `openafr/recaller.py` (deterministic caller) · `scripts/recall_erg11.py`
(`fill --random --seed`) · `openafr/backtest.py` (`prevalence`, Wilson CI) ·
`scripts/backtest_earlywarning.py prevalence`.
**Provenance:** `data/earlywarning/runlog/fill.jsonl` run_id `e4a0def6586a`
(started 2026-08-19T23:38Z, finished 2026-08-20T16:29Z; code `d0a131d`).

## Snapshot + sample

- **Snapshot:** NCBI release **PDG000000067.675**, 29,294 isolates, 96.6% SRA-callable.
  Empty-snapshot content sha256 `d25ba2a3…`; filled sha256 `320482612a…`.
- **Sample:** `fill --random --seed 0 --limit 200` — a **seeded, reproducible random
  sample** of the 28,298 callable-pending isolates (NOT the first-N-in-snapshot-order,
  which correlates with resistance and gave a biased 100% at n=20). 199 called, 1 partial
  (excluded), 0 failed. ~17 h wall time on a GCP e2-standard-2 (download-bound).
- **Sanity gate (same tools/params) passed 4/4 first:** B11220 WT clean, B11221→F126L,
  B11243/B11245→Y132F, all 524/524 residues covered.

Reproduce:

```bash
python scripts/snapshot_ncbi_auris.py pull                      # -> PDG000000067.675
python scripts/recall_erg11.py fill --random --seed 0 --limit 200 <snapshot>
python scripts/backtest_earlywarning.py prevalence <snapshot>
```

## Headline

> **Azole event frequency = 160 / 199 = 80.4%** (resolved ERG11 panel carrying a known
> resistance mutation). **Wilson 95% CI [74.3%, 85.3%] (±5.5%).**

Excluded from the denominator (not counted azole-negative): 1 partial, 0 failed/refused,
29,094 pending. **39 of the 199 resolved isolates carry only *non-panel* substitutions**
→ azole-negative (novelty is `emergence.py`'s job, not the frequency's).

### By panel mutation (of 199 resolved)

| mutation | count | share |
|---|---|---|
| V125A | 70 | 35.2% |
| F126L | 70 | 35.2% |
| K143R | 47 | 23.6% |
| Y132F | 43 | 21.6% |

(An isolate can carry more than one; V125A+F126L co-occur.)

## What it means

1. **The biased n=20 (100%) was a clustering artifact, confirmed.** The representative
   base rate is **~80%, not ~100%** — there is a real, measurable ~20% of the population
   without a known panel mutation.
2. **Still a high, near-saturated baseline.** Azole-panel resistance is already present in
   roughly four of five sequenced *C. auris* isolates.
3. **The FKS1 / echinocandin (v2) redirection is now confirmed by data, not hunch.** An
   early-warning system for a mutation already in ~80% [74–85%] of isolates has little to
   warn about — azole *emergence* is a weak signal against this baseline. The dynamic,
   clinically scary axis is echinocandin/FKS1 resistance, which is still rare. The v2 data
   gate (TODOS.md) has **fired**.
</content>
</invoke>
