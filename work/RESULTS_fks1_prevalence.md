# Echinocandin (FKS1) event frequency in *C. auris* — the v2 measured deliverable

**The number the FKS1/echinocandin early-warning track (v2) was built to produce.** How common
is a known FKS1 hot-spot echinocandin-resistance substitution (S639F/P/Y) in the sequenced
*C. auris* population? Until now this printed `NOT MEASURED YET` — `fks1_call` was empty for every
isolate because NCBI runs no AMR pipeline on *C. auris* and the FKS1 re-caller had never been
*run* on real reads. This is the first representative measurement, and the mirror-image question
to the ERG11/azole one ([`RESULTS_prevalence.md`](RESULTS_prevalence.md), ~80% saturated).

**Code:** `openafr/fks1_caller.py` (deterministic windowed caller) · `scripts/recall_fks1.py`
(`fill --random --seed`, `prevalence`) · Wilson CI in `openafr/backtest.py`.
**Provenance:** `data/earlywarning/runlog/fill-fks1.jsonl` run_id `037ca9a103ea`
(code `9e53e8f`, GCP e2-standard-2, us-central1-a; VM deleted after the run).

## Snapshot + sample

- **Snapshot:** NCBI release **PDG000000067.686**, 29,699 isolates, 28,690 SRA-callable (96.6%).
  Empty-snapshot content sha256 `a74eeb26…`; filled sha256 `1f1876866ba4…`.
  FKS1 reference `fks1_cds.fasta` sha256 `da0cf3eb…`.
- **Sample:** `fill --random --seed 0 --limit 450` — a **seeded, reproducible random sample** of
  the callable-pending isolates (not first-N-in-snapshot-order, which clusters by study). Of 450
  attempted: **443 resolved**, 4 partial (window uncovered → excluded, never counted susceptible),
  3 failed (fetch/align error → left for retry, excluded). ~1 day wall time, download-bound.
- **Sanity gate (same tools/params) passed 4/4 first:** B11220 (clade II) wild-type clean
  (no S639 hit), B11211→S639F, B12136→S639P, B12663→S639Y. Tools: fasterq-dump 3.1.1,
  minimap2 2.28, samtools 1.21.

Reproduce:

```bash
python scripts/snapshot_ncbi_auris.py pull                          # -> PDG000000067.686
python scripts/recaller_sanity_fks1.py                              # 4/4 PASS gate first
python scripts/recall_fks1.py fill <snapshot> --limit 450 --random --seed 0
python scripts/recall_fks1.py prevalence <snapshot>
```

## Headline

> **Echinocandin event frequency = 10 / 443 = 2.3%** (resolved FKS1 panel carrying a validated
> S639F/P/Y resistance substitution). **Wilson 95% CI [1.2%, 4.1%].**

Excluded from the denominator (not counted echinocandin-negative): 4 partial, 3 failed, 29,249
pending.

### By panel substitution (of 443 resolved)

| substitution | count | note |
|---|---|---|
| S639F | 9 | validated HS1 panel token |
| S639Y | 1 | validated HS1 panel token |

### Non-panel hot-spot changes — reported verbatim, NOT counted as resistance (13 isolates)

The caller emits any hot-spot-window substitution, but only S639F/P/Y are the pre-validated,
panel-tagged resistance calls. These 13 are honestly *reported but excluded* from the frequency
(novelty adjudication is `emergence.py`'s job, not the prevalence's):

| change | count | region |
|---|---|---|
| D642Y | 6 | HS1 (near S639, not a validated token) |
| R1354H | 4 | HS2 (reported, not yet panel-tagged) |
| L640V | 2 | HS1 |
| F635C | 1 | HS1 |

## What it means

1. **A low, non-saturated echinocandin baseline — exactly the useful case.** At **2.3%
   [1.2%, 4.1%]**, FKS1 resistance is still rare in the sequenced population. Contrast the ERG11
   azole baseline of **~80%**: an early-warning system has real headroom here — there is a low
   background rate for a genuine *rise* to stand out against.
2. **The v2 premise is now confirmed by data, not hunch.** The redirection from azole to
   echinocandin (TODOS v2) was justified on the argument that echinocandin is the dynamic,
   clinically-scary axis while azole resistance is already ubiquitous. Both halves are now
   measured: azole ~80% (no headroom), echinocandin ~2% (headroom). The dynamic axis is FKS1.
3. **Honest scope reminder — detection only.** Echinocandins coordinate no metal, so the
   CYP51/heme-iron structural moat does not transfer to a glucan synthase. This run produces a
   *prevalence*, not a structural verdict (TODOS v2 step 2, deferred by design). The three failed
   isolates are a per-isolate fetch/align issue, not a systematic one (well under the abort
   threshold), and the WT sanity control was clean, so the calls are trustworthy.
