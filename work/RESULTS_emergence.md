# Emergence detection (issue #22)

**What:** the core method NCBI does not provide. Given the snapshot store (#21),
it diffs a snapshot's recent window against its baseline and flags ERG11 (CYP51)
azole-resistance substitutions that are **first-appearing** or **rising** — the
signal a metadata feed cannot give you (spike #20).

**Code:** `openafr/emergence.py` (library) · `scripts/detect_emergence.py` (CLI:
`run` on a real snapshot, `demo` on synthetic calls) · `tests/test_emergence.py`
(10 offline tests). Reproduce: `python scripts/detect_emergence.py demo`.

## The two signals, defined honestly

Split a snapshot by NCBI visibility date (`target_creation_date`) into a
**baseline** (created ≤ `as_of`) and a **recent window** (`as_of` < created ≤
`as_of` + `window_days`). Then, over *called* isolates only, per substitution:

- **FIRST-APPEARANCE** — present among recent isolates, never in the baseline.
- **RISING** — recent share of called isolates exceeds baseline share by
  ≥ `min_delta`.

Both require ≥ `min_count` recent carriers, so one sequenced isolate cannot raise
an alarm. Defaults: `window_days=180`, `min_count=3`, `min_delta=0.05`,
`backlog_frac=0.5`. Every verdict echoes its thresholds so it is self-describing.

## Three honesty constraints that shaped the method

1. **The denominator is *called* isolates, never all isolates.** The store leaves
   `erg11_call` empty until the re-caller fills it and records *why* (#21). This
   detector inherits that: an empty call is **not** counted as susceptible.
   Proportions are taken over called isolates and the calling coverage is
   reported next to every verdict, so "6% carry Y132F" is never read off 2%-called
   data.

2. **Undated isolates are set aside, not guessed.** An isolate with no usable
   `target_creation_date` cannot be placed in the baseline-vs-recent split, so it
   is returned in a separate `undated` bucket and counted, never bucketed by
   assumption (same discipline as `baseline_as_of`).

3. **The backlog artifact has a guard.** A lab depositing a years-old backlog
   makes old isolates appear *now* — a spike in the recent window that is really
   archival. The tell is **deposit latency**, not sample age: a carrier is
   archival if it became visible ≥ `max_lag_years` (default 2) after it was
   *collected*. If ≥ `backlog_frac` of a mutation's carriers are archival, the
   signal is annotated `possible_backlog` **with the reason** and sorted below
   genuine signals — flagged, not silently dropped, not blindly trusted.

## What it does on real data today: an honest "cannot detect yet"

Run against the live snapshot (`PDG000000067.671`, 29,085 isolates), baseline cut
2025-01-01:

```
coverage: baseline 0/15432 called, recent 0/6019 called, 0 undated
NOTE: no called isolates in the recent window -- cannot detect mutation
emergence yet. Calls are pending the ERG11 re-caller; an empty call is honestly
uncalled, never assumed susceptible.
no emergence signals over the configured thresholds.
```

There are 6,019 isolates in the recent window and 15,432 in the baseline, but
**zero** carry a resistance call, because NCBI runs no AMR pipeline on C. auris
(#20) and the re-caller has not yet populated `erg11_call`. The correct output is
*not* "all clear" — it is "cannot detect yet." The method is complete and
unblocked; its **input** (resistance calls) is what remains blocked, and it says
so rather than emitting a false negative.

## Proving the method fires: the synthetic demo

Because real calls are pending, `demo` builds a synthetic called snapshot that
exercises every branch (and doubles as a fixture for when calls arrive):

```
4 signal(s) (3 genuine, 1 possible-backlog):
  F126L   [FIRST-APPEARANCE, RISING]        baseline 0/6 (0.0%) -> recent 4/13 (30.8%)  delta +30.8%
  Y132F   [RISING]                          baseline 5/6 (83.3%) -> recent 13/13 (100.0%) delta +16.7%
  K143R   [RISING]                          baseline 1/6 (16.7%) -> recent 4/13 (30.8%)  delta +14.1%
  TR34    [FIRST-APPEARANCE, RISING, possible-backlog, non-canonical-token]
                                            baseline 0/6 (0.0%) -> recent 5/13 (38.5%)  delta +38.5%
          backlog: 5/5 carriers deposited >=2y after collection (archival 100%); likely a backlog
```

F126L (a real C. auris azole substitution) is caught as a genuine
first-appearance; Y132F and K143R as rising; and the TR34 spike — every carrier
collected in 2019 but deposited in 2025 — is correctly demoted as a
`possible_backlog`. TR34 also shows the parser keeping a non-canonical token
(it is a promoter allele, not a point substitution) rather than dropping a call
the re-caller emitted.

## Issue #22 checklist

- [x] Flag **first-appearance** mutations (unseen in the baseline window).
- [x] Flag **rising** mutations (recent share up vs baseline by ≥ `min_delta`).
- [x] Thresholds defined and documented (`min_count`, `min_delta`, `window_days`,
      `backlog_frac`, `max_lag_years`), echoed in every verdict.
- [x] Backlog guard (old-collection / fresh-visibility spike flagged with reason).

## What this unblocks / hands off

- **#23 (mutation → pocket mapping)** consumes each flagged substitution
  (`signal["mutation"]`, canonical tokens filtered by `is_substitution`) and maps
  it onto the modeled C. auris ERG11 structure. That is where the moat starts.
- **The real signal stays blocked on the ERG11 re-caller.** Until `erg11_call` is
  populated from the linked SRA reads, `run` on live snapshots will keep
  returning the honest empty state. The method is ready; the calls are the gap.
