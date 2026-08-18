# ERG11 re-caller — run tracking (sanity + fill)

**What:** the durable, git-diffable record of *every time the re-caller runs on real
data* — the tool+network paths that print to a terminal (sanity) or overwrite a snapshot
in place (fill), and whose outcomes therefore vanished before this existed. The re-caller
science lives in `openafr/recaller.py` (deterministic, tested); this doc is about
**provenance**: proving what ran, on what, with which tools and code.

**Code:** `openafr/runlog.py` (append-only JSONL logger) · wired into
`scripts/recaller_sanity.py`, `scripts/recall_erg11.py {fill,recall}` ·
`tests/test_runlog.py` + `tests/test_recall_fill_logging.py` (18 offline tests).
The logs live in `data/earlywarning/runlog/{sanity,fill,recall}.jsonl` (committed;
see that dir's README for the field schema).

## What "tracking every run" means here

Every invocation of the three un-reproducible paths appends **one JSON line** capturing:

- **when** — `started_at` / `finished_at` (UTC), and a `run_id`;
- **what code** — `code_commit` short-SHA + `code_dirty` (was the tree modified);
- **what tools** — `fasterq-dump` / `minimap2` / `samtools` `--version` strings, so a
  call is tied to the exact aligner/consensus build that produced it;
- **what reference** — the ERG11 CDS filename + its sha256, so a call is tied to what it
  was called *against*;
- **the outcome** — a run-level `summary` (sanity: pass/fail/inconclusive; fill:
  called/partial/failed/attempted) and a per-unit `items` list (per strain / per isolate);
- for **fill**, the snapshot's `snapshot_sha256_before` / `_after` and whether it `wrote` —
  so a re-`fill` that overwrites calls leaves a trail, and a crashed mid-batch run
  (`status: "error"`, equal digests, `wrote: false`) is honestly shown to have changed
  nothing.

The record is written **on exit even if the run crashes** — a failed run is logged with
whatever accumulated plus the error, never lost. Logging is best-effort and can never fail
or alter the science run (a write failure warns to stderr and moves on).

## Reading the history

```python
from openafr import runlog
for r in runlog.read_runs("sanity"):     # oldest first
    print(r["started_at"], r["status"], r.get("summary"), r["code_commit"])
```

## Sanity-run ledger (curated)

The machine-readable truth is `data/earlywarning/runlog/sanity.jsonl`; this table is the
human summary a reviewer skims. Append a row after each run.

| date (UTC) | code | tools (samtools/minimap2) | result | notes |
|---|---|---|---|---|
| _pending first run_ | — | — | — | Needs a Linux/x86 host with the bioconda tools (osx-arm64 has no bioconda build — see `scripts/README_recaller_sanity.md`). Run `python scripts/recaller_sanity.py`, then paste the `sanity.jsonl` summary here. |

`INCONCLUSIVE` (a panel residue was uncovered, so the caller honestly refused to guess)
does not fail the run — widen coverage and retry; it is recorded, not hidden.

## Fill history

`data/earlywarning/runlog/fill.jsonl` is the per-fill ledger. Because a `fill` overwrites
`erg11_call` in place, **the before/after snapshot digests in the log are the only record
of what a given fill changed** — pair them with `snapshots/INDEX.tsv` (the content hash of
the committed snapshot) to reconstruct the state at any point. The measured azole-event
frequency that gates the MVP (TODOS.md) should always be read *alongside* the `fill` run
that produced the calls it was measured on.

## Why this shape (not a second RESULTS write-up per run)

The other pipeline stages have one-time `RESULTS_*.md` analyses. The re-caller is
different: it is meant to run *repeatedly* as new isolates land, so the record has to be
**append-only and machine-readable**, mirroring the existing `snapshots/INDEX.tsv`
convention — small enough to commit, structured enough to diff and to query. This doc is
the human companion; the `.jsonl` files are the memory.
