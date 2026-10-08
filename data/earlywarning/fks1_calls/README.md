# FKS1 call tables — the durable home for the billable echinocandin calls (issue #171)

The resistance-weather page is a thin view over a live NCBI *C. auris* pull. That pull carries
**no echinocandin/FKS1 call** — NCBI runs no AMR pipeline on *C. auris* — so every row arrives
in an honest `pending:` state and the page shows the calm day-0 "watching since…" picture.

The calls are ours to manufacture, out-of-band, with the FKS1 re-caller
(`scripts/recall_fks1.py`), which runs **billably on GCP** (SRA download + align is heavy; see
the runbook). Those calls are **not regenerable without re-paying**. The multi-MB snapshot blob
they are computed in is gitignored (regenerable NCBI metadata), so the calls need their own
small, committed home — this directory — or the public page would snap back to day-0 on the
next weekly live pull.

## Layout

```
data/earlywarning/fks1_calls/
  README.md                  # this file
  INDEX.tsv                  # committed: one provenance row per fill
  PDG000000067.<N>.tsv       # committed: the compact call table for release N
```

## Call table schema (one row per RESOLVED isolate)

`isolate_key`, `fks1_call`, `fks1_resistance_source`, `run_acc`, `target_creation_date`.

"Resolved" = the re-caller has spoken (a non-`pending:` source), whatever its verdict
(`called` / `partial` / `refused` / `failed`). A `pending:` isolate the re-caller has not
reached is simply absent. `isolate_key` is the version-stripped `target_acc`, so the table
overlays cleanly onto a later live pull even after NCBI re-processes an isolate.

## How it is produced and consumed

- **Produced** by `scripts/recall_fks1.py fill <snapshot>` (non-dry-run): after writing the
  gitignored snapshot it (re)writes `PDG…<N>.tsv` from the snapshot's resolved rows and appends
  an `INDEX.tsv` provenance row (release tag, as-of, window, n_called, sha256, timestamp).
  Idempotent — a re-fill re-emits the union of calls-so-far, never duplicates.
- **Consumed** by `scripts/render_weather.py` and `scripts/deliver_alert.py`
  (`openafr.earlywarning.overlay_fks1_calls`): on a live pull / `--tsv`, the newest table here is
  overlaid by `isolate_key`, so the page/alert reflects the real echinocandin picture. An isolate
  with no committed call stays in its honest pending state — never silently wild-type.

This file and an empty-but-header INDEX are committed so the first real GCP fill only **appends**
— nothing else to wire. See `work/RUNBOOK_fks1_fill.md` for the end-to-end fill runbook.
