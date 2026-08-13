# Snapshot store + baseline (issue #21)

**What:** the persistence layer for the C. auris azole early-warning track. Each
NCBI Pathogen Detection pull becomes a normalized, versioned, byte-stable
snapshot, giving emergence detection (#22) a well-defined "what we saw before"
plus the two derived views it needs — a baseline as of a date, and a diff between
two pulls.

**Code:** `openafr/earlywarning.py` (library) · `scripts/snapshot_ncbi_auris.py`
(CLI) · `tests/test_earlywarning_snapshot.py` (9 offline tests) ·
`data/earlywarning/` (store + README). Reproduce: `python
scripts/snapshot_ncbi_auris.py pull`.

## What the spike forced on the design

The spike (#20, `RESULTS_earlywarning_spike.md`) proved NCBI is a metadata +
genome-pointer feed, **not** an AMR feed: it calls zero azole determinants for
C. auris. So the store carries isolate metadata and the SRA genome pointer now,
and an **explicitly pending** resistance field for the ERG11 re-caller (#23) to
fill. `erg11_call` stays empty; `resistance_source` records why —
`pending:sra-recaller` when reads are linked, `pending:no-reads` when not. An
uncalled isolate is honestly uncalled, never silently "wild type". This is the
one design choice that keeps #22 from mistaking "not yet analyzed" for
"susceptible".

## Design decisions (load-bearing)

1. **Stable isolate identity = `target_acc` minus its `.N` version.** NCBI bumps
   the version when it re-processes the same isolate; the base accession persists
   across releases, so that is the diff key.
2. **Snapshots are byte-stable per release.** Rows are sorted by isolate key and
   the pull timestamp is kept *out* of the data. Re-pulling a release yields an
   identical file → idempotent, and `git diff` between releases shows only real
   change. Verified: two re-pulls of 671 produced the identical SHA-256.
3. **Provenance lives in a small committed INDEX, blobs don't.** Each snapshot is
   ~3 MB and regenerable, so the `.tsv` blobs are gitignored. `INDEX.tsv`
   (release tag, pull time, isolate count, content hash, source URL) is committed
   — that is the durable, diffable memory. A hash in the INDEX lets anyone
   rebuild a byte-identical snapshot and verify it.
4. **Stdlib only** (urllib + csv + json + hashlib), matching the probe. No pandas
   into an rdkit-only dependency surface.

## Live run (2026-08-13)

Latest release **PDG000000067.671**, **29,085 isolates**, 0 unkeyable rows
dropped. Numbers reproduce the spike:

| Metric | Value |
|--------|-------|
| Isolates with SRA reads (ERG11-callable by #23) | 28,093 / 29,085 = **96.6%** |
| Isolates with a creation date (baseline-placeable) | 29,085 / 29,085 = **100%** |
| Top country | USA 25,680 (**88%**) — US-centric, as the spike flagged |
| `erg11_call` populated | **0** (correct — nothing to call yet; #23's job) |

**Baseline as of 2025-01-01:** 15,432 isolates across 43 countries (0 undated,
i.e. every isolate is placeable in time). This is the concrete answer to "what
had NCBI seen by date D".

**Diff 670 → 671 (the emergence-detection input):** 14 new isolates, 40 removed,
0 version-revised. The 14 new are all USA, created 2026-08-12, each with SRA
reads — exactly the substrate #22 will scan and #23 will genotype. The 40 removed
show NCBI does withdraw/rename isolates between releases; the diff surfaces that
churn honestly rather than hiding it.

## What this unblocks / hands off

- **#22 (emergence detection):** consumes `diff_snapshots` (new/rising isolates)
  and `baseline_as_of`. The method — first-appearance vs. rising, artifact guards
  — is now the only missing piece; the memory it needs exists.
- **#23 (ERG11 re-caller):** the critical path per the spike. It reads `run_acc`
  from a snapshot and writes back `erg11_call` / `resistance_source`. The schema
  column is already reserved for it.

## Honest limits

- The store records *isolate arrival*, not *resistance*. Until #23 lands, every
  `erg11_call` is empty by construction — this layer alone cannot warn on
  anything; it is memory, not signal.
- Lead-time (collection_date → target_creation_date lag) is not yet computed here;
  the spike flagged it for #27. This store has the raw fields to measure it.
- Coverage is 88% USA. Any "global" framing downstream (#24/#25) must be
  qualified, as the spike already recorded.
