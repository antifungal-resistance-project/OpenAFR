# Early-warning snapshot store (issue #21)

Persistence layer for the C. auris azole early-warning track. Turns each NCBI
Pathogen Detection pull into a normalized, versioned snapshot so emergence
detection (#22) has a well-defined "what we saw before".

## Layout

```
data/earlywarning/
  README.md                      # this file (committed)
  snapshots/
    INDEX.tsv                    # committed: provenance of every pull
    PDG000000067.<N>.tsv         # gitignored: the per-release snapshot blobs (~3 MB each)
```

**Only `INDEX.tsv` is committed.** The snapshots themselves are large and fully
regenerable, so they are gitignored (see repo `.gitignore`). `INDEX.tsv` is the
durable, git-diffable memory of the store: one row per pulled release with the
release tag, pull timestamp, isolate count, and the snapshot's content SHA-256.
Anyone can rebuild a byte-identical snapshot from the release tag and check it
against the recorded hash.

## Schema (one row per isolate)

`isolate_key`, `target_acc`, `biosample_acc`, `run_acc`, `asm_acc`,
`collection_date`, `target_creation_date`, `country`, `geo_loc_name`, `lat_lon`,
`erg11_call`, `resistance_source`.

`isolate_key` is `target_acc` with its `.N` version stripped — the identity that
persists when NCBI re-processes an isolate.

### `erg11_call` is deliberately empty

The spike (#20, `work/RESULTS_earlywarning_spike.md`) proved NCBI runs no AMR
pipeline on C. auris: there is no azole-resistance call to read off the feed. So
`erg11_call` stays blank and `resistance_source` records *why* —
`pending:sra-recaller` when SRA reads are linked (the ERG11 re-caller, #23, can
fill it) or `pending:no-reads` when they are not. An uncalled isolate is honestly
uncalled, never silently "wild type".

## Commands

```bash
# Pull the latest release, write a snapshot, update INDEX:
python scripts/snapshot_ncbi_auris.py pull

python scripts/snapshot_ncbi_auris.py pull --release 671   # pin a release
python scripts/snapshot_ncbi_auris.py pull --dry-run       # report, don't write

# Baseline: isolates NCBI had as of a date (by target_creation_date):
python scripts/snapshot_ncbi_auris.py baseline snapshots/PDG000000067.671.tsv --as-of 2025-01-01

# New isolates between two pulls (the emergence-detection input, #22):
python scripts/snapshot_ncbi_auris.py diff snapshots/PDG000000067.670.tsv snapshots/PDG000000067.671.tsv
```

Library: `openafr/earlywarning.py`. Tests: `tests/test_earlywarning_snapshot.py`.
