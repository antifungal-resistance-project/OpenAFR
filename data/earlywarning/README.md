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
  runlog/                        # committed: append-only provenance of every re-caller run
    {sanity,fill,recall}.jsonl   #   one JSON line per run (see runlog/README.md)
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
`erg11_call` starts blank and `resistance_source` records *why* —
`pending:sra-recaller` when SRA reads are linked (the ERG11 re-caller can fill it)
or `pending:no-reads` when they are not. An uncalled isolate is honestly uncalled,
never silently "wild type".

The **ERG11 re-caller** (`openafr/recaller.py` + `scripts/recall_erg11.py`) is what
turns `pending:sra-recaller` into an actual call — the piece that makes the whole
early-warning track fire on real data instead of an honest null. It manufactures the
call NCBI does not provide, from the isolate's own ERG11 sequence, against a pinned
B8441 reference (`erg11_reference/`, accession XM_085597798.1, sha256-verified with a
numbering check on the V125/F126/Y132/K143 hotspots). Once it runs, `resistance_source`
becomes one of:

- `sra-recaller:called` — substitution(s) called (`erg11_call` populated),
- `sra-recaller:wild-type` — panel fully covered, no substitution,
- `sra-recaller:partial(panel-uncalled:…)` — a hotspot could not be covered (empty
  call, but *never* silently wild type),
- `sra-recaller:failed(…)` / `:refused(…)` — fetch/align failed, or the consensus was
  out of frame (indels can't be point substitutions). The isolate is marked, not dropped.

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

# --- ERG11 re-caller: fill the pending:sra-recaller calls ---

# Deterministic (no network, no external tools): call from a consensus CDS FASTA.
python scripts/recall_erg11.py call --consensus-fasta isolate_erg11.fasta

# Orchestrated (needs fasterq-dump + minimap2 + samtools>=1.13, and network):
python scripts/recall_erg11.py recall SRR1234567           # one isolate
python scripts/recall_erg11.py fill snapshots/PDG000000067.671.tsv --limit 5   # a snapshot, in place
```

The re-caller is split on purpose: the **call logic** (consensus CDS → tokens) is
stdlib-only, deterministic, and unit-tested against the pinned reference; the
**reads → consensus** orchestration shells out to bioinformatics tools and is gated on
their presence (exactly like the docking scripts gate on `vina`), so it is honest about
being the un-reproducible-here layer rather than faking a call when the tools are absent.

Libraries: `openafr/earlywarning.py`, `openafr/recaller.py`.
Tests: `tests/test_earlywarning_snapshot.py`, `tests/test_recaller.py`.
