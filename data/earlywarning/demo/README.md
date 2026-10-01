# Re-caller demo samples

A tiny, self-contained sample dataset for the **one-command reproducible run** of both
re-callers (issue C1, #162). It exercises the deterministic `call` path — stdlib-only, no
conda env, no network, no bioinformatics tools, any platform (including osx-arm64).

## Run it

From the repo root:

```bash
python scripts/demo_recallers.py          # run all four samples, print a table
python scripts/demo_recallers.py --check   # also verify each against EXPECTED.txt
```

Each sample is a single-record consensus **CDS FASTA** — exactly what the re-caller's
`call` subcommand consumes:

```bash
python scripts/recall_erg11.py call --consensus-fasta data/earlywarning/demo/erg11_Y132F.fasta
python scripts/recall_fks1.py  call --consensus-fasta data/earlywarning/demo/fks1_S639F.fasta
```

## The samples

Two per gene — a wild-type susceptible control and one canonical resistance mutation.
Each is derived from the project's **pinned reference CDS** by a single documented codon
edit (the same way `tests/test_recaller.py` synthesises inputs), so every byte traces back
to a committed, hash-pinned source:

| file | gene | derivation from reference | expected call |
|---|---|---|---|
| `erg11_wildtype.fasta` | ERG11 | verbatim `erg11_reference/erg11_cds.fasta` | *(none — wild-type)* |
| `erg11_Y132F.fasta` | ERG11 | codon 132 `TAC → TTC` (Tyr→Phe) | `Y132F` |
| `fks1_wildtype.fasta` | FKS1 | verbatim `fks1_reference/fks1_cds.fasta` | *(none — wild-type)* |
| `fks1_S639F.fasta` | FKS1 | codon 639 `TCx → TTx` (Ser→Phe) | `S639F` |

`Y132F` is the dominant global azole-resistance mutation; `S639F` is a canonical HS1
echinocandin-resistance mutation. The wild-type samples are susceptible controls that must
emit **no** panel hit — a spurious call on them would mean the caller, not the sample,
is wrong.

## `EXPECTED.txt`

The committed oracle: the `call` / `source` / `panel` lines each sample must produce.
`demo_recallers.py --check` compares the live callers against it and exits non-zero on any
mismatch, and `tests/test_demo_recallers.py` runs that check under CI — so the samples and
the callers can never silently drift apart.

## Scope

This demo covers only the reproducible-anywhere `call` path. The full reads → consensus
orchestration (`recall` / `fill` / sanity) needs the bioconda toolchain and multi-GB SRA
downloads; it is documented in `work/RUNBOOK_fks1_run.md` and the ERG11 runbook, not here.
