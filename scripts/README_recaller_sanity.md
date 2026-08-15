# ERG11 re-caller — known-genotype sanity run

The cheap de-risking slice before the ERG11 re-caller is trusted on a real 29k-isolate
snapshot (TODOS.md, "Build the ERG11 re-caller" → OPEN step 1). It runs the full
reads → consensus → call orchestration on a handful of *published* known-genotype
*C. auris* strains and checks the emitted panel call against the literature.

A PASS here vouches for the real `scripts/recall_erg11.py fill` path, because the sanity
harness reuses that script's exact orchestration (`reads_to_consensus`) — same tools,
same alignment preset, same consensus parameters.

## What you need

The orchestration downloads SRA reads and aligns them, so it needs the bioconda tools
(now pinned in `environment.yml`) plus network and a few GB of scratch per strain:

```bash
conda env create -f environment.yml     # or: mamba env create -f environment.yml
conda activate openafr
# provides fasterq-dump (sra-tools), minimap2, samtools >=1.13
```

Platform note: bioconda builds these for **linux-64 / osx-64**, not osx-arm64. Run the
sanity harness (and any real `fill`) on a Linux/x86 host — the deterministic `call` path
and the tested core stay stdlib-only and run anywhere. The harness gates on the three
binaries and exits with an install message if they are absent (it never fakes a call).

## Run it

```bash
# cheapest smoke test first — the reference strain against its own reference:
python scripts/recaller_sanity.py --only B8441

# then the full set:
python scripts/recaller_sanity.py

# inspect an alignment workdir if a row surprises you:
python scripts/recaller_sanity.py --only B11221 --keep-tmp
```

Exit code is nonzero iff a row **FAIL**s. `INCONCLUSIVE` (a panel residue was uncovered,
so the caller honestly refused to guess) does not fail the run — widen coverage and retry.

## The truth set

`../data/earlywarning/recaller_sanity/known_genotypes.tsv` — the canonical clade
reference strains from Lockhart et al. 2017 (*Clin Infect Dis* 64:134), whose SRA runs
are in BioProject **PRJNA328792**:

| strain | clade | run | expected panel | why it's a good check |
|---|---|---|---|---|
| B8441  | I (Ia) | SRR3883469 | *wild-type* | the reference CDS (XM_085597798.1) is B8441-derived — its own reads **must** call wild-type; spurious panel hits mean the alignment/consensus is wrong, not the science |
| B11220 | II  | SRR3883452 | *wild-type* | clade II is azole-susceptible; confirms non-panel clade SNPs are emitted but not mis-tagged as resistance |
| B11221 | III | SRR7701905 | F126L | clade III (V125A may co-occur as the VF125AL haplotype) |
| B11243 | IV  | SRR3883464 | Y132F | clade IV |
| B11245 | IV  | SRR3883466 | Y132F | second, independent clade-IV confirmation |

The `expected_panel` value is a literature-derived **hypothesis the run confirms**, not an
oracle: a mismatch means look at *both* the caller and the expectation. In particular, if
B11221 emits `F126T` instead of `F126L`, that is the documented nomenclature subtlety
(the re-caller's panel tags `F126L`), not a bug — the harness prints the actual token so a
human can judge.

## After it passes

Widen to a real snapshot:

```bash
python scripts/recall_erg11.py fill data/earlywarning/snapshots/<release>.tsv --limit 20
```

then re-read the azole event frequency via `scripts/backtest_earlywarning.py replay` —
the measured number that gates both the azole MVP result and the FKS1/v2 decision
(TODOS.md).
