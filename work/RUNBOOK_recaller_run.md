# Runbook — the ERG11 re-caller's first real run

The operational companion to TODOS.md ("Build the ERG11 re-caller" → OPEN steps 1–2).
Everything here has been de-risked offline; this is the order to run it in on the day,
and how to read each step. **Run on a Linux/x86 host** — the reads→consensus tools are
bioconda linux-64/osx-64 only (see `environment.yml`), and the step pulls multi-GB SRA
downloads. The deterministic core and all tests stay stdlib-only and run anywhere.

## 0. Prerequisites (once, on the run host)

```bash
git checkout main && git pull            # must include PRs #42, #43, #44
conda env create -f environment.yml      # or: mamba env create -f environment.yml
conda activate openafr
python -m pytest -q                       # sanity: the tested core is green here
samtools --version | head -1              # must be >= 1.13 (samtools consensus)
```

- **Merge PRs #42 #43 #44 first.** They provide, respectively: the samtools version gate +
  offline `plan`; the `prevalence` deliverable command; the all-Illumina sanity fixture.
  Running against a `main` missing any of them reintroduces a failure this runbook assumes
  is closed.
- **Disk:** budget ~5 GB scratch for the sanity set, and for `fill` roughly the per-isolate
  read size × `--limit` (Illumina *C. auris* WGS runs are ~0.5–1.5 GB each). `/tmp` must
  hold one isolate's reads at a time; each isolate's workdir is deleted after it is called.
- The orchestration **gates on the tools**: if `fasterq-dump` / `minimap2` / `samtools`
  are absent or samtools is < 1.13, `recall`/`fill`/sanity exit up front with an install
  message — before any download. That gate is deliberate; do not bypass it.

## 1. Pull a fresh snapshot (metadata only — no big downloads)

```bash
python scripts/snapshot_ncbi_auris.py pull        # latest release; writes a snapshot TSV + updates INDEX
```

- Verified live 2026-08-18: release **PDG000000067.673**, 29,153 isolates, **96.6%
  SRA-callable**. The snapshot TSV is multi-MB and gitignored; the small INDEX manifest is
  committed as the durable record.
- Note the written path (e.g. `data/earlywarning/snapshots/PDG000000067.673.tsv`); every
  step below takes it as `<snapshot>`.

## 2. Preflight — know the download scope *before* committing (offline, no tools)

```bash
python scripts/recall_erg11.py plan <snapshot> --list-runs | head
```

Reports how many rows a `fill` would re-call (= SRA downloads), how many are skipped for
lacking a run accession, and the distinct accessions. This runs anywhere and touches no
network — use it to size the run and pick a starting `--limit`. **Do not run `fill`
unbounded**: ~28k callable isolates is ~28k downloads.

## 3. Sanity run — prove the orchestration on known genotypes first

```bash
python scripts/recaller_sanity.py --only B11220     # wild-type control first (must emit NO panel hit)
python scripts/recaller_sanity.py                    # then the full Illumina set
```

- The fixture (`data/earlywarning/recaller_sanity/known_genotypes.tsv`) is four Illumina
  paired-WGS strains from Lockhart 2017 (BioProject PRJNA328792): **B11220** (wild-type
  control), **B11221** (F126L), **B11243** and **B11245** (Y132F) — a wild-type + three
  resistance calls across clades II/III/IV. See `scripts/README_recaller_sanity.md`.
- **Read it:** a spurious panel hit on **B11220** means the *orchestration*
  (alignment/consensus), not the science, is wrong — stop and inspect (`--keep-tmp` keeps
  the workdir). `PASS` on the resistance strains vouches for the real `fill` path (same
  tools/preset/params). `INCONCLUSIVE` (a panel residue uncovered) is honest, not a fail —
  it means depth fell below the caller's floor (`MIN_DEPTH = 10`, a constant in
  `recall_erg11.py`, not a CLI flag); a deeper run fixes it, or lower that constant only
  with care. `FAIL` (expected token missing, or a hit on a WT) is the real stop condition.
- If B11221 emits `F126T` instead of `F126L`, that is the documented nomenclature subtlety,
  not a bug — the harness prints the actual token to judge.
- **Do not proceed to `fill` until the sanity run is green.**

## 4. Fill — re-call real isolates, start small then widen

```bash
python scripts/recall_erg11.py fill <snapshot> --limit 20         # start here
# widen once the first batch looks right:
python scripts/recall_erg11.py fill <snapshot> --limit 500
```

- `fill` writes `erg11_call` + `resistance_source` **in place** for pending rows, and logs
  every run (append-only) to `data/earlywarning/runlog/fill.jsonl` with before/after
  snapshot digests — so a re-run can never silently overwrite without a trace. `--dry-run`
  re-calls but does not write (note: it still downloads).
- Honesty to expect in the output: a covered isolate with no substitution →
  `sra-recaller:wild-type`; a panel residue uncovered → `sra-recaller:partial(...)` with an
  **empty** call (never a silent wild-type); a fetch/align failure → `sra-recaller:failed`
  (left for retry, not dropped); a frameshift → `sra-recaller:refused`. A long-read isolate
  (rare on NCBI) degrades safely to `partial`, never a wrong call.

## 5. The deliverable — measure the azole event frequency

```bash
python scripts/backtest_earlywarning.py prevalence <snapshot>
```

- **This is the number the whole run exists to produce** (TODOS step 2) and the gate on the
  FKS1/v2 decision. It reports panel-positive / **resolved** panel, with partial / failed /
  pending **excluded** (never counted azole-negative — the same "uncalled is not wild type"
  rule the caller enforces).
- Before `fill` populates calls it prints **"NOT MEASURED YET"** (honest: cannot measure,
  not a zero signal). After `fill` it prints the measured frequency + per-mutation breakdown.
- Also run the walk-forward null and the arrival budget for the write-up:
  `backtest_earlywarning.py replay <snapshot>` and `... lag <snapshot>`.

## 6. Record the result

- Commit the appended `runlog/*.jsonl` line(s) — the durable provenance of what was run.
- Add the measured event frequency to `work/RESULTS_backtest.md` (and the sanity ledger in
  `work/RESULTS_recaller_sanity.md`).
- **Then the FKS1/v2 fork (TODOS):** a *measured* near-zero azole frequency justifies the
  echinocandin redirection; a non-trivial frequency keeps the azole MVP as the live track.
  Either way it is now a measured number, not a measured-zero-because-unmeasured null.

## Abort conditions (stop and investigate, don't push through)

- samtools < 1.13, or any gate/tool error → fix the environment; do not bypass the gate.
- Sanity **B11220** shows a panel hit, or any strain **FAIL**s → orchestration bug; a wider
  `fill` would multiply it across the snapshot.
- `fill` marks a large fraction `failed` (not `partial`/`wild-type`/`called`) → likely a
  systematic fetch/align problem, not per-isolate biology. Re-check on one isolate with
  `recall <SRR> --keep-tmp`.
