# Runbook — the FKS1 (echinocandin) re-caller's first real run

The v2 companion to [`RUNBOOK_recaller_run.md`](RUNBOOK_recaller_run.md) (ERG11/azole) and to
TODOS.md ("Widen early-warning detection to echinocandin/FKS1 resistance (v2)" → OPEN step 1:
*a real FKS1 `fill` → measured prevalence*). The FKS1 detection pipeline is built end-to-end
(PRs #111–#114) and de-risked offline; what it has never had is a run on real reads, so
`recall_fks1.py prevalence` still prints **"NOT MEASURED YET"**. This is the order to run it in
on the day, and how to read each step.

**Run on a Linux/x86 host** — the reads→consensus tools (`fasterq-dump` / `minimap2` /
`samtools ≥ 1.13`) are bioconda linux-64/osx-64 only (no osx-arm64 build), and the step pulls
multi-GB SRA downloads. The deterministic core and all tests stay stdlib-only and run anywhere.
This is the workload that needs the Google Cloud (or equivalent Linux/x86) VM.

## Why this run matters (read before renting the box)

Unlike the ERG11 run, this one has real scientific payoff. The measured azole prevalence came
back ~80% saturated ([`RESULTS_prevalence.md`](RESULTS_prevalence.md)), which is exactly why v2
exists: azole *emergence* is a weak early-warning signal at that baseline. FKS1/echinocandin
resistance is the genuinely dynamic axis, and it is expected to be **low-prevalence** — which is
what would make it a *useful* trigger. So the number this run produces (a low, non-saturated
echinocandin frequency, or an alarming rising one) is the deliverable the whole v2 track was
built to measure.

Honest scope reminder: FKS1 is **detection-only**. Echinocandins coordinate no metal, so the
CYP51/heme-iron structural moat does not transfer to a glucan synthase — there is no structural
verdict for this half, by design (TODOS v2 step 2, deferred). This run produces a *prevalence*,
not a structural so-what.

## 0. Prerequisites (once, on the run host)

```bash
git checkout main && git pull            # must include PRs #111 #112 #113 #114
conda env create -f environment.yml      # or: mamba env create -f environment.yml
conda activate openafr
python -m pytest -q                       # sanity: the tested core is green here
samtools --version | head -1              # must be >= 1.13 (samtools consensus)
```

- **Merge PRs #111–#114 first.** They provide, respectively: the snapshot `fks1_call` schema +
  migration + batch `plan`/`fill`/`prevalence`; the alert compose/render surface; the scheduled
  detection-only delivery; and the turnkey wiring. A `main` missing any of them reintroduces a
  gap this runbook assumes is closed.
- **Disk:** budget ~5 GB scratch, and for `fill` roughly per-isolate read size × `--limit`
  (Illumina *C. auris* WGS runs are ~0.5–1.5 GB each). `/tmp` must hold one isolate's reads at a
  time; each isolate's workdir is deleted after it is called. FKS1 aligns the same reads to a
  different reference and calls only two ~9-codon windows, so the compute per isolate is trivial
  — **the cost is entirely the SRA downloads**, identical to ERG11.
- The orchestration **gates on the tools**: if `fasterq-dump` / `minimap2` / `samtools` are
  absent or samtools is < 1.13, `recall`/`fill` exit up front with an install message — before
  any download. That gate is deliberate; do not bypass it.

## 1. Pull a fresh snapshot (metadata only — no big downloads)

```bash
python scripts/snapshot_ncbi_auris.py pull        # latest release; writes a snapshot TSV + updates INDEX
```

- A freshly pulled snapshot **already carries the `fks1_call` / `fks1_resistance_source`
  columns** natively (the v2 schema, PR #113), each seeded to an honest pending state from
  `run_acc`. Note the written path (e.g. `data/earlywarning/snapshots/PDG000000067.NNN.tsv`);
  every step below takes it as `<snapshot>`.
- **Only if you reuse a pre-v2 snapshot** (one written before #113), migrate it in place first —
  idempotent, offline, no tools:
  ```bash
  python scripts/recall_fks1.py migrate <snapshot>   # adds fks1 columns; a current snapshot is left unchanged
  ```

## 2. Preflight — know the download scope *before* committing (offline, no tools)

```bash
python scripts/recall_fks1.py plan <snapshot> --list-runs | head
# preview the exact representative sample you intend to fill:
python scripts/recall_fks1.py plan <snapshot> --random --seed 0 --limit 400
```

Reports how many rows a `fill` would re-call (= SRA downloads), how many are skipped for lacking
a run accession, and the distinct accessions. Runs anywhere, touches no network — use it to size
the run and pick a starting `--limit`. **Do not run `fill` unbounded**: ~28k callable isolates is
~28k downloads. The `run_acc` set is shared with ERG11, so this scope should match the ERG11
preflight.

## 3. Sanity run — prove the orchestration on known genotypes first

This is now turnkey, exactly like ERG11's. A committed FKS1 sanity fixture
(`data/earlywarning/recaller_sanity/fks1_known_genotypes.tsv`) and harness
(`scripts/recaller_sanity_fks1.py`) run four published-genotype strains through the *same*
reads→window-consensus orchestration the real `fill` uses, so a PASS vouches for the `fill`
path (same tools, preset, per-window `samtools consensus -r`).

```bash
python scripts/recaller_sanity_fks1.py --only B11220     # wild-type control first (must emit NO S639 hit)
python scripts/recaller_sanity_fks1.py                    # then the full set (WT + one strain per S639 token)
```

- The fixture is a **wild-type negative** control plus **one positive per panel token**:
  - **B11220** (`SRR3883452`, clade II) — echinocandin-*susceptible*; its FKS1 windows must
    emit **no S639 panel hit**. A spurious `S639*` there means the orchestration
    (alignment/consensus), not the biology, is wrong. (Same strain/run as the ERG11 WT control.)
  - **B11211** (`SRR28270963`) → **S639F**, **B12136** (`SRR10461137`) → **S639P**,
    **B12663** (`SRR7909222`) → **S639Y** — the three positive controls, one per panel token,
    from the published FKS1 benchmark dataset (PMC12323592). All four runs are ENA-verified
    Illumina PAIRED WGS (see the fixture header), matching the real `fill` path's `-ax sr` preset.
- **Read it:** `PASS` on the resistant strains vouches for the real `fill` path. `FAIL` (expected
  token missing, or a hit on the WT) is the real stop condition. `INCONCLUSIVE` is honest, not a
  fail — the panel window was uncalled (coverage gap), refused (an in-window indel), or missing;
  the caller correctly refused to guess. A deeper run or `--keep-tmp` (keeps the workdir to
  inspect) helps. `MIN_DEPTH = 10` (a constant in `recall_fks1.py`, not a CLI flag) is the
  per-position depth floor below which a residue is uncalled.
- Every sanity run is logged append-only to `data/earlywarning/runlog/` with per-strain outcomes,
  same as ERG11 — a durable PASS/FAIL record, not a scrolled-past line.

Panel to expect (from `data/earlywarning/fks1_reference/PROVENANCE.md`): **HS1 = S639F/P/Y**
(the dominant, panel-tagged call); **HS2 (R1354)** is reported but *not yet* panel-tagged, so a
non-S639 window difference is emitted verbatim, not flagged as resistance. A window uncovered
below `MIN_DEPTH` is refused for *that window* with a reason and left uncalled — never a silent
wild-type, and never fatal to the other window.

**Do not proceed to a wide `fill` until the sanity run is green** (WT clean + at least one
positive token confirmed).

## 3b. Decide the sample size *before* you look (choose `--limit`)

Same decide-first discipline as the ERG11 runbook and the pre-registrations: pick `--limit` from
the precision you want to quote, *before* seeing the point estimate. But the FKS1 question is the
mirror image of ERG11's — you expect a **low** prevalence, so the relevant precision is the
**upper** confidence bound, not the worst-case (p=0.5) half-width.

- If FKS1 resistance is genuinely rare, `k` positives out of `n` resolved is small, and the
  Wilson CI is *tighter* than the p=0.5 table — but you need `n` large enough to see any
  positives at all. Rule of thumb: if the true rate is ~1%, `n=400` expects ~4 positives; `0/n`
  gives a ~`3/n` upper 95% bound (n=400 → "≤ ~0.75%").
- A pragmatic first target: **~400 resolved** distinguishes a near-zero echinocandin rate from a
  few-percent emerging one, and matches the ERG11 representative-sample scale. Read `n` as
  **resolved** isolates (partial/failed/pending don't count toward the denominator), so oversize
  `--limit` by ~10–15% for the unresolved fraction. `plan --random --seed S --limit N` previews
  the exact sample.
- Worst-case (p=0.5) half-widths, for reference — conservative here since real FKS1 p is low:
  n=200 → ±6.9%, n=384 → ±5.0%, n=500 → ±4.4%, n=1000 → ±3.1%. Regenerate with
  `python3 -c "from openafr import backtest as bt; print(bt.wilson_halfwidth_worst_case(400))"`.

## 4. Fill — re-call real isolates, start small then widen

```bash
python scripts/recall_fks1.py fill <snapshot> --limit 20                       # fast orchestration check
# widen once the first batch looks right -- sample RANDOMLY, not in snapshot order:
python scripts/recall_fks1.py fill <snapshot> --limit 500 --random --seed 0    # representative
```

- **Use `--random` for any prevalence you intend to quote.** Without it, `--limit N` re-calls
  the first N pending rows in *snapshot order*, which tracks NCBI submission order and clusters
  by study — a biased head, not the population rate (the ERG11 run's first non-random `--limit
  20` came back 18/18 from one outbreak). `--random` draws a seeded representative sample
  (`--seed`, default 0, reproducible). A small leading `--limit 20` *without* `--random` is still
  the right *first* step — it's a fast orchestration check, not a frequency estimate.
- `fill` writes `fks1_call` + `fks1_resistance_source` **in place** for pending rows, alongside
  (not overwriting) the azole `erg11_call` columns, and logs every run append-only to
  `data/earlywarning/runlog/` with before/after snapshot digests, so a re-run can never silently
  overwrite without a trace. `--dry-run` re-calls but does not write (it still downloads).
- Honesty to expect in the output, per window: a covered window with no substitution →
  `sra-recaller:wild-type`; a window residue uncovered → `partial(...)` with an **empty** call
  for that window (never a silent wild-type); a fetch/align failure → `failed` (left for retry,
  not dropped); a within-window indel/frameshift → `refused` for that window, without discarding
  the other. A long-read isolate degrades safely to `partial`, never a wrong call.

## 5. The deliverable — measure the echinocandin event frequency

```bash
python scripts/recall_fks1.py prevalence <snapshot>
```

- **This is the number the whole v2 run exists to produce** (TODOS v2 step 1). Note it is a
  *distinct* subcommand from ERG11 — `recall_fks1.py prevalence`, **not**
  `backtest_earlywarning.py prevalence`. It reports panel-positive (S639F/P/Y) / **resolved**
  panel, with partial / failed / pending **excluded** (never counted echinocandin-negative — the
  same "uncalled is not wild type" rule the caller enforces).
- Before `fill` populates calls it prints **"NOT MEASURED YET"** (honest: cannot measure, not a
  zero signal). After `fill` it prints the measured frequency, its **95% Wilson CI**, and the
  per-substitution breakdown. **Quote the interval, not the point estimate** — at a low true
  rate, the upper bound is the defensible statement ("echinocandin resistance ≤ X% in this
  representative NCBI sample").

## 6. Record the result

- Commit the appended `runlog/` line(s) — the durable provenance of what was run.
- Write the measured echinocandin frequency into a `work/RESULTS_fks1_prevalence.md` (mirror of
  `RESULTS_prevalence.md`), and record the sanity outcome (which controls, what they emitted).
- Update TODOS.md: v2 step 1 moves from OPEN to DONE-with-measured-number, exactly as the ERG11
  step 2 did on 2026-08-20.
- Interpretation to state honestly: a low, non-saturated FKS1 frequency **confirms** the v2
  premise (echinocandin is the dynamic axis with early-warning headroom, unlike ~80%-saturated
  azole); a surprisingly high one is itself a finding. Either way it is now a measured number,
  not a measured-zero-because-unmeasured null.

## 7. (Optional, same host) Grade the caller against the full external benchmark

The `fill` above measures *prevalence* trusting a caller whose accuracy is only a 4-strain
control (step 3). While the box is up, also run the pre-registered **concordance** validation
— the caller's measured sensitivity/specificity against the 100-isolate PMC12323592 benchmark
(TODOS v2 step 3; `work/PREREGISTRATION_fks1_concordance.md`, frozen sha `da8e3827…`). This is
the "widen n → measured CI" move that retired preprint Limitation #1, applied to the caller.

**Prerequisite done off-host:** transcribe PMC12323592 Table 1 (100 isolates: SRR, FKS1
genotype, R/S) into `data/earlywarning/recaller_sanity/fks1_benchmark_100.tsv` (schema in the
script header) and pin its hash. The paper is **not open-access** — do NOT auto-scrape it (an
LLM scrape produced a self-contradictory duplicate row); copy from the PDF and verify each SRR
live against ENA `read_run` (Illumina/PAIRED/WGS), then:

```bash
shasum -a 256 data/earlywarning/recaller_sanity/fks1_benchmark_100.tsv \
  >> work/PREREG_fks1_concordance.sha256          # freeze the fixture hash
python scripts/validate_fks1_concordance.py --limit 3   # cheap smoke test first
python scripts/validate_fks1_concordance.py             # full benchmark + pre-registered verdict
```

The grader re-checks the prereg + fixture hashes and refuses to certify if either moved. Read
the verdict: **PASS** → fold the measured accuracy into `RESULTS_fks1_prevalence.md`;
**UNDERPOWERED** (resolved < 80%) → more reads / lower `--min-depth`; **FAIL** → a named,
reproducible caller defect (which token, which direction) that blocks promoting the FKS1 track
past detection-of-record. The scoring engine (`openafr/concordance.py`) is unit-tested offline,
so a surprise here is the *reads/orchestration*, not the metric.

## Abort conditions (stop and investigate, don't push through)

- samtools < 1.13, or any gate/tool error → fix the environment; do not bypass the gate.
- The **B11220** negative control emits any `S639*` panel hit → orchestration bug; a wider
  `fill` would multiply a false-positive resistance call across the snapshot. Stop.
- `fill` marks a large fraction `failed` (not `partial`/`wild-type`/`called`) → likely a
  systematic fetch/align problem, not per-isolate biology. Re-check one isolate with
  `recall <SRR> --keep-tmp`.
- **Every** window comes back `partial` (uncovered) → the FKS1 reference/window coordinates
  aren't matching the assembled consensus; inspect a kept workdir before spending on more
  downloads. (The offline tests cover the caller logic, but the binaries end-to-end are exactly
  what this run is proving.)
