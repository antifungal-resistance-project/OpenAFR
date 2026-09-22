# Runbook — one GCP session, two tracks: FKS1 concordance + #137 echinocandin accuracy

The FKS1 truth fixture is frozen on `main` (PR #148 — `data/earlywarning/recaller_sanity/fks1_benchmark_100.tsv`,
98 scored isolates, hash pinned in `work/PREREG_fks1_concordance.sha256`). That single fixture is
the truth set for **both** open validations:

- **FKS1 concordance** (`scripts/validate_fks1_concordance.py`) — is the caller's *token* right vs
  the paper's genotype? (`work/PREREGISTRATION_fks1_concordance.md`, sha `da8e3827…`)
- **#137 verdict accuracy** (`scripts/validate_diagnostic_accuracy.py`, PR #146) — when the engine
  emits a *verdict*, how often does it agree with the isolate's *echinocandin R/S phenotype*, in
  CLSI-M23 error names (VME/ME)? (`work/PREREGISTRATION_diagnostic_accuracy.md`)

Both are gated on the **same** reads→consensus→`call_windows` pass over the 98 isolates (the
Linux/x86 + `fasterq-dump`/`minimap2`/`samtools≥1.13` + multi-GB SRA workload). Running that pass
**once** and harvesting its output feeds both graders — so this is one cloud session, not two.

## v2 RE-MEASURE — the current task (FKS1 HS3 panel extension, #154 landed)

> The v1 pass below already ran (2026-09-20): concordance PASSed, but the #137 echinocandin arm
> **FAILed** (VME 4/40 = 10.0%). PR #154 then extended the caller with an HS3 window + widened
> HS1/HS2 tags, exactly as frozen in `work/PREREGISTRATION_diagnostic_accuracy_v2.md` (sha in
> `work/PREREG_diagnostic_accuracy_v2.sha256`), and characterised the FAIL as pure marker coverage —
> so this is a **re-run of the same one pass on `main` at/after #154**, testing the pre-committed
> projection **VME 1/40 = 2.5% [0.4, 12.9]**. Follow sections 0–A–B–C below, with these deltas:

1. **Fixture path — do NOT overwrite the v1 fixture.** The v2 prereg requires the v1 fixture
   (`fks1_accuracy_98.tsv`, frozen FAIL, still consumed by `scripts/characterize_coverage_ceiling.py`)
   to be preserved. Re-harvest to a **new** path:
   ```bash
   python scripts/validate_fks1_concordance.py \
       --emit-accuracy-fixture data/earlywarning/recaller_sanity/fks1_accuracy_98_v2.tsv
   ```
2. **Pin into the v2 file.** Freeze the re-harvested fixture in `PREREG_diagnostic_accuracy_v2.sha256`
   (NOT the v1 file). The grader now reads pins from **both** files and re-checks the v2 prereg sha,
   so the exact freeze command the emit step prints already targets the v2 file for a `_v2` fixture:
   ```bash
   shasum -a 256 data/earlywarning/recaller_sanity/fks1_accuracy_98_v2.tsv \
     >> work/PREREG_diagnostic_accuracy_v2.sha256
   ```
3. **Grade the v2 fixture.** `python scripts/validate_diagnostic_accuracy.py --fixture
   data/earlywarning/recaller_sanity/fks1_accuracy_98_v2.tsv` — same frozen bar (inherited verbatim
   from v1). Report PASS/FAIL/UNDERPOWERED per the v2 prereg's committed interpretation (both
   directions), and which isolates, if any, stayed missed (coverage vs. mechanism).
4. **Concordance must still PASS.** The new HS3 window must not regress HS1/HS2 token concordance —
   confirm the concordance verdict from the same pass is still PASS before trusting the accuracy arm.

Everything else (host setup §0, the reads→caller pass §A, the abstention/UNRESOLVED handling §C) is
unchanged. The v1 commands in §A–C stay as the historical record of the FAIL.

## Offline pre-flight — already green (2026-09-17, this workspace)

Verified locally before renting the box, on `origin/main` merged with the #146 branch
(`jajjer/next-work-item-v20`, merges clean, no conflicts):

- `pytest test_accuracy test_concordance test_verdict test_fks1_caller test_recall_fks1_orchestration test_accuracy_harvest`
  → **78 passed** (incl. the new harvest logic).
- `validate_fks1_concordance.py --limit 1` → integrity gate PASSES both frozen hashes, then stops
  exactly at the missing bioinformatics tools (the only thing the cloud host adds).
- `validate_diagnostic_accuracy.py` (no `--fixture`) → cleanly reports **BLOCKED** (its paired
  fixture is produced by step A on the host).

So nothing offline is left: the harvest is implemented and unit-tested; the only missing
ingredient is compute (the reads→caller pass over the 98 isolates).

## 0. Host setup

Follow `work/RUNBOOK_fks1_run.md` §0 (same conda env, same tool versions, `samtools ≥ 1.13`).
`git checkout main && git pull` — must include PR #148 (fixture) and a merge of PR #146 (the
accuracy machinery: `openafr/accuracy.py`, `scripts/validate_diagnostic_accuracy.py`,
`work/PREREGISTRATION_diagnostic_accuracy.md` + its `.sha256`). Merge #146 to main first, or check
out a branch that has both.

## A. FKS1 concordance run (Track 1) — and harvest the #137 fixture in the SAME pass

```bash
python scripts/validate_fks1_concordance.py --limit 3   # cheap smoke test first
python scripts/validate_fks1_concordance.py \
    --emit-accuracy-fixture data/earlywarning/recaller_sanity/fks1_accuracy_98.tsv
```

Re-checks the prereg + fixture hashes and refuses to certify if either moved. Read the
pre-registered verdict:
- **PASS** → fold the measured sensitivity/specificity/identity (with Wilson CIs) into
  `work/RESULTS_fks1_prevalence.md`, and retire the caller's "4-strain sanity control only" caveat.
- **UNDERPOWERED** (resolved < 80%) → coverage too thin; neither pass nor fail. More reads / lower
  min-depth, per that runbook's step 7.

`--emit-accuracy-fixture` writes the #137 paired fixture **from this same reads→caller pass** — no
second download. Run it **without `--limit`** so all 98 isolates are harvested. Its `called_verdict`
column is `echinocandin_verdict(call=result)` on each isolate's live caller output (one of
`RESISTANCE_MARKER_DETECTED` / `UNCHARACTERIZED_VARIANT` / `NO_KNOWN_MARKER` / `UNRESOLVED`), and its
`susceptibility` is copied verbatim from the frozen concordance fixture — so it satisfies the #137
prereg's *"filled by the orchestration, NOT hand-authored"* rule. The row-assembly logic is
unit-tested offline in `tests/test_accuracy_harvest.py`; only its execution needs this host.

## B. Freeze the harvested #137 fixture

```bash
shasum -a 256 data/earlywarning/recaller_sanity/fks1_accuracy_98.tsv \
  >> work/PREREG_diagnostic_accuracy.sha256
```

(The grader prints these exact freeze + run commands at the end of step A.)

## C. #137 accuracy run (Track 2)

```bash
python scripts/validate_diagnostic_accuracy.py --fixture data/earlywarning/recaller_sanity/fks1_accuracy_98.tsv
```

Reports VME (phenotype R, verdict `NO_KNOWN_MARKER`) and ME (phenotype S, verdict
`RESISTANCE_MARKER_DETECTED`) with CIs, per the frozen bar. UNRESOLVED isolates are excluded, never
counted as a negative. Read the echinocandin arm's PASS/UNDERPOWERED/named-error outcome into the
diagnostics MVP write-up.

## What this session still does NOT unblock (say so plainly)

- **Azole / ERG11 accuracy arm of #137** — Lockhart 2017 gives only 4 clade reference strains
  (a sanity control, not an error rate). Declared separately **underpowered** until a larger paired
  *C. auris* ERG11 collection (or the pooled *C. albicans* option-C route) is transcribed + pinned.
- **Calibration (#136)** — needs the option-C paired *C. albicans* ERG11 panel, which does not
  exist publicly. Unchanged by this session. [[calibration-track-blocked]]
- **ERG11 prevalence `fill`** — its own cloud run (disk ceiling), see
  `work/RUNBOOK_recaller_run.md`. [[recaller-run-environment]]
