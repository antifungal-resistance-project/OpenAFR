# Supplementary Appendix S1 — Claims→artifacts audit and reproducibility

**Companion to** `work/MANUSCRIPT_fks1_diagnostics_JCM.md` (JCM Full-Length Text). The same
frozen artifacts back the companion preprint `work/PREPRINT_diagnostics_engine.md`; this
appendix is written against the manuscript's claims and applies to both.

**Purpose.** Every quantitative claim in the manuscript resolves here to the exact
hash-frozen artifact behind it, and every number is reproducible from committed inputs by
the recipe in §S1.3. This is the integrity pass for track-A (#161): no claim without an
artifact, no artifact without a hash, no number without a way to re-derive it.

**Audit status (verified 2026-09-30).** All 8 frozen hashes cited by the manuscript
(4 pre-registrations, 4 fixtures) were recomputed on disk and **all 8 match**. No claim in
the manuscript lacks a backing artifact; no number required correction.

---

## S1.1 Claims→artifacts traceability matrix

One row per quantitative claim. "RESULTS" is the human-readable record; "Pre-reg" is the
falsifiable projection frozen *before* the measurement; "Fixture" is the scored data; hashes
are 8-char prefixes of the sha256 in §S1.2. Code modules are the implementations the number
is produced by.

### Caller concordance (the load-bearing prerequisite)

| Claim | Manuscript locus | RESULTS | Pre-reg (sha) | Fixture (sha) | Code |
|---|---|---|---|---|---|
| Panel confusion TP=37, TN=61, FP=0, FN=0; 98/98 resolved | Abstract, Results §"token-accurate", Table 1 | `RESULTS_fks1_concordance.md` | `da8e3827` | `68eef59f` (benchmark) | `concordance.py`, `fks1_caller.py` |
| Sensitivity 100% (90.6–100), specificity 100% (94.1–100), exact-token identity 100% | Abstract, Table 1 | `RESULTS_fks1_concordance.md` | `da8e3827` | `68eef59f` | `concordance.py` |
| Per-token sensitivity S639F 16/16, S639P 7/7, S639Y 14/14 | Results §"token-accurate" | `RESULTS_fks1_concordance.md` | `da8e3827` | `68eef59f` | `concordance.py` |
| Resistance-detection ceiling 36/46 = 78.3% (64.4–87.7) — reported, not gated | Results §"token-accurate" | `RESULTS_fks1_concordance.md` | `da8e3827` | `68eef59f` | `concordance.py` |

### Verdict accuracy (v1 → v2 → v3)

| Claim | Manuscript locus | RESULTS | Pre-reg (sha) | Fixture (sha) | Code |
|---|---|---|---|---|---|
| v1 (narrow S639): VME 4/40 = 10.0% (4.0–23.1) FAIL; ME 1/47 = 2.1%; abst 11/98; confusion TP=36 FN=4 FP=1 TN=46 | Results §v1, Table 2 | `RESULTS_diagnostic_accuracy.md` | `70d5d2ae` | `1aeffea5` | `interpret.py`, `accuracy.py` |
| v2 (+HS3): VME 0/45 = 0.0% (0.0–7.9); ME 7/52 = 13.5% (6.7–25.3) FAIL; abst 1/98; confusion TP=45 FN=0 FP=7 TN=45 | Results §v2, Table 3 | `RESULTS_diagnostic_accuracy.md` | `014094c6` | `a7c5056f` | `interpret.py`, `accuracy.py` |
| v3 (+PPV tier): VME 0/42 = 0.0% (0.0–8.4); ME 1/46 = 2.2% (0.4–11.3); sens 100%; spec 97.8%; cat. agreement 87/88 = 98.9%; abst 10/98 = 10.2% PASS; confusion TP=42 FN=0 FP=1 TN=45 | Abstract, Results §v3, Table 4 | `RESULTS_diagnostic_accuracy.md` | `f3f9686e` | `02fc37e1` | `interpret.py`, `accuracy.py`, `harvest_diagnostic_accuracy_v3.py` |
| Pass bar: VME pt ≤3% & upper ≤15%; ME ≤5%; abst ≤30% — identical across v1/v2/v3 | Methods §"pass bar", Tables 2–4 | (frozen in each pre-reg) | `70d5d2ae`/`014094c6`/`f3f9686e` | — | `accuracy.py` |

### Prevalence framing (Introduction)

| Claim | Manuscript locus | RESULTS | Pre-reg (sha) | Fixture (sha) | Code |
|---|---|---|---|---|---|
| Azole marker event frequency 160/199 = 80.4% (74.3–85.3) | Introduction | `RESULTS_prevalence.md` | — (descriptive snapshot) | — | `recaller.py` |
| Echinocandin event frequency 10/443 = 2.3% (1.2–4.1) | Introduction | `RESULTS_fks1_prevalence.md` | — (descriptive snapshot) | — | `fks1_caller.py` |

### External-literature claims (not repo artifacts)

| Claim | Manuscript locus | Source | Nature |
|---|---|---|---|
| Benchmark = 98 scored isolates with genotype + echinocandin R/S | Methods, throughout | PMC12323592 (Misas et al., *Microbiol Spectr* 2025) — hand-transcribed to fixture `68eef59f` | external ref + transcribed fixture |
| D642Y treated as wild-type by benchmark expert panel; not resistance-defining in surveillance | Discussion | PMC12323592 expert panel; CDC EID 25-0760 | external ref |
| D642Y carriers split 2 R / 5 S (~29% PPV) on this cohort | Results §v2, Discussion | derived from fixture `a7c5056f` | repo-derived from external fixture |
| W691L CRISPR-confirmed as resistance-causing | Discussion | Jacobs et al., *AAC* 2023 (PMC10269051) | external ref |

---

## S1.2 Frozen-artifact inventory + verification

The frozen surface is pinned by six `work/PREREG_*.sha256` files — each pins both its
pre-registration markdown **and** the fixture that measurement scored. The standing drift
gate `scripts/repro.py` (CI: `.github/workflows/repro.yml`) verifies every one on each PR
and fails loudly (exit 1) on any drift.

| Artifact | sha256 | Pinned by |
|---|---|---|
| `work/PREREGISTRATION_fks1_concordance.md` | `da8e38276229cc838b467c7a1915529f9c413d40e7f4bfd1e8d4e8711b92fb22` | `PREREG_fks1_concordance.sha256` |
| `data/earlywarning/recaller_sanity/fks1_benchmark_100.tsv` | `68eef59f2c82e1133efac82331706cdaccb4f44c36771002e281e44f2b76ade3` | `PREREG_fks1_concordance.sha256` |
| `work/PREREGISTRATION_diagnostic_accuracy.md` | `70d5d2aeb1d7863998829766a7ed473e5f8722b157dcb312cce4603fe253fa93` | `PREREG_diagnostic_accuracy.sha256` |
| `data/earlywarning/recaller_sanity/fks1_accuracy_98.tsv` | `1aeffea5834430b3505b9aad97a37bcc30ccc663f2a2f1b3ec16ec4c4bb71029` | `PREREG_diagnostic_accuracy.sha256` |
| `work/PREREGISTRATION_diagnostic_accuracy_v2.md` | `014094c654a75fcb5853397bbb1484f1803a053bbbd904f70eac48c9801e85e7` | `PREREG_diagnostic_accuracy_v2.sha256` |
| `data/earlywarning/recaller_sanity/fks1_accuracy_98_v2.tsv` | `a7c5056f8528d711bd0398dc27e64871d555dcc8953f19e8df71889f6a26668d` | `PREREG_diagnostic_accuracy_v2.sha256` |
| `work/PREREGISTRATION_diagnostic_accuracy_v3.md` | `f3f9686efd63d3df47d8e5fae9f3a8f7cc2bd7abbfb14909eac3b5743dd7430a` | `PREREG_diagnostic_accuracy_v3.sha256` |
| `data/earlywarning/recaller_sanity/fks1_accuracy_98_v3.tsv` | `02fc37e10dadc88f8f2dc5528f75df080ad8f5ef1150d1f62a79b34a07fd9f15` | `PREREG_diagnostic_accuracy_v3.sha256` |

Re-verify all eight at once (from the repo root):

```bash
shasum -a 256 -c work/PREREG_fks1_concordance.sha256 \
                 work/PREREG_diagnostic_accuracy.sha256 \
                 work/PREREG_diagnostic_accuracy_v2.sha256 \
                 work/PREREG_diagnostic_accuracy_v3.sha256
```

Or verify the whole frozen surface (protocol + every pre-registration + downstream
derivations) with the standing gate:

```bash
python scripts/repro.py        # exit 0 = no drift
```

All eight confirmed matching on 2026-09-30.

---

## S1.3 How to reproduce each number

Two tiers. Tier A reproduces the **verdict-accuracy** numbers (Tables 2–4) offline from
committed fixtures — no network, no SRA, minutes on a laptop — because each accuracy
fixture already carries the caller's per-isolate verdict. Tier B is the live
reads→caller→score path: it produces the **caller-concordance** numbers (Table 1) and
regenerates the accuracy fixtures from raw reads, and needs SRA access plus the aligner
stack. Table 1's grader has no offline fixture-scoring mode by design — it re-reads the
genomes every run, which is exactly the fidelity it certifies.

### Tier A — offline, from committed fixtures (no SRA pull)

Tables 2–4 reproduce without re-reading a single genome, because the caller output is
frozen in the accuracy fixtures above.

```bash
# Verdict accuracy v1 / v2 / v3 → Tables 2 / 3 / 4.
# The version is selected by the --fixture path (the grader verifies that fixture's hash
# against the union of the three PREREG_diagnostic_accuracy*.sha256 pins files).
python scripts/validate_diagnostic_accuracy.py --fixture data/earlywarning/recaller_sanity/fks1_accuracy_98.tsv     # v1
python scripts/validate_diagnostic_accuracy.py --fixture data/earlywarning/recaller_sanity/fks1_accuracy_98_v2.tsv  # v2
python scripts/validate_diagnostic_accuracy.py --fixture data/earlywarning/recaller_sanity/fks1_accuracy_98_v3.tsv  # v3

# v3 PPV re-tiering derivation (no SRA re-read): re-runs echinocandin_verdict() over the
# frozen v2 caller output, asserting verdicts may only move
# RESISTANCE_MARKER_DETECTED → UNCHARACTERIZED_VARIANT for low-PPV-only calls.
python scripts/harvest_diagnostic_accuracy_v3.py
```

Each grader re-checks, at run time, the frozen pre-registration hash(es) and the given
fixture's hash, and refuses to certify if any moved — so a green run is also a
hash-integrity check. (Confirm exact flags with `--help`.)

### Tier B — full SRA → caller re-harvest (regenerates the fixtures)

This is what produced the fixtures and Table 1, and was run on cloud hardware. It pulls
reads from SRA, builds windowed consensus over the FKS1 hot-spots, and calls panel tokens —
the identical orchestration (`reads_to_window_consensus`) used for production `fill`. It
needs `fasterq-dump` (sra-tools), `minimap2`, and `samtools` on PATH.

```bash
# Caller concordance → Table 1 (sens/spec/identity, confusion, 36/46 ceiling).
# Re-reads all 98 benchmark genomes; --limit N for a cheap smoke test.
python scripts/validate_fks1_concordance.py

# Regenerate the accuracy fixtures from reads (SRA accessions → caller tokens).
python scripts/recall_fks1.py
```

**Run environment of record** (concordance + v1/v2): Google Cloud `e2-standard-4`,
`us-central1-a`, VM deleted after each run; `fasterq-dump` 3.1.1, `minimap2` 2.28,
`samtools` 1.21; concordance `run_id a8fdb3327ae1`, 2026-09-20
(`data/earlywarning/runlog/concordance-fks1.jsonl`). v3 added no SRA read — it is the
Tier-A offline re-tiering above.

### What is offline- vs. SRA-reproducible

| Number | Tier A (offline) | Tier B (SRA) |
|---|---|---|
| Table 1 — caller concordance | ❌ (grader re-reads genomes by design) | ✅ `validate_fks1_concordance.py` |
| Tables 2–4 — verdict accuracy | ✅ from committed fixtures | ✅ regenerates fixtures first |
| Prevalence 160/199, 10/443 (Introduction) | — (descriptive snapshots) | ✅ via the re-callers over their cohorts |
| Frozen-hash integrity (§S1.2) | ✅ `shasum -c` / `repro.py` | n/a |
