# FKS1 re-caller concordance vs the 100-isolate PMC12323592 benchmark — measured accuracy

**The FKS1 echinocandin re-caller's accuracy, measured at last.** Until now the caller's error
rate was only ever a 4-strain PASS/FAIL smoke test (`scripts/recaller_sanity_fks1.py`); its
measured *prevalence* (2.3%, [`RESULTS_fks1_prevalence.md`](RESULTS_fks1_prevalence.md)) trusted a
caller whose real-reads accuracy was asserted, not quantified. This run closes that gap — the same
"widen n → measured CI" move that retired preprint Limitation #1 for the docking gate, applied to
the caller — by grading it against a **fixed, published, purpose-built external truth set**.

**Pre-registration:** `work/PREREGISTRATION_fks1_concordance.md` (frozen 2026-09-07, sha
`da8e3827…`), pass bar + both-direction interpretation fixed before any isolate beyond the 4-strain
control was scored. The grader re-checked this document's SHA-256, the pinned FKS1 CDS reference
hash, and the truth-set fixture hash, and would refuse to certify if any moved — **all three
matched at run time.**

**Code:** `openafr/fks1_caller.py` (windowed caller) · `scripts/recall_fks1.py`
(`reads_to_window_consensus`, the identical orchestration the real `fill` uses) ·
`openafr/concordance.py` (pinned metric, unit-tested in `tests/test_concordance.py`) ·
`scripts/validate_fks1_concordance.py` (grader).
**Provenance:** `data/earlywarning/runlog/concordance-fks1.jsonl` run_id `a8fdb3327ae1`,
finished 2026-09-20T05:30Z; GCP e2-standard-4, us-central1-a (VM deleted after the run). Tools:
fasterq-dump 3.1.1, minimap2 2.28, samtools 1.21.

## Truth set

**PMC12323592** — *"A benchmark dataset for validating FKS1 mutations in Candida auris"* — Table 1,
transcribed by hand (paywalled; not scrapable) and hash-pinned in
`data/earlywarning/recaller_sanity/fks1_benchmark_100.tsv`
(`work/PREREG_fks1_concordance.sha256`). **98 isolates scored** (the duplicated/self-contradictory
row footnoted out during freeze, PR #148). Each row carries the SRA run accession, the FKS1
genotype, and the echinocandin R/S interpretation.

## Headline — PASS on every gated metric

> **All 98 isolates resolved (0 unresolved).** Panel confusion (resolved): **TP=37 TN=61 FP=0
> FN=0** — zero false calls in either direction.

| Metric (gated) | Point | Wilson 95% | Point bar | Lower bar | Verdict |
|---|---|---|---|---|---|
| Panel-detection **sensitivity** | 37/37 = 100.0% | [90.6%, 100%] | ≥90% | ≥75% | ✓ |
| Panel-detection **specificity** | 61/61 = 100.0% | [94.1%, 100%] | ≥95% | ≥90% | ✓ |
| Exact-token **identity** | 37/37 = 100.0% | [90.6%, 100%] | ≥95% | — | ✓ |
| **Resolved fraction** | 98/98 = 100.0% | — | ≥80% | — | ✓ |

Overall accuracy 98/98 = 100.0% [96.2%, 100%]. **Pre-registered verdict: PASS.**

### Per-token sensitivity (of expected carriers)

| token | detected | Wilson 95% |
|---|---|---|
| S639F | 16/16 = 100.0% | [80.6%, 100%] |
| S639P | 7/7 = 100.0% | [64.6%, 100%] |
| S639Y | 14/14 = 100.0% | [78.5%, 100%] |

### Resistance-detection ceiling — reported, NOT gated

> Of resolved **phenotypically resistant** isolates, the HS1 caller flags **36/46 = 78.3%**,
> Wilson 95% **[64.4%, 87.7%]**.

Per the pre-registration this is a *measured ceiling, not a pass/fail*: it sits below panel
sensitivity **by design**, because resistant isolates exist whose mechanism is not an HS1 S639
substitution (F635, D642, HS2 R1354, HS3 W691, and off-target/efflux mechanisms) — resistance an
HS1 caller cannot and must not claim to see. This number honestly bounds what a detection-only HS1
pipeline can achieve against phenotype, and is the direct input to the #137 verdict-accuracy result.

## What it means

1. **The caller is token-accurate — the prevalence number can be trusted.** On real reads over a
   purpose-built external benchmark, the caller emits the right panel token with **100% sensitivity,
   100% specificity, 100% exact-token identity**, and refused to guess on nothing (0/98 unresolved).
   The 4-strain-sanity-control-only caveat on the prevalence result is **retired**: the caller's
   error rate is now measured and CI-bounded.
2. **Perfect token accuracy is not the same as perfect resistance detection.** The
   resistance-detection ceiling (36/46 = 78.3%) is the honest gap: the caller reads the genome
   correctly, but ~1 in 5 phenotypically-resistant isolates carries resistance the HS1 panel does
   not cover. This is a property of the *panel's biological scope*, not of the caller's fidelity —
   and it is exactly what drives the #137 verdict-accuracy FAIL against a clinical VME bar
   (`work/RESULTS_diagnostic_accuracy.md`, `docs/DIAGNOSTIC_GO_NO_GO.md` C2/Rung A).
3. **Detection-of-record is now certified; a clinical diagnostic is not.** A PASS here licenses
   the FKS1 track as a trustworthy *detection* engine of record. It does **not** promote it to a
   susceptibility diagnostic — that claim is separately governed by #137 and its gate ladder, and
   it does not clear.
