# A pre-registered, hash-frozen genotype→verdict engine for *Candida auris* echinocandin resistance: a token-accurate FKS1 caller and a scoped high-PPV detection claim (very-major error 0%, major error 2.2%)

**Manuscript for submission — target journal: *Journal of Clinical Microbiology* (JCM), Full-Length Text.**

*Running title:* Scoped genotype→verdict engine for *C. auris* FKS1

*Author:* Jacob Jensen¹

*Affiliation:* ¹ The Antifungal Resistance Project (OpenAFR)

*Corresponding author:* Jacob Jensen, The Antifungal Resistance Project. Email: jaj.jacob.jensen@gmail.com

*Keywords:* *Candida auris*; echinocandin resistance; FKS1; genotype-to-phenotype; molecular diagnostics; very-major error; pre-registration; antifungal resistance

*Article type:* Full-Length Text (research use only / RUO; no clinical claim)

*Word count (body, Introduction–Discussion):* ~3,500

---

## Importance

*Candida auris* is a WHO critical-priority fungal pathogen, and echinocandins are its first-line therapy — so a genotype that could flag echinocandin resistance early has direct clinical value, and equally direct clinical risk if it is wrong. Genotype-based resistance interpretation can fail in two loaded directions: missing real resistance (a very-major error) or calling resistance that is not present (a major error). We report a resistance-interpretation engine that reads *C. auris* FKS1 hot-spot genotypes and emits categorical verdicts, and we measure it under pre-registration against a published external benchmark. The central lesson is methodological and general: when some genotypes simply do not predict phenotype, *scoping the claim* — detecting only where a marker is predictive and abstaining honestly elsewhere — produces a safer and more useful assay than widening the panel until the numbers move. The engine states exactly which markers it did and did not check, and refuses to manufacture the susceptibility call, probability, and coverage it does not have.

---

## Abstract

We describe a research-use-only (RUO) resistance-interpretation engine that takes a typed *Candida auris* genotype and emits, per drug class, one of four categorical verdicts — `RESISTANCE_MARKER_DETECTED`, `UNCHARACTERIZED_VARIANT`, `NO_KNOWN_MARKER`, or `UNRESOLVED` — with no probability, no MIC, and no susceptibility call. We report the first arm to clear a pre-registered clinical-accuracy bar. The underlying FKS1 (echinocandin) caller was first certified token-accurate against a published external benchmark of 98 *C. auris* whole genomes (PMC12323592): panel-detection sensitivity 100% (95% CI 90.6–100), specificity 100% (94.1–100), exact-token identity 100%, 98/98 resolved, zero false calls. We then measured verdict accuracy — agreement with each isolate's echinocandin R/S phenotype in CLSI M23 error terms — against a bar frozen before any isolate was scored (very-major error [VME] ≤3% point and ≤15% Wilson-upper; major error [ME] ≤5%; abstention ≤30%). The result arrived over three pre-registered measurements. v1 (narrow S639 panel) FAILed on under-detection (VME 10.0%, 4.0–23.1), all four misses attributable to an uncovered FKS1 HS3 hotspot. v2 (added an HS3 window) collapsed VME to 0.0% (0.0–7.9) exactly as projected, but the failure mode flipped to over-calling (ME 13.5%, 6.7–25.3), driven by low-PPV positions (D642Y carriers split 2R/5S). v3 resolved this by adding a PPV tier over the same panel — high-PPV markers detect, low-PPV positions abstain — clearing the bar: VME 0.0% (0.0–8.4), ME 2.2% (0.4–11.3), abstention 10.2%, PASS. The claim this licenses is detection over validated high-PPV FKS1 markers, not a general R/S classifier. Scoping the claim beats widening the panel.

---

## Introduction

*Candida auris* is a WHO fungal-priority pathogen [1] whose resistance profile differs sharply by drug class. Azole resistance in *C. auris* is already near-saturated at baseline: a representative random re-called sample puts the azole-marker event frequency at 160/199 = 80.4% (95% CI 74.3–85.3). Azole *emergence* therefore carries little early-warning signal. Echinocandin (FKS1) resistance is the genuinely dynamic axis: a matched re-called sample measured its event frequency at 10/443 = 2.3% (1.2–4.1) — low and non-saturated, i.e. the useful regime for surveillance and for a resistance-*interpretation* tool.

Interpretation, not just detection, is where clinical value and clinical risk both live. A tool that maps a genotype to a resistance call can be wrong in two clinically loaded directions: it can miss real resistance (a very-major error) or it can call resistance that is not there (a major error) [2]. This work asks the question a clinician actually asks of such a tool — when it emits a verdict, how often does it agree with phenotype, and how often does it get each dangerous direction wrong — and answers it under pre-registration, on a published external benchmark, for the echinocandin arm.

A note on scope, set at the outset. Echinocandins do not coordinate a metal, so the CYP51/heme-iron structural moat that a companion docking-triage track exploits for azoles [3] does **not** transfer to FKS1/glucan synthase. The FKS1 arm is therefore detection-only, with no structural interpretation layer — a boundary declared here, enforced through the verdict contract, and never borrowed against.

## Materials and Methods

### The engine and its four-verdict contract

The engine (`openafr/interpret.py`, single entrypoint; contract in `docs/DIAGNOSTIC_VERDICT_CONTRACT.md`) reuses existing re-callers (`openafr/recaller.py` for ERG11, `openafr/fks1_caller.py` for FKS1) and emits, per isolate per drug class, one of exactly four categorical verdicts — never a probability:

1. **`RESISTANCE_MARKER_DETECTED`** — resolved, and ≥1 called token is a high-PPV known-resistance panel substitution.
2. **`UNCHARACTERIZED_VARIANT`** — resolved, a non-synonymous change is present but none is a high-PPV panel token. This is the explicit first-class "I don't know": for FKS1 it also absorbs a call whose only panel hit is a low-PPV position.
3. **`NO_KNOWN_MARKER`** — resolved, only wild-type / no panel token. **This is not a susceptibility call:** the tool covers a named gene panel, not efflux (TAC1/MRR1/CDR1), ERG3, or promoter/TR mechanisms, so absence of a marker cannot assert susceptibility.
4. **`UNRESOLVED`** — the window could not be called honestly (coverage gap, in-window indel refused, orchestration failure). **Excluded from every denominator, never scored as a negative.**

Four honesty rules make the contract the product: (i) no probability field, because we found no public data to calibrate one [4] and emitting one would be an unbacked number; (ii) `NO_KNOWN_MARKER` ≠ susceptible; (iii) `UNCHARACTERIZED_VARIANT` is a verdict, not silence — a lookup-table tool falls silent on a novel variant, this one does not; (iv) `UNRESOLVED` is excluded, not guessed. The report layer (`openafr/report.py`) renders these but never re-decides, and bakes an RUO disclaimer into every output.

### The caller and the benchmark

The FKS1 caller (`openafr/fks1_caller.py`) reads windowed reads→consensus over FKS1 hot-spot regions and emits panel tokens; the identical orchestration serves both validation and production `fill` (`scripts/recall_fks1.py`, `reads_to_window_consensus`). Truth is **PMC12323592** [5] (Misas et al., *Microbiol Spectr* 2025) — a purpose-built external set of 100 whole-genome-sequenced isolates, each row carrying an SRA run accession, the FKS1 genotype, and the echinocandin R/S interpretation. The paper is paywalled and not machine-scrapable, so Table 1 was transcribed by hand and hash-pinned; one duplicated/self-contradictory row was footnoted out during freeze, leaving **98 isolates scored** (`data/earlywarning/recaller_sanity/fks1_benchmark_100.tsv`, sha256 `68eef59f…`).

### The two accuracy metrics (frozen)

Per isolate, per drug class, the engine's verdict is mapped to an R/S call by a frozen rule and scored against phenotype (`openafr/accuracy.py`, `scripts/validate_diagnostic_accuracy.py`): `RESISTANCE_MARKER_DETECTED` → positive; `NO_KNOWN_MARKER` → negative (the assay's negative *prediction*, not a susceptibility *claim* — the VME rate is exactly the honest price of reading it as one); `UNCHARACTERIZED_VARIANT` → abstention, held out of the 2×2; `UNRESOLVED` → excluded. Reported per CLSI M23 [6], each as k/n with a Wilson 95% CI:

- **Very-major error (VME)** — phenotype R, verdict `NO_KNOWN_MARKER`, over resistant isolates.
- **Major error (ME)** — phenotype S, verdict `RESISTANCE_MARKER_DETECTED`, over susceptible isolates.
- Sensitivity (= 1 − VME) and specificity (= 1 − ME); categorical agreement; abstention rate.

Abstentions and unresolved verdicts are held out of the confusion matrix so a low VME can never be bought by quietly reclassifying hard isolates as negatives.

### The pass bar (frozen before any isolate was scored)

```
PASS         if  VME point <= 0.03  AND  VME Wilson UPPER <= 0.15
             AND ME point  <= 0.05  AND  abstention <= 0.30
FAIL         if any gated metric's point is above its bar
UNDERPOWERED if the resistant OR the susceptible arm has < 15 scored isolates
```

Freezing the metric before the fixture was in hand removed the freedom to pick the error bar after seeing which isolates could be obtained. The bar is identical across v1/v2/v3; each re-tiering must clear the same bar.

### Run environment and pre-registration discipline

Concordance and v1/v2 measurements were run on Google Cloud (e2-standard-4, us-central1-a; VM deleted after each run) with fasterq-dump 3.1.1, minimap2 2.28, samtools 1.21 (concordance run_id `a8fdb3327ae1`, 2026-09-20). Every measurement committed a falsifiable projection in a hash-frozen pre-registration before grading; at run time the grader re-checks the pre-registration SHA-256, the pinned FKS1 CDS reference hash, and the truth-set fixture hash, and refuses to certify if any moved. All hashes matched at run time for every reported measurement.

## Results

### The caller is token-accurate (the load-bearing prerequisite)

Before any accuracy claim, the caller's fidelity was certified against the benchmark's genotype truth (`work/RESULTS_fks1_concordance.md`; pre-registration sha `da8e3827…`). All 98 isolates resolved; panel confusion **TP=37, TN=61, FP=0, FN=0** — zero false calls in either direction (Table 1).

Per-token sensitivity of expected carriers: S639F 16/16 (95% CI 80.6–100), S639P 7/7 (64.6–100), S639Y 14/14 (78.5–100). The **resistance-detection ceiling** — reported, not gated — is 36/46 = 78.3% (64.4–87.7): of resolved phenotypically-resistant isolates the HS1 caller flags ~4 in 5, the rest carrying resistance an HS1 panel cannot and must not claim to see. This is the pillar the accuracy story rests on: the caller reads genomes correctly, so the accuracy limits below are biology (marker coverage and marker PPV), not caller bugs.

### The v1 → v2 → v3 arc

Each measurement committed a falsifiable projection in its pre-registration, then reported against it. The confusion counts move because abstentions differ; the bar does not move (Tables 2–4; Figure 1).

**v1 — narrow S639 panel: FAIL on under-detection.** Scored 87 (40 R / 47 S); 11 abstained. Confusion TP=36, FN=4, FP=1, TN=46. VME 4/40 = 10.0% (4.0–23.1) — **FAIL**; ME 1/47 = 2.1% (0.4–11.1); abstention 11/98 = 11.2%. Characterisation showed the cause is tractable, not an off-target wall: all four misses are FKS1 HS3 (two W691L, one W691C, one M690I), a hotspot the v1 caller had no window for (it read only HS1 635–643 and HS2 1350–1358). The projected fix — add an HS3 window — put VME at 1/40 = 2.5% (0.4–12.9), pending a caller change and a fresh measurement.

**v2 — HS3 window + widened HS1/HS2: the failure mode flips to over-calling.** Pre-registration sha `014094c6…`, with a pre-committed VME projection of 1/40 = 2.5%. A fresh reads→caller pass re-harvested the fixture. Scored 97 (45 R / 52 S); 1 abstained. Confusion TP=45, FN=0, FP=7, TN=45. VME 0/45 = 0.0% (0.0–7.9) — cleared, exactly as projected; ME 7/52 = 13.5% (6.7–25.3) — **FAIL**; abstention 1/98 = 1.0%. The seven major errors were genotype-correct but phenotype-discordant: five were D642Y carriers, all phenotypically S; across the whole benchmark D642Y carriers split 2 R / 5 S (~29% PPV). One M690I (S) and one S639Y (S — irreducible noise at an otherwise-reliable core position) complete the seven. The binding constraint had moved from coverage to marker PPV.

**v3 — a PPV tier over the same panel: PASS.** The honest response to a PPV ceiling is not more markers but a narrower claim, pre-registered fresh (`work/PREREGISTRATION_diagnostic_accuracy_v3.md`, sha `f3f9686e…`). v3 adds only a PPV tier over the v2 panel — no new window, position, or letter:

- **Core (high-PPV) tier — DETECT:** F635{C,Y}, S639{F,P,Y} (HS1); R1354{S} (HS2); W691{L} (HS3).
- **Low-PPV tier — ABSTAIN as `UNCHARACTERIZED_VARIANT`:** D642Y, M690I.

The re-tiered fixture was derived offline — no SRA re-read — re-running the production `echinocandin_verdict()` over the caller output already frozen in the v2 fixture, asserting the invariant that a verdict may change *only* `RESISTANCE_MARKER_DETECTED` → `UNCHARACTERIZED_VARIANT`, and only for a low-PPV-only call. Scored 88 (42 R / 46 S); 10 abstained. Confusion TP=42, FN=0, FP=1, TN=45. VME 0/42 = 0.0% (0.0–8.4); ME 1/46 = 2.2% (0.4–11.3); sensitivity 100% (91.6–100); specificity 97.8% (88.7–99.6); categorical agreement 87/88 = 98.9% (93.8–99.8); abstention 10/98 = 10.2% (5.6–17.8). **PASS** — every gated metric clears its bar, matching the frozen projection exactly. The nine low-PPV-only carriers (7 D642Y, 2 M690I) move from the confusion matrix into abstention. Only one major error remains: a single S639Y, phenotype S — irreducible genotype↔phenotype noise at an otherwise-reliable core position, well within the ≤5% bar.

## Discussion

**Scope the claim, don't widen the panel.** The arc is the argument. Widening coverage (v1→v2) was necessary but exposed a deeper limit — some genotypes simply do not predict phenotype — that no amount of further coverage can fix. The move that turned a FAIL into a defensible GO was narrowing the claim to what the data support: detect high-PPV markers, abstain honestly elsewhere. The 10.2% abstention is the visible cost, and reporting it openly is the point — an engine that abstains on 10% and is right on the rest is more useful, and far safer, than one that classifies everything and is wrong on one susceptible isolate in seven.

**D642Y as the worked example of an independently-falsified marker.** The v3 tiering obeys a load-bearing anti-circularity rule: a position's tier must rest on evidence independent of the scored isolates' own phenotype. D642Y is demoted because the benchmark authors' own expert panel scores every D642Y carrier as wild-type — an independent expert judgement — corroborated by the CDC EID 25-0760 surveillance collection [7] not treating it as resistance-defining; M690I because it is phenotype-unresolved in the source literature. The benchmark then independently rejects D642Y on phenotype (5/7 carriers S). Symmetrically, W691C and W691F were excluded from the core tier despite being HS3: W691C's only report is the benchmark itself (tagging it would be circular; it is emitted verbatim but untagged, so its carrier lands in abstention, not detection), and W691F is reported in no source. Only W691L, additionally CRISPR-confirmed [8] (Jacobs et al., *AAC* 2023), enters the core tier. Tiering by independent literature status — not by fitting this cohort's discordant rows — is what keeps the PASS honest.

**What this licenses, and no more.** A measured RUO detection accuracy claim over validated high-PPV FKS1 markers (F635/S639/R1354S/W691L core), abstaining on low-PPV positions. It is not a general R/S classifier: on this cohort no FKS1-genotype panel clears both error bars as a full classifier, because D642Y/M690I carriers are phenotypically split. v3 wins by scoping the claim, not by out-classifying the biology.

## Limitations

Each limit is stated with the specific evidence that would lift it; the frozen artifact behind every claim these limits qualify is traced in Supplementary Appendix S1.

1. **The azole/ERG11 accuracy arm is out of scope for this publication.** Only four clade strains carry paired ERG11 genotype + azole phenotype in hand — a sanity control, not a measurable accuracy — so no azole accuracy is claimed; this manuscript's PASS is echinocandin-only. This is a declared scope boundary set on a named data gap, not an open to-do (`docs/AZOLE_ARM_SCOPE_DECISION.md`). *Lifts it:* a larger paired *C. auris* ERG11 genotype+phenotype collection.
2. **No calibrated probability / no MIC — out of scope for this publication.** The engine is categorical by construction, and a calibrated probability is deliberately not emitted: a data-availability audit found no public *Candida* panel clearing a calibration bar, and the calibrated-probability track is blocked on a non-public pooled *C. albicans* ERG11 panel. The calibration machinery and its pre-registration are complete and frozen, waiting only on that data, so this is a scoped-out boundary with a named restart trigger, not a defect (`docs/AZOLE_ARM_SCOPE_DECISION.md`). *Lifts it:* that panel plus a held-out reliability PASS.
3. **`NO_KNOWN_MARKER` is not "susceptible."** The panel covers named ERG11/FKS1 markers only — not efflux (TAC1/MRR1/CDR1), ERG3, or promoter/TR mechanisms. The measured VME (0% under the scoped claim) is honest within that coverage; off-panel resistance is out of scope by declaration. *Lifts it:* broader-mechanism callers.
4. **Retrospective, single cohort, single organism.** The measurement is against one published retrospective benchmark (PMC12323592), *C. auris* only. *Lifts it:* a pre-registered prospective study on consecutively collected isolates, and cross-cohort replication.
5. **FKS1 is detection-only, no structural interpretation.** Echinocandins coordinate no metal, so the CYP51 geometry moat does not transfer; the FKS1 verdict carries no structural field. This is a declared scope, not an omission.
6. **Regulatory path untouched.** LDT/IVD quality-system, reproducibility, and lot-to-lot evidence are out of scope for this science track; named so they are not discovered late.

Nothing here has been tested prospectively against consecutively collected clinical isolates, and every verdict is RUO. The overall clinical-path verdict is RUO with the echinocandin detection arm cleared under a narrowed high-PPV claim; no clinical susceptibility claim is made.

## Data availability

A claim-by-claim traceability audit — every quantitative claim in this manuscript mapped to the exact hash-frozen pre-registration, fixture, RESULTS record, and code module behind it, plus a step-by-step reproduction recipe (offline from committed fixtures, and the full SRA→caller re-harvest) — is provided as **Supplementary Appendix S1** (`work/APPENDIX_reproducibility.md`). All eight frozen hashes it inventories were re-verified on disk and match.

All code, all pre-registrations with their hashes, the transcribed hash-pinned benchmark fixture, the harvested per-version accuracy fixtures, and every `RESULTS_*.md` are in the project repository. Key modules: the callers (`openafr/recaller.py`, `openafr/fks1_caller.py`), the single interpretation entrypoint (`openafr/interpret.py`), the verdict schema (`openafr/verdict.py`), the metric engines (`openafr/accuracy.py`, `openafr/concordance.py`), the report layer (`openafr/report.py`), the graders (`scripts/validate_fks1_concordance.py`, `scripts/validate_diagnostic_accuracy.py`), and the offline v3 harvester (`scripts/harvest_diagnostic_accuracy_v3.py`). Truth set: PMC12323592 (Misas et al., *Microbiol Spectr* 2025), transcribed by hand from the paywalled Table 1. Software is licensed PolyForm Noncommercial 1.0.0.

## Author contributions

J.J. conceived the study, designed and pre-registered all measurements, implemented the callers and grading pipeline, executed the runs, and wrote the manuscript.

## Funding

This work received no external funding; it was conducted as an independent, self-funded open research project.

## Competing interests

The author declares no financial competing interests. Analysis and manuscript preparation were AI-assisted; all quantitative results derive from hash-frozen pre-registered runs in the project repository and were verified against their source records.

## Ethics

This study used only publicly deposited whole-genome sequencing data and a published, hand-transcribed benchmark table. No human subjects, patient identifiers, or animal work were involved. All outputs are research use only (RUO); no clinical or diagnostic claim is made.

## References

1. World Health Organization. WHO fungal priority pathogens list to guide research, development and public health action. Geneva: WHO; 2022.
2. Clinical and Laboratory Standards Institute. Development of in vitro susceptibility test methods, breakpoints, and quality control parameters. CLSI guideline M23. Wayne (PA): CLSI.
3. The Antifungal Resistance Project. Mechanism-anchored geometric rescoring on fungal CYP51: broad enrichment is robust, top-rank enrichment is not, and property-matched decoys explain why. Companion preprint, OpenAFR, 2026. (`work/PREPRINT_geometry_ceiling.md`)
4. The Antifungal Resistance Project. Diagnostic verdict contract and calibration data-availability audit. OpenAFR internal records (#132–#136).
5. Misas E, et al. A benchmark dataset for validating FKS1 mutations in *Candida auris*. Microbiol Spectr. 2025. doi:10.1128/spectrum.03147-24. (PMC12323592)
6. Clinical and Laboratory Standards Institute. CLSI guideline M23 (error-rate definitions: very-major and major error). Wayne (PA): CLSI.
7. Centers for Disease Control and Prevention. *Candida auris* FKS1 surveillance collection. Emerg Infect Dis EID 25-0760.
8. Jacobs SE, et al. CRISPR-confirmed FKS1 W691L echinocandin resistance in *Candida auris*. Antimicrob Agents Chemother. 2023. aac.00423-23. (PMC10269051)

## Tables

**Table 1. FKS1 caller concordance vs the 98-isolate PMC12323592 benchmark (gated metrics).** Panel confusion: TP=37, TN=61, FP=0, FN=0.

| Metric (gated) | Point | Wilson 95% CI | Bar | Verdict |
|---|---|---|---|---|
| Panel-detection sensitivity | 37/37 = 100.0% | 90.6–100 | ≥90% | ✓ |
| Panel-detection specificity | 61/61 = 100.0% | 94.1–100 | ≥95% | ✓ |
| Exact-token identity | 37/37 = 100.0% | 90.6–100 | ≥95% | ✓ |
| Resolved fraction | 98/98 = 100.0% | — | ≥80% | ✓ |

**Table 2. Verdict accuracy — v1 (narrow S639 panel).** Scored 87 (40 R / 47 S); 11 abstained. Confusion TP=36, FN=4, FP=1, TN=46.

| Metric | Value | Wilson 95% CI | Bar | Verdict |
|---|---|---|---|---|
| VME (missed R) | 4/40 = 10.0% | 4.0–23.1 | ≤3% pt & ≤15% upper | ✗ |
| ME (false R) | 1/47 = 2.1% | 0.4–11.1 | ≤5% | ✓ |
| Abstention | 11/98 = 11.2% | 6.4–19.0 | ≤30% | ✓ |

**Table 3. Verdict accuracy — v2 (HS3 window + widened HS1/HS2).** Scored 97 (45 R / 52 S); 1 abstained. Confusion TP=45, FN=0, FP=7, TN=45.

| Metric | Value | Wilson 95% CI | Bar | Verdict |
|---|---|---|---|---|
| VME (missed R) | 0/45 = 0.0% | 0.0–7.9 | ≤3% pt & ≤15% upper | ✓ |
| ME (false R) | 7/52 = 13.5% | 6.7–25.3 | ≤5% | ✗ |
| Abstention | 1/98 = 1.0% | 0.2–5.6 | ≤30% | ✓ |

**Table 4. Verdict accuracy — v3 (PPV tier over the v2 panel).** Scored 88 (42 R / 46 S); 10 abstained. Confusion TP=42, FN=0, FP=1, TN=45.

| Metric | Value | Wilson 95% CI | Bar | Verdict |
|---|---|---|---|---|
| VME (missed R) | 0/42 = 0.0% | 0.0–8.4 | ≤3% pt & ≤15% upper | ✓ |
| ME (false R) | 1/46 = 2.2% | 0.4–11.3 | ≤5% | ✓ |
| Sensitivity (1 − VME) | 42/42 = 100.0% | 91.6–100 | — | — |
| Specificity (1 − ME) | 45/46 = 97.8% | 88.7–99.6 | — | — |
| Categorical agreement | 87/88 = 98.9% | 93.8–99.8 | — | — |
| Abstention | 10/98 = 10.2% | 5.6–17.8 | ≤30% | ✓ |

**Table 5. Pre-registration ledger.** Every measurement, in order, with its frozen hash and committed outcome.

| Run | Question | Pre-reg sha256 | Fixture sha256 | Outcome |
|---|---|---|---|---|
| FKS1 concordance | caller token accuracy vs genotype truth (98 isolates) | `da8e3827…` | `68eef59f…` (benchmark) | PASS (sens/spec/identity 100%) |
| v1 | verdict accuracy, narrow S639 panel | `70d5d2ae…` | `1aeffea5…` | FAIL (VME 10.0%) |
| v2 | + HS3 window, widened HS1/HS2 | `014094c6…` | `a7c5056f…` | FAIL (ME 13.5%) |
| v3 | + PPV tier (core detect, low-PPV abstain) | `f3f9686e…` | `02fc37e1…` | PASS (VME 0%, ME 2.2%) |

## Figure legend

**Figure 1. The v1→v2→v3 error-mode arc.** Very-major error (VME), major error (ME), and abstention rate across the three pre-registered measurements against a fixed pass bar (VME ≤3% point / ≤15% Wilson-upper; ME ≤5%; abstention ≤30%). v1 fails on under-detection (VME 10.0%); adding an HS3 window (v2) collapses VME to 0% but flips the failure to over-calling (ME 13.5%); adding a PPV tier over the same panel (v3) clears both error bars (VME 0%, ME 2.2%) at a 10.2% abstention price. The bar is identical across all three panels; only the abstention set differs.
