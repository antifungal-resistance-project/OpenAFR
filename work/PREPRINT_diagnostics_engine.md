# A genotype→verdict engine for *Candida auris* echinocandin resistance that scopes its claim to high-PPV markers: measured very-major error 0%, major error 2.2%, at a 10% honest-abstention price

**Status:** preprint / internal writeup, 2026-09-27. Companion to
`work/PREPRINT_geometry_ceiling.md` (the docking-triage track). Research use only (RUO);
no clinical claim. Every quantitative statement below is drawn verbatim from a frozen
`work/RESULTS_*.md` or pre-registration in this repository; nothing here is recomputed.

---

## Abstract

We describe an RUO resistance-interpretation engine that takes a typed *C. auris* genotype
and emits, per drug class, one of four **categorical** verdicts —
`RESISTANCE_MARKER_DETECTED` / `UNCHARACTERIZED_VARIANT` / `NO_KNOWN_MARKER` / `UNRESOLVED` —
with no probability, no MIC, and no susceptibility call. The product is the honesty of that
contract: the engine states exactly which markers it did and did not check, and refuses to
manufacture the parts it cannot back.

We report the first arm of that engine to clear a pre-registered clinical-accuracy bar. The
underlying FKS1 (echinocandin) caller was first certified token-accurate against a published
external benchmark of 98 *C. auris* whole genomes (PMC12323592): panel-detection sensitivity
100% [90.6, 100], specificity 100% [94.1, 100], exact-token identity 100%, 98/98 resolved,
zero false calls. We then measured the *verdict* accuracy — agreement with each isolate's
echinocandin R/S phenotype in CLSI-M23 error terms — against a bar frozen before any isolate
was scored (very-major error ≤3% point and ≤15% Wilson-upper; major error ≤5%; abstention
≤30%; ≥15 isolates per arm).

The result arrived over three pre-registered measurements, each committing a falsifiable
projection in advance. **v1** (narrow S639 panel) FAILed on *under-detection*: very-major
error 10.0% [4.0, 23.1], because 4 of 40 resistant isolates carried resistance the panel had
no window for — all four FKS1 HS3, a tractable coverage gap, not an off-target wall. **v2**
(added an HS3 window, widened HS1/HS2) collapsed the very-major error to 0.0% [0.0, 7.9]
exactly as projected — but the failure mode *flipped* to *over-calling*: major error rose to
13.5% [6.7, 25.3], because two newly-tagged positions (D642Y, M690I) are known FKS1 changes
whose carriers split R/S (D642Y 2 R / 5 S; ~29% PPV) and so do not predict the phenotype.
More markers could not fix this: dropping D642Y re-opens the very-major error; tagging it
over-calls. **v3** resolved it by adding a **PPV tier** over the same panel — high-PPV core
markers (F635/S639/R1354S/W691L) DETECT; low-PPV positions (D642Y/M690I) ABSTAIN
(`UNCHARACTERIZED_VARIANT`, never resistance and never susceptible) — clearing the bar:
very-major error 0/42 = 0.0% [0.0, 8.4], major error 1/46 = 2.2% [0.4, 11.3], abstention
10.2%, PASS, matching its frozen projection exactly.

The central finding is methodological: for a genotype→phenotype engine, *scoping the claim*
(detect where the marker is predictive, abstain where it is not) beats *widening the panel*.
The 10.2% abstention is the honest price, reported openly. The claim this licenses is
detection over validated high-PPV FKS1 markers, not a general R/S classifier — no
FKS1-genotype panel clears both error bars as a full classifier on this cohort. The azole
(ERG11) arm remains underpowered and the calibrated-probability track remains blocked on
non-public data; both are carried as declared limits, not glossed.

---

## 1. Background

*Candida auris* is a WHO priority fungal pathogen. Azole resistance in *C. auris* is already
near-saturated at baseline — a representative random re-called sample puts the azole-marker
event frequency at 160/199 = 80.4% [74.3, 85.3] (`work/RESULTS_prevalence.md`) — so azole
*emergence* carries little early-warning signal. Echinocandin (FKS1) resistance is the
genuinely dynamic axis: a matched re-called sample measured its event frequency at 10/443 =
2.3% [1.2, 4.1] (`work/RESULTS_fks1_prevalence.md`) — low and non-saturated, i.e. the useful
regime for surveillance and for a resistance-*interpretation* tool.

Interpretation, not just detection, is where clinical value and clinical risk both live. A
tool that maps a genotype to a resistance call can be wrong in two clinically loaded
directions: it can miss real resistance (a very-major error) or it can call resistance that
is not there (a major error). This work asks the question a clinician actually asks of such a
tool — *when it emits a verdict, how often does it agree with phenotype, and how often does it
get each dangerous direction wrong* — and answers it under pre-registration, on a published
external benchmark, for the echinocandin arm.

A note on scope, set at the outset. Echinocandins do not coordinate a metal, so the
CYP51/heme-iron structural moat that the docking track exploits for azoles
(`work/PREPRINT_geometry_ceiling.md`) does **not** transfer to FKS1/glucan synthase. The FKS1
arm is therefore *detection-only, no structural so-what* — declared, enforced through the
verdict contract, and never borrowed against.

## 2. The engine and its contract

The engine (`openafr/interpret.py`, single entrypoint; contract in
`docs/DIAGNOSTIC_VERDICT_CONTRACT.md`, #133) reuses the existing re-callers
(`openafr/recaller.py` for ERG11, `openafr/fks1_caller.py` for FKS1) and emits, per isolate
per drug class, one of exactly four **categorical** verdicts — never a probability:

1. **`RESISTANCE_MARKER_DETECTED`** — resolved, and ≥1 called token is a high-PPV
   known-resistance panel substitution.
2. **`UNCHARACTERIZED_VARIANT`** — resolved, a non-synonymous change is present but none is a
   high-PPV panel token. The explicit first-class **"I don't know"**: for ERG11 it carries a
   calibrated-*low* structural best-guess (#135); for FKS1 it now also absorbs a call whose
   only panel hit is a **low-PPV** position (§5).
3. **`NO_KNOWN_MARKER`** — resolved, only wild-type / no panel token. **This is not a
   susceptibility call** (rule 2 below): the tool covers a named gene panel, not efflux
   (TAC1/MRR1/CDR1), ERG3, or promoter/TR mechanisms, so absence of a marker cannot assert S.
4. **`UNRESOLVED`** — the window could not be called honestly (coverage gap, in-window indel
   refused, orchestration failure). **Excluded from every denominator, never scored as a
   negative.**

Four load-bearing honesty rules make the contract the product: (1) no probability field day
one — #132 found no public data to calibrate one, so emitting one would be an unbacked number;
(2) `NO_KNOWN_MARKER` ≠ susceptible; (3) `UNCHARACTERIZED_VARIANT` is a verdict, not silence
(a lookup-table tool falls silent on a novel variant — this one does not); (4) `UNRESOLVED` is
excluded, not guessed. The report layer (`openafr/report.py`) renders these but never
re-decides, and bakes an RUO disclaimer into every output.

## 3. Methods

### 3.1 The caller and the benchmark

The FKS1 caller (`openafr/fks1_caller.py`) reads windowed reads→consensus over FKS1 hot-spot
regions and emits panel tokens; the identical orchestration serves both the validation runs
and production `fill` (`scripts/recall_fks1.py`, `reads_to_window_consensus`). Truth is
**PMC12323592** (Misas et al., *"A benchmark dataset for validating FKS1 mutations in Candida
auris,"* Microbiol Spectr 2025; doi:10.1128/spectrum.03147-24) — a purpose-built external set
of 100 WGS isolates, each row carrying an SRA run accession, the FKS1 genotype, and the
echinocandin R/S interpretation. The paper is paywalled and not scrapable, so Table 1 was
transcribed by hand and hash-pinned; one duplicated/self-contradictory row was footnoted out
during freeze (#148), leaving **98 isolates scored**
(`data/earlywarning/recaller_sanity/fks1_benchmark_100.tsv`, sha `68eef59f…`).

### 3.2 The two accuracy metrics (frozen)

Per isolate, per drug class, the engine's verdict is mapped to an R/S call by a frozen rule
and scored against phenotype (`openafr/accuracy.py`, `scripts/validate_diagnostic_accuracy.py`):
`RESISTANCE_MARKER_DETECTED` → positive; `NO_KNOWN_MARKER` → negative (the assay's negative
*prediction*, not a susceptibility *claim* — the very-major-error rate is exactly the honest
price of reading it as one); `UNCHARACTERIZED_VARIANT` → **abstention, held out of the 2×2**;
`UNRESOLVED` → excluded. Reported, each as k/n with a Wilson 95% CI:

- **Very-major error (VME)** — phenotype R, verdict `NO_KNOWN_MARKER`, over resistant isolates.
- **Major error (ME)** — phenotype S, verdict `RESISTANCE_MARKER_DETECTED`, over susceptible.
- Sensitivity (= 1 − VME) and specificity (= 1 − ME); categorical agreement; abstention rate.

Abstentions and unresolved verdicts are held out of the confusion matrix so a low VME can
never be bought by quietly reclassifying hard isolates as negatives.

### 3.3 The pass bar (frozen before any isolate was scored)

```
PASS         if  VME point <= 0.03  AND  VME Wilson UPPER <= 0.15
             AND ME point  <= 0.05  AND  abstention <= 0.30
FAIL         if any gated metric's point is above its bar
UNDERPOWERED if the resistant OR the susceptible arm has < 15 scored isolates
```

Freezing the metric before the fixture was in hand removed the freedom to pick the error bar
after seeing which isolates could be obtained. The bar is identical across v1/v2/v3; each
re-tiering must clear the *same* bar.

## 4. The caller is token-accurate (the load-bearing prerequisite)

Before any accuracy claim, the caller's *fidelity* was certified against the benchmark's
*genotype* truth (`work/RESULTS_fks1_concordance.md`; prereg sha `da8e3827…`; run_id
`a8fdb3327ae1`, 2026-09-20, GCP e2-standard-4; fasterq-dump 3.1.1 / minimap2 2.28 / samtools
1.21). All 98 isolates resolved; panel confusion **TP=37 TN=61 FP=0 FN=0** — zero false calls
in either direction.

| Metric (gated) | Point | Wilson 95% | Verdict |
|---|---|---|---|
| Panel-detection sensitivity | 37/37 = 100.0% | [90.6%, 100%] | ✓ (bar ≥90%) |
| Panel-detection specificity | 61/61 = 100.0% | [94.1%, 100%] | ✓ (bar ≥95%) |
| Exact-token identity | 37/37 = 100.0% | [90.6%, 100%] | ✓ (bar ≥95%) |
| Resolved fraction | 98/98 = 100.0% | — | ✓ (bar ≥80%) |

Per-token sensitivity: S639F 16/16 [80.6, 100], S639P 7/7 [64.6, 100], S639Y 14/14 [78.5, 100].
The **resistance-detection ceiling** — reported, *not* gated — is 36/46 = 78.3% [64.4, 87.7]:
of resolved phenotypically-resistant isolates the HS1 caller flags ~4 in 5, the rest carrying
resistance an HS1 panel cannot and must not claim to see. This is the pillar the accuracy
story rests on: **the caller reads genomes correctly; the accuracy limits below are biology
(marker coverage and marker PPV), not caller bugs.**

## 5. Results — the v1 → v2 → v3 arc

Each measurement committed a falsifiable projection in its pre-registration, then reported
against it. The confusion counts move because abstentions differ; the bar does not move.

### 5.1 v1 — narrow S639 panel: FAIL on under-detection

Scored 87 (40 R / 47 S); 11 abstained. Confusion **TP=36 FN=4 FP=1 TN=46**.

| Metric | Value | Wilson 95% | Bar | Verdict |
|---|---|---|---|---|
| **VME** (missed R) | 4/40 = 10.0% | [4.0%, 23.1%] | ≤3% pt **and** ≤15% upper | ✗ |
| ME (false R) | 1/47 = 2.1% | [0.4%, 11.1%] | ≤5% | ✓ |
| Abstention | 11/98 = 11.2% | [6.4%, 19.0%] | ≤30% | ✓ |

**FAIL** on VME. Characterisation (`work/RESULTS_diagnostic_coverage_ceiling.md`) showed the
cause is tractable, not an off-target wall: **all 4 misses are FKS1 HS3** — B19617 (W691L),
B19618 (W691L), B22769 (W691C), B21978 (M690I, unresolved even by the source paper) — a
hotspot the v1 caller had **no window for** (it read only HS1 635–643 and HS2 1350–1358). The
11 abstentions were mostly HS1 substitutions the caller *already reads but did not tag*
(7 × D642Y, F635C/Y, plus HS2 R1354S). The projected fix — add an HS3 window — put VME at
1/40 = 2.5% [0.4, 12.9] (clears bar), pending a caller change and a fresh measurement.

### 5.2 v2 — HS3 window + widened HS1/HS2 (#154): the failure mode flips to over-calling

Prereg sha `014094c6…`, with a pre-committed VME projection of 1/40 = 2.5%. A fresh GCP
reads→caller pass re-harvested the fixture (`fks1_accuracy_98_v2.tsv`, sha `a7c5056f…`).
Scored 97 (45 R / 52 S); 1 abstained. Confusion **TP=45 FN=0 FP=7 TN=45**.

| Metric | Value | Wilson 95% | Bar | Verdict | v1 → v2 |
|---|---|---|---|---|---|
| **VME** (missed R) | 0/45 = 0.0% | [0.0%, 7.9%] | ≤3% pt / ≤15% upper | ✓ | 10.0% → **0.0%** |
| **ME** (false R) | 7/52 = 13.5% | [6.7%, 25.3%] | ≤5% | ✗ | 2.1% → **13.5%** |
| Abstention | 1/98 = 1.0% | [0.2%, 5.6%] | ≤30% | ✓ | 11.2% → 1.0% |

The HS3 extension did exactly what it projected — VME collapsed to 0% — but the panel now
**over-calls**: **FAIL** on ME. The 7 major errors were genotype-correct but
phenotype-discordant: 5 of 7 were **D642Y** carriers, all phenotypically S; across the whole
benchmark D642Y carriers split **2 R / 5 S** (~29% PPV). There is no genotype rule that wins —
tag D642Y and you get 5 major errors; drop it and its 2 R-only carriers (B20717, B21114)
become very-major misses. One M690I (B20326, S) and one S639Y (B20673, S — irreducible noise
at an otherwise-reliable core position) complete the seven. The binding constraint had moved
from **coverage** to **marker PPV**.

### 5.3 v3 — a PPV tier over the same panel: PASS

The honest response to a PPV ceiling is not more markers but a *narrower claim*, pre-registered
fresh (`work/PREREGISTRATION_diagnostic_accuracy_v3.md`, sha `f3f9686e…`). v3 adds **only** a
PPV tier over the v2 panel — no new window, position, or letter:

- **Core (high-PPV) tier — DETECT:** F635{C,Y}, S639{F,P,Y} (HS1); R1354{S} (HS2); W691{L} (HS3).
- **Low-PPV tier — ABSTAIN as `UNCHARACTERIZED_VARIANT`:** D642Y, M690I.

`openafr/fks1_caller.py` exposes `panel_tier`/`core_panel_hits`; `openafr/verdict.py`
`_classify` drives `RESISTANCE_MARKER_DETECTED` off the core hits, a non-core panel hit
falling through to abstention. The re-tiered fixture was **derived offline** — no SRA re-read,
no GCP: `scripts/harvest_diagnostic_accuracy_v3.py` re-ran the production
`echinocandin_verdict()` over the caller output already frozen in the v2 fixture
(`fks1_accuracy_98_v3.tsv`, sha `02fc37e1…`), asserting the invariant that a verdict may change
*only* `RESISTANCE_MARKER_DETECTED` → `UNCHARACTERIZED_VARIANT`, and only for a low-PPV-only call.

Scored 88 (42 R / 46 S); 10 abstained. Confusion **TP=42 FN=0 FP=1 TN=45**.

| Metric | Value | Wilson 95% | Bar | Verdict | v2 → v3 |
|---|---|---|---|---|---|
| **VME** (missed R) | 0/42 = 0.0% | [0.0%, 8.4%] | ≤3% pt / ≤15% upper | ✓ | 0.0% → 0.0% |
| **ME** (false R) | 1/46 = 2.2% | [0.4%, 11.3%] | ≤5% | ✓ | 13.5% → **2.2%** |
| Sensitivity (1 − VME) | 42/42 = 100.0% | [91.6%, 100%] | — | — | 100% → 100% |
| Specificity (1 − ME) | 45/46 = 97.8% | [88.7%, 99.6%] | — | — | 86.5% → **97.8%** |
| Categorical agreement | 87/88 = 98.9% | [93.8%, 99.8%] | — | — | 92.8% → **98.9%** |
| Abstention | 10/98 = 10.2% | [5.6%, 17.8%] | ≤30% | ✓ | 1.0% → 10.2% |

**PASS** — every gated metric clears its bar, matching the frozen projection (VME 0/42, ME
1/46, abstention 10/98) exactly. The 9 low-PPV-only carriers (7 D642Y, 2 M690I) move from the
confusion matrix into abstention: the 3 resistant among them (D642Y-only B20717/B21114; M690I
B21978) become honest `UNCHARACTERIZED_VARIANT` — not very-major misses — and the 6 susceptible
among them abstain instead of being over-called R, which is what fixes the v2 major-error
failure. Only **B20673** (S639Y, phenotype S) remains a single major error — irreducible
genotype↔phenotype noise at an otherwise-reliable core position, well within the ≤5% bar.

## 6. Discussion

**Scope the claim, don't widen the panel.** The arc is the argument. Widening coverage (v1→v2)
was necessary but exposed a deeper limit — some genotypes simply do not predict phenotype —
that no amount of further coverage can fix. The move that turned a FAIL into a defensible GO
was *narrowing the claim to what the data support*: detect high-PPV markers, abstain honestly
elsewhere. The 10.2% abstention is the visible cost, and reporting it openly is the point — an
engine that abstains on 10% and is right on the rest is more useful, and far safer, than one
that classifies everything and is wrong on 1 susceptible isolate in 7.

**D642Y as the worked example of an independently-falsified marker.** The v3 tiering obeys a
load-bearing anti-circularity rule: a position's tier must rest on evidence *independent of the
scored isolates' own phenotype*. D642Y is demoted because the benchmark authors' own expert
`expected_panel` scores every D642Y carrier as `-` (wild-type) — an independent expert
judgement — corroborated by the CDC EID 25-0760 surveillance collection not treating it as
resistance-defining; M690I because it is phenotype-unresolved in the source literature. The
benchmark then *independently* rejects D642Y on phenotype (5/7 carriers S). Symmetrically,
W691C and W691F were **excluded from the core tier** despite being HS3: W691C's only report is
the benchmark itself (tagging it would be circular; it is emitted verbatim but untagged, so its
carrier B22769 lands in abstention, not detection), and W691F is reported in no source. Only
W691L, additionally CRISPR-confirmed (Jacobs et al., AAC 2023, aac.00423-23 / PMC10269051),
enters the core tier. Tiering by independent literature status — not by fitting this cohort's
discordant rows — is what keeps the PASS honest.

**What this licenses, and no more.** A measured RUO *detection* accuracy claim over validated
high-PPV FKS1 markers (F635/S639/R1354S/W691L core), abstaining on low-PPV positions. It is
**not** a general R/S classifier: on this cohort no FKS1-genotype panel clears both error bars
as a full classifier, because D642Y/M690I carriers are phenotypically split. v3 wins by scoping
the claim, not by out-classifying the biology.

## 7. Limitations (quantified and carried, not retired early)

Same discipline as the geometry preprint's §6 and the go/no-go register
(`docs/DIAGNOSTIC_GO_NO_GO.md`): each limit is stated with the specific evidence that would
lift it.

1. **The azole/ERG11 accuracy arm is out of scope for this publication.** Only 4 Lockhart clade
   strains carry paired ERG11 genotype + azole phenotype in hand — a sanity control, not a
   measurable accuracy — so no azole accuracy is claimed and this preprint's PASS is
   echinocandin-only. This is a declared scope boundary on a named data gap, not an open to-do
   (`docs/AZOLE_ARM_SCOPE_DECISION.md`). *Lifts it:* a larger paired *C. auris* ERG11
   genotype+phenotype collection (the go/no-go's stated single most valuable next data acquisition).
2. **No calibrated probability / no MIC — out of scope for this publication.** The engine is
   categorical by construction: #132 found no public *Candida* panel that clears the calibration
   bar, and the calibrated-probability track #136 is blocked on a non-public pooled *C. albicans*
   ERG11 panel. The calibration machinery (`openafr/calibration.py`) and its pre-registration are
   complete and frozen, waiting only on that data — a scoped-out boundary with a named restart
   trigger, not a defect (`docs/AZOLE_ARM_SCOPE_DECISION.md`). *Lifts it:* that option-C panel +
   a held-out reliability PASS (Rung B).
3. **`NO_KNOWN_MARKER` is not "susceptible."** The panel covers named ERG11/FKS1 markers only —
   not efflux (TAC1/MRR1/CDR1), ERG3, or promoter/TR mechanisms. The measured VME (0% under the
   scoped claim) is honest *within* that coverage; off-panel resistance is out of scope by
   declaration. *Lifts it:* broader-mechanism callers (separate work, not required for the
   current claim).
4. **Retrospective, single cohort, single organism.** The measurement is against one published
   retrospective benchmark (PMC12323592), *C. auris* only; the caller windows FKS1 hot-spots on
   real reads but was validated on one truth set. *Lifts it:* a pre-registered *prospective*
   study on consecutively collected isolates (Rung D), and cross-cohort replication.
5. **FKS1 is detection-only, no structural so-what.** Echinocandins coordinate no metal, so the
   CYP51 geometry moat does not transfer; the FKS1 verdict carries no `structural` field. This
   is a declared scope, not an omission.
6. **Regulatory path untouched.** LDT/IVD quality-system, reproducibility, and lot-to-lot
   evidence (Rung E) are out of scope for the science track; named so they are not discovered
   late.

Nothing here has been tested prospectively against consecutively collected clinical isolates,
and every verdict is RUO. The overall clinical-path verdict is **RUO, Rung A now cleared for
echinocandin under a narrowed high-PPV detection claim; NO-GO at every rung above.**

## 8. Pre-registration ledger

Every measurement of the diagnostics track, in order, with its frozen hash and committed
outcome:

| run | question | pre-reg sha256 | fixture sha256 | outcome |
|---|---|---|---|---|
| FKS1 concordance | caller token accuracy vs genotype truth (98 isolates) | `da8e3827…` | `68eef59f…` (benchmark) | **PASS** (sens/spec/identity 100%) |
| #137 v1 | verdict accuracy, narrow S639 panel | `70d5d2ae…` | `1aeffea5…` | **FAIL** (VME 10.0% [4.0, 23.1]) |
| #137 v2 | + HS3 window, widened HS1/HS2 (#154) | `014094c6…` | `a7c5056f…` | **FAIL** (ME 13.5% [6.7, 25.3]) |
| #137 v3 | + PPV tier (core detect, low-PPV abstain) | `f3f9686e…` | `02fc37e1…` | **PASS** (VME 0%, ME 2.2%) |

Each pre-registration fixed its metric, bar, code change, anti-circularity rule, re-derivation
protocol, and a falsifiable projection **before** grading; the grader re-checks the prereg,
caller-reference, and fixture hashes at run time and refuses to certify if any moved. v1→v2 each
required a fresh GCP reads→caller pass (2026-09-20 concordance harvest; 2026-09-24 v2
re-measure); v3 was derived offline from the frozen v2 caller output (no SRA re-read), a change
allowed to move a verdict only from detection to abstention.

## 9. Data and code availability

All code, all pre-registrations with their hashes, the transcribed hash-pinned benchmark
fixture, the harvested per-version accuracy fixtures, and every `RESULTS_*.md` are in this
repository. Key modules: the callers (`openafr/recaller.py`, `openafr/fks1_caller.py`), the
single interpretation entrypoint (`openafr/interpret.py`), the verdict schema
(`openafr/verdict.py`), the metric engines (`openafr/accuracy.py`, `openafr/concordance.py`),
the report layer (`openafr/report.py`), the graders
(`scripts/validate_fks1_concordance.py`, `scripts/validate_diagnostic_accuracy.py`), the
offline v3 harvester (`scripts/harvest_diagnostic_accuracy_v3.py`), and the coverage-ceiling
characterisation (`scripts/characterize_coverage_ceiling.py`). Truth set: PMC12323592 (Misas et
al., Microbiol Spectr 2025), transcribed by hand from the paywalled Table 1. The GCP re-callers
run on a lean linux/amd64 host; SRA/tool details are in `work/RUNBOOK_fks1_run.md` and
`work/RUNBOOK_fks1_concordance_and_accuracy.md`. Licensed PolyForm Noncommercial 1.0.0.

---

### Note on framing

The ingredients here are not new. FKS1 hot-spot genotyping and CLSI error-rate reporting are
standard. What this repository contributes is a pre-registered, hash-frozen execution against a
published external benchmark that reports its two *failures* — under-detection, then
over-calling — as first-class results, and then earns its PASS by *narrowing the claim to what
the data support* rather than by widening a panel until the numbers move. The engine's
differentiator is not a better number; it is that it says exactly what it did and did not check,
abstains rather than guess, and refuses to manufacture the susceptibility call, the probability,
and the coverage it does not have. The honest boundary is the product.
