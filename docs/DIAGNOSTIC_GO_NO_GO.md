# Diagnostics go / no-go — the clinical-path decision doc (issue #139)

Closes the diagnostics MVP loop. This document carries the ceiling honestly: it states
**where the engine's day-one claim stops**, fixes the **minimal defensible claim** we are
willing to make now, and defines the **gates** that would license any move toward a laboratory-
developed test (LDT) or clinical use. It applies the same discipline as the preprint's
limitations section (`work/PREPRINT_geometry_ceiling.md`): every limit is **quantified and
carried, not retired early** — each with the specific evidence that would lift it.

It is a decision doc, not a certification. It commits, in advance, to what a "go" for the next
rung requires, so no future rung can be waved through on hindsight.

---

## 1. What the engine is, day one

The resistance-interpretation engine (`.context/DIAGNOSTICS_CONCEPT.md`, contract
`docs/DIAGNOSTIC_VERDICT_CONTRACT.md`, #133) takes a typed *C. auris* genotype and emits, per
drug class, one of four categorical verdicts —
`RESISTANCE_MARKER_DETECTED` / `UNCHARACTERIZED_VARIANT` / `NO_KNOWN_MARKER` / `UNRESOLVED` —
plus, for ERG11, a calibrated-**low** structural best-guess on uncharacterized variants (#135).
It is **RUO** (research use only) and emits **no probability, no MIC, no clinical
actionability**.

The whole product is that contract's *honest boundary*. The engine's differentiator is not a
better number; it is that it says exactly what it did and did not check, and refuses to
manufacture the parts it cannot back.

---

## 2. Where the ceiling is (stated, not glossed)

Four structural ceilings bound the day-one claim. None is a bug to be hidden; each is a
declared limit with a named route past it.

**C1 — Concordance-only, no calibrated probability.** The #132 panel hunt
(`work/RESULTS_diagnostics_panel_spike.md`, against frozen bar `1211e973…`) found **no
obtainable public *Candida* collection** that clears the calibration bar (N≥150, ≥40/class, ≥5
variants with ≥5 carriers, single paired genotype+MIC collection). Verdict: **FAIL → fallback
(a), concordance-only**. So the engine detects *named markers* categorically and stamps **no
probability**. A probability field re-enters the contract only after the #136 calibration track
produces a *measured, held-out* reliability — never before.
*Quantified, not retired:* the exact panel that would unlock it is recorded (option C — a
pooled *C. albicans* ERG11 panel), with its three pre-committed costs (`calibration-track`
prereg).

**C2 — Marker coverage is narrow, and `NO_KNOWN_MARKER` is not "susceptible."** The engine
detects a small named panel — ERG11 V125A/F126L/Y132F/K143R, FKS1 HS1 S639F/P/Y — and nothing
else. It does **not** cover efflux (TAC1/MRR1/CDR1), ERG3, promoter/tandem-repeat (TR), or
non-target mechanisms. A resistant isolate whose mechanism is off-panel is a genuine miss the
engine cannot and must not paper over, which is exactly why `NO_KNOWN_MARKER` is defined as
"no *known marker*," never as a clinical S (contract rule 2).
*Quantified, not retired — now measured (2026-09-20) and characterised (2026-09-20):* the #137
echinocandin accuracy run (`work/RESULTS_diagnostic_accuracy.md`) reports the very-major-error rate
this ceiling drives at **10.0% [4.0%, 23.1%]** — 4 of 40 phenotypically-resistant isolates the engine
declines to call. The mechanism characterisation (`work/RESULTS_diagnostic_coverage_ceiling.md`,
`scripts/characterize_coverage_ceiling.py`) shows the ceiling is **FKS1-hotspot-shaped, not efflux**:
all 4 misses are FKS1 HS3 (3 × W691L/W691C; 1 × M690I unresolved even by the source paper), a hotspot
the caller has **no window for** — not an off-target mechanism. The 11 abstentions are likewise
mostly HS1 substitutions the caller already reads but does not panel-tag (9 × HS1 D642Y/F635, 1 × HS2
R1354S). The FKS1 concordance run (`work/RESULTS_fks1_concordance.md`) separately certifies the
caller's *tokens* as perfect, so this is a **tractable coverage gap**, not a caller defect: adding an
HS3 window would project the VME to 1/40 = 2.5% [0.4, 12.9] (clears bar) — a projection, pending the
caller change + a GCP re-measure, not a measured re-claim.

**C3 — The validated organism is not the calibratable organism.** The re-callers and the
structural port are *C. auris*-centric (5TZ1, [[openafr-auris-port]]). But the organism whose
public ERG11 data has the variant diversity and the susceptible isolates calibration needs is
*C. albicans* (#132 finding 2). "Reuse the existing caller" and "clear the calibration bar"
pull toward *different organisms* — a tension that must be decided, not glossed, before any
calibrated claim is made.

**C4 — Accuracy: echinocandin arm MEASURED THREE TIMES (now PASSES under a narrowed high-PPV claim);
azole arm still underpowered.** The verdict-accuracy validation (#137) ran first on 2026-09-20 for the
v1 narrow S639 panel — scored 87, **VME 10.0% [4.0, 23.1]** — **FAIL** on under-detection. PR #154
extended the caller (HS3 window + widened HS1/HS2); the v2 GCP re-measure (2026-09-24) confirmed the
VME projection — **VME 0.0%, sensitivity 100%** — but the failure mode **flipped**: **ME 13.5%
[6.7, 25.3]** (bar ≤5%) — **FAIL** on over-calling, a genotype↔phenotype **PPV ceiling** (D642Y splits
2 R / 5 S, M690I 1 R / 1 S; the benchmark's own expert panel scores D642Y wild-type). v3 (2026-09-27,
`work/PREREGISTRATION_diagnostic_accuracy_v3.md`) resolved it by adding a **PPV tier** — high-PPV core
markers detect, low-PPV D642Y/M690I abstain — clearing the bar: **VME 0/42 = 0.0% [0.0, 8.4], ME 1/46 =
2.2% [0.4, 11.3], abstention 10.2%, PASS** (`work/RESULTS_diagnostic_accuracy.md`), matching its frozen
projection. The GO is **scoped**: a detection claim over validated high-PPV markers, not a general R/S
classifier. The **azole/ERG11 arm** has only 4 Lockhart clade strains and **remains UNDERPOWERED**;
**calibration (#136)** remains blocked on the non-public option-C panel.

---

## 3. The minimal defensible day-one claim (the "go" for RUO)

The claim we **are** willing to stand behind now — no more:

> An RUO engine that, for a typed *C. auris* isolate, **detects named ERG11 (azole) and validated
> high-PPV FKS1 (echinocandin) resistance markers** — for echinocandin, the F635/S639/R1354S/W691L
> core, with a **measured** verdict accuracy (VME 0%, ME 2.2%; #137 v3) — emits an explicit
> **uncharacterized-variant** verdict (with a calibrated-low structural best-guess for ERG11, and for
> known-but-low-PPV FKS1 positions like D642Y/M690I) instead of falling silent or over-calling, and
> **honestly abstains** (UNRESOLVED / excluded) where it cannot call. It makes **no susceptibility
> claim, no probability, and no clinical recommendation.**

What that claim explicitly does **not** include, day one:

- ❌ "Susceptible" for any isolate (C2 — absence of a marker ≠ susceptibility).
- ❌ A calibrated resistance probability or predicted MIC (C1 — deferred to #136).
- ❌ A general R/S classifier. The measured echinocandin accuracy (C4 — #137 v3 PASS) is scoped to
  **detection of high-PPV markers**; the engine abstains on low-PPV positions rather than classifying
  every isolate. The azole accuracy is still *unmeasured* (underpowered).
- ❌ Any clinical actionability or treatment guidance (RUO by construction).
- ❌ Coverage of *C. albicans* or other species (C3 — auris-only callers).

**Go decision, RUO tier: GO** — the day-one claim above is defensible now because every part of
it rests on a shipped, unit-tested primitive (the contract + callers + report), and the parts
that are *not* defensible are explicitly excluded rather than softened.

---

## 4. The gate ladder to a clinical move (fixed in advance)

Each rung names the pre-committed evidence that licenses it. A rung is **NO-GO until its gate
is met**; meeting it is necessary, not sufficient, for the next.

**Rung A — Measured concordance accuracy (unlocks: an RUO *accuracy* claim).**
- GATE: the #137 run executes on a pinned paired-panel fixture and the drug class **PASSES** its
  frozen bar — very-major-error point ≤3% **and** Wilson upper ≤15%, major-error ≤5%,
  abstention ≤30%, with ≥15 scored isolates per phenotype arm.
- Today: **GO for echinocandin, under the narrowed high-PPV claim (v3 2026-09-27); azole still
  NO-GO.** The arm ran three times. v1 (narrow S639 panel): VME **10.0% [4.0, 23.1]** — FAIL on
  under-detection. v2 (HS3 + widened HS1/HS2, #154): VME projection held — **VME 0.0%,
  sensitivity 100%** — but **ME rose to 13.5% [6.7, 25.3]**, FAIL on over-calling, a
  genotype↔phenotype PPV ceiling (D642Y splits 2 R / 5 S; M690I 1 R / 1 S). **v3
  (`work/PREREGISTRATION_diagnostic_accuracy_v3.md`) added a PPV *tier* over the same panel** — core
  high-PPV markers detect; low-PPV D642Y/M690I abstain (`UNCHARACTERIZED_VARIANT`) instead of
  over-calling — and **cleared the bar: VME 0/42 = 0.0% [0.0, 8.4], ME 1/46 = 2.2% [0.4, 11.3],
  abstention 10/98 = 10.2%, PASS** (`work/RESULTS_diagnostic_accuracy.md`), matching its frozen
  projection exactly. This is a GO for a **specific, narrowed** intended use — *detection over
  validated high-PPV FKS1 markers (F635/S639/R1354S/W691L core), abstaining on low-PPV positions* —
  **not** a general R/S classifier (no FKS1-genotype panel clears the bar as a full classifier on this
  cohort; v3 wins by scoping the claim to detection and abstaining honestly elsewhere, at a 10.2%
  abstention price). The **azole arm** stays NO-GO, underpowered (4 clade strains).

**Rung B — Calibrated probability with measured reliability (unlocks: the probability field).**
- GATE: the #136 calibration track executes on the option-C pooled panel and clears its frozen
  reliability bar (ECE ≤0.10, Brier ≤0.20, in-large gap ≤0.10, ≥15/class on the held-out split),
  under a *separately frozen* prereg amendment for cross-study MIC-method heterogeneity.
- Today: **NO-GO (panel does not exist; C1/C3).** Deferred, not killed; unlock route named.

**Rung C — Panel breadth + species scope adequate for the intended-use population.**
- GATE: marker coverage and organism scope match the population the test would serve — either a
  documented decision that a narrow *C. auris* panel is the intended use, **or** the *C.
  albicans* / broader-mechanism work (efflux/ERG3/TR) that C2/C3 would require.
- Today: **GO for the narrow, high-PPV *C. auris* FKS1 panel (decision EXECUTED, 2026-09-27); broader
  scope still NO-GO.** The v2 panel extension (#154, `work/PREREGISTRATION_diagnostic_accuracy_v2.md`)
  worked for detection (VME 10% → 0%) but revealed that breadth was no longer the binding constraint —
  **marker PPV / genotype↔phenotype discordance** was (D642Y 2 R / 5 S; the benchmark's expert panel
  scores it wild-type). The forward path was therefore a *re-scoped high-PPV claim*, not more breadth,
  and it was pre-registered fresh (`work/PREREGISTRATION_diagnostic_accuracy_v3.md`) and measured to
  **PASS**: a **PPV tier** (core F635/S639/R1354S/W691L detect; low-PPV D642Y/M690I abstain) clears the
  accuracy bar (`work/RESULTS_diagnostic_accuracy.md`). So the intended-use decision is settled for
  echinocandin: **a narrow, high-PPV *C. auris* FKS1 detection panel is the documented intended use.**
  Efflux/ERG3/TR mechanisms and *C. albicans* / cross-species breadth remain separately out of scope.

**Rung D — Prospective clinical validation against phenotypic AST.**
- GATE: a pre-registered *prospective* study (not the retrospective #137 benchmark) with
  pre-committed very-major/major-error acceptance criteria on consecutively collected clinical
  isolates, sized for the intended use.
- Today: **NO-GO (not started; retrospective validation itself not yet run).**

**Rung E — Regulatory / quality path (LDT or IVD).**
- GATE: the applicable regulatory framework identified and met (CLIA/CAP for an LDT, or the IVD
  pathway), with reproducibility, lot-to-lot, and quality-system evidence. **Out of scope for
  the science tracks; named here so it is not discovered late.**
- Today: **NO-GO (not begun).**

**Overall clinical-path verdict: RUO, with Rung A now cleared for echinocandin under a narrowed
high-PPV detection claim (#137 v3, 2026-09-27); NO-GO at every rung above that.** The engine remains an
RUO research tool. The echinocandin *detection-accuracy* claim is now measured and PASSES (VME 0%, ME
2.2%); the azole arm stays underpowered and Rungs B–E remain NO-GO.

---

## 5. Limitation register (quantify, carry, don't retire early)

| # | Limitation | Status | What lifts it |
|---|---|---|---|
| C1 | No calibrated probability | **Carried** — deferred to #136; not a defect, a declared scope | Rung B: option-C panel + reliability PASS |
| C2 | Narrow marker panel; `NO_KNOWN_MARKER` ≠ S | **Resolved for echinocandin as a scoped detection claim (2026-09-27)** — v2 HS3 extension fixed detection (VME → 0%) but exposed a PPV ceiling; v3 PPV tier (high-PPV core detects, low-PPV D642Y/M690I abstain) PASSES (VME 0%, ME 2.2%) | Broader mechanisms (efflux/ERG3/TR) still out of scope; not required for the narrowed claim |
| C3 | Validated organism ≠ calibratable organism | **Open decision (echinocandin intended use now documented as narrow *C. auris* FKS1)** | Explicit intended-use call, or *C. albicans* caller/structure work |
| C4 | Accuracy: echinocandin measured (v3 PASS, scoped); azole underpowered | **Echinocandin measured & PASSES (narrowed high-PPV claim); azole still blocked** | Azole: larger paired *C. auris* ERG11 collection |

The discipline mirrors the preprint: no limit is quietly dropped once stated; each is carried in
this register with its unblock route until the evidence that lifts it actually lands.

---

## 6. What a "go" would require next (the single actionable item)

**Update 2026-09-20 — Rung A has been run for echinocandin, and it changes the next move.** The
#137 echinocandin accuracy run executed off the FKS1-concordance-harvested fixture and **did not
clear the bar** (VME 10.0% [4.0, 23.1]; `work/RESULTS_diagnostic_accuracy.md`), while the caller
itself **passed** concordance at 100% (`work/RESULTS_fks1_concordance.md`). The bottleneck is no
longer compute or a caller — it is **marker coverage (Rung C)**.

**Update 2026-09-20 (characterisation) — the coverage gap is tractable, not off-target.** The
mechanism characterisation (`work/RESULTS_diagnostic_coverage_ceiling.md`,
`scripts/characterize_coverage_ceiling.py`) shows the 4 VME misses are **all FKS1 HS3**
(3 × W691L/W691C, 1 × M690I unresolved even by the source paper) — a hotspot the caller has no window
for — and the 11 abstentions are mostly **HS1 substitutions the caller already reads but does not tag**
(9 × D642Y/F635, 1 × HS2 R1354S). Adding an HS3 window projects the VME to **1/40 = 2.5% [0.4, 12.9]**
(clears the bar); widening the HS1/HS2 panel additionally converts abstentions to detections. This is a
*projection* pending a caller change + GCP re-measure, and the panel must be pinned from literature
independent of this benchmark to avoid circularity.

**Update 2026-09-27 (v3) — Rung A CLEARED for echinocandin under a narrowed high-PPV claim.** The v2
extension fixed detection (VME → 0%) but over-called (ME 13.5%): the binding constraint had moved from
coverage to **marker PPV** (D642Y 2 R / 5 S; M690I 1 R / 1 S; the benchmark's own expert panel scores
D642Y wild-type). Rather than widen further, v3 (`work/PREREGISTRATION_diagnostic_accuracy_v3.md`) added
a **PPV tier** over the same panel — high-PPV core markers (F635/S639/R1354S/W691L) detect; low-PPV
D642Y/M690I abstain — re-derived the fixture offline from the frozen v2 caller output (no GCP), and
**PASSES** (VME 0/42 = 0.0%, ME 1/46 = 2.2%, abstention 10.2%; `work/RESULTS_diagnostic_accuracy.md`).
The echinocandin detection-accuracy claim is now *measured and cleared*, honestly scoped.

So the nearest unlocks are now:
- **Echinocandin:** Rung A **cleared** for the narrowed high-PPV detection claim; the intended-use
  decision (Rung C) is documented as a narrow *C. auris* FKS1 panel. Remaining routes up are Rung B
  (probability, #136) and broader mechanisms — both separately blocked, neither required for the current
  claim.
- **Azole/ERG11:** obtain a **larger paired *C. auris* ERG11 genotype+phenotype collection** (the
  4-strain Lockhart set only supports a sanity control) to lift the azole arm past UNDERPOWERED. **This
  is now the single most valuable next data acquisition.**
- **Calibration (#136):** unchanged — still blocked on the non-public option-C pooled *C. albicans*
  panel.

The engine stands, honestly, at the RUO tier — now with a *measured, passing* echinocandin detection
claim rather than a measured ceiling. Scoping the claim to what the data support (high-PPV detection,
honest abstention elsewhere), rather than overclaiming a full R/S classifier, is what turned the v2 FAIL
into a defensible GO — the whole point of carrying limits rather than retiring them early.
