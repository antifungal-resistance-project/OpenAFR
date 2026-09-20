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
*Quantified, not retired — now measured (2026-09-20):* the #137 echinocandin accuracy run
(`work/RESULTS_diagnostic_accuracy.md`) reports the very-major-error rate this ceiling drives at
**10.0% [4.0%, 23.1%]** — 4 of 40 phenotypically-resistant isolates carry off-panel (non-HS1)
resistance the engine correctly declines to call. The FKS1 concordance run
(`work/RESULTS_fks1_concordance.md`) separately certifies the caller's *tokens* as perfect, so this
is a coverage ceiling, not a caller defect.

**C3 — The validated organism is not the calibratable organism.** The re-callers and the
structural port are *C. auris*-centric (5TZ1, [[openafr-auris-port]]). But the organism whose
public ERG11 data has the variant diversity and the susceptible isolates calibration needs is
*C. albicans* (#132 finding 2). "Reuse the existing caller" and "clear the calibration bar"
pull toward *different organisms* — a tension that must be decided, not glossed, before any
calibrated claim is made.

**C4 — Accuracy: echinocandin arm now MEASURED (and FAILs the clinical bar); azole arm still
underpowered.** The verdict-accuracy validation (#137, `work/PREREGISTRATION_diagnostic_accuracy.md`)
ran on 2026-09-20 for the echinocandin arm, off the fixture harvested by the FKS1 concordance pass
(`work/RESULTS_diagnostic_accuracy.md`): scored 87 isolates (40 R / 47 S), **VME 10.0%
[4.0%, 23.1%]**, ME 2.1%, categorical agreement 94.3%, abstention 11.2% — **FAIL** (VME point and
Wilson upper both above bar; ME and abstention pass). The FAIL is C2's mechanism-coverage ceiling,
not a caller error (the concordance run certified the tokens at 100%/100%/100%,
`work/RESULTS_fks1_concordance.md`). The **azole/ERG11 arm** has only 4 Lockhart clade strains and
**remains UNDERPOWERED**; **calibration (#136)** remains blocked on the non-public option-C panel.
So the echinocandin error rate is now known and does not meet a diagnostic bar; the azole rate is
still unmeasured.

---

## 3. The minimal defensible day-one claim (the "go" for RUO)

The claim we **are** willing to stand behind now — no more:

> An RUO engine that, for a typed *C. auris* isolate, **detects named ERG11 (azole) and FKS1
> HS1 (echinocandin) resistance markers**, emits an explicit **uncharacterized-variant** verdict
> (with a calibrated-low structural best-guess for ERG11) instead of falling silent on novel
> changes, and **honestly abstains** (UNRESOLVED / excluded) where it cannot call. It makes
> **no susceptibility claim, no probability, and no clinical recommendation.**

What that claim explicitly does **not** include, day one:

- ❌ "Susceptible" for any isolate (C2 — absence of a marker ≠ susceptibility).
- ❌ A calibrated resistance probability or predicted MIC (C1 — deferred to #136).
- ❌ A measured accuracy figure (C4 — pending the #137 runs; today's number is *unmeasured*).
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
- Today: **NO-GO — echinocandin arm RAN and did not clear the bar (2026-09-20).** Scored 87
  isolates (40 R / 47 S ≥15 each), VME **10.0% [4.0%, 23.1%]** (bar ≤3% / upper ≤15%) — FAIL; ME
  2.1% and abstention 11.2% pass (`work/RESULTS_diagnostic_accuracy.md`). This is a *measured*
  no-go driven by C2 marker coverage, not a blocked-on-data one — the caller itself passed
  concordance at 100% (`work/RESULTS_fks1_concordance.md`). The **azole arm** stays NO-GO,
  underpowered (4 clade strains). Lifting Rung A for echinocandin now requires **Rung C breadth**
  (off-panel mechanisms), not more compute.

**Rung B — Calibrated probability with measured reliability (unlocks: the probability field).**
- GATE: the #136 calibration track executes on the option-C pooled panel and clears its frozen
  reliability bar (ECE ≤0.10, Brier ≤0.20, in-large gap ≤0.10, ≥15/class on the held-out split),
  under a *separately frozen* prereg amendment for cross-study MIC-method heterogeneity.
- Today: **NO-GO (panel does not exist; C1/C3).** Deferred, not killed; unlock route named.

**Rung C — Panel breadth + species scope adequate for the intended-use population.**
- GATE: marker coverage and organism scope match the population the test would serve — either a
  documented decision that a narrow *C. auris* panel is the intended use, **or** the *C.
  albicans* / broader-mechanism work (efflux/ERG3/TR) that C2/C3 would require.
- Today: **NO-GO (auris-only, HS1/ERG11-hotspot-only).**

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

**Overall clinical-path verdict: NO-GO at every rung above RUO.** The engine is an RUO research
tool. Rung A is the immediate next unlock and is blocked only on obtaining/transcribing a
paired fixture and a compute run — not on any new science.

---

## 5. Limitation register (quantify, carry, don't retire early)

| # | Limitation | Status | What lifts it |
|---|---|---|---|
| C1 | No calibrated probability | **Carried** — deferred to #136; not a defect, a declared scope | Rung B: option-C panel + reliability PASS |
| C2 | Narrow marker panel; `NO_KNOWN_MARKER` ≠ S | **Quantified (2026-09-20)** — echinocandin VME 10.0% [4.0, 23.1] | Rung C for breadth (off-panel/efflux mechanisms) |
| C3 | Validated organism ≠ calibratable organism | **Open decision** | Explicit intended-use call, or *C. albicans* caller/structure work |
| C4 | Accuracy: echinocandin measured (FAIL); azole underpowered | **Echinocandin measured & carried; azole still blocked** | Azole: larger paired *C. auris* ERG11 collection; echinocandin: Rung C breadth |

The discipline mirrors the preprint: no limit is quietly dropped once stated; each is carried in
this register with its unblock route until the evidence that lifts it actually lands.

---

## 6. What a "go" would require next (the single actionable item)

**Update 2026-09-20 — Rung A has been run for echinocandin, and it changes the next move.** The
#137 echinocandin accuracy run executed off the FKS1-concordance-harvested fixture and **did not
clear the bar** (VME 10.0% [4.0, 23.1]; `work/RESULTS_diagnostic_accuracy.md`), while the caller
itself **passed** concordance at 100% (`work/RESULTS_fks1_concordance.md`). The bottleneck is no
longer compute or a caller — it is **marker coverage (Rung C)**: ~1 in 5 phenotypically-resistant
isolates carry off-panel (non-HS1 / efflux) mechanisms this engine does not detect.

So the nearest unlocks are now:
- **Echinocandin:** a **Rung C** decision — either document a narrow-panel intended use that owns the
  10% VME honestly, or extend coverage to the off-panel mechanisms (efflux / non-HS1 FKS1 hotspots).
- **Azole/ERG11:** obtain a **larger paired *C. auris* ERG11 genotype+phenotype collection** (the
  4-strain Lockhart set only supports a sanity control) to lift the azole arm past UNDERPOWERED.
- **Calibration (#136):** unchanged — still blocked on the non-public option-C pooled *C. albicans*
  panel.

The engine stands, honestly, at the RUO tier. The echinocandin accuracy is now a *measured* ceiling,
not an asserted one — which is the whole point of carrying limits rather than retiring them early.
