# Scope decision: the azole/ERG11 arm and calibrated probability are out of scope for the current publication (issue #160 / A2)

**Decided 2026-09-30.** This record closes the open question carried by issue **#160 (A2)** —
*"power the azole diagnostics arm, or document it out-of-scope with a stated reason"* — and
resolves, for the manuscript, the #136 calibration block. It is a scope decision, not a
retraction: the science is paused on a named data gap, not abandoned. It applies the same
discipline as `docs/DIAGNOSTIC_GO_NO_GO.md` — state where the claim stops, and name exactly what
would move it.

---

## Decision

For the current publication (the JCM manuscript `work/MANUSCRIPT_fks1_diagnostics_JCM.md` and the
companion diagnostics preprint `work/PREPRINT_diagnostics_engine.md`), the following are **declared
out of scope**, firmly and in advance — not deferred to a later rung of this paper:

1. **The azole/ERG11 resistance-accuracy arm.** No measured accuracy is claimed for the azole
   drug class. The published accuracy claim is **echinocandin / FKS1 detection only**.
2. **The calibrated-probability track (#136).** The engine emits categorical verdicts with no
   probability field. No calibrated resistance probability or predicted MIC is claimed.

These join the already-declared out-of-scope boundaries (off-panel mechanisms, structural
interpretation for FKS1, the regulatory path). The point of declaring them here is that they are
**named now, so they are not discovered late** at review.

## Why — this is a data-availability decision, not a code or effort decision

Powering either arm this cycle is infeasible because the required data does not publicly exist.
The engine and its pre-registration are finished; only the input is missing.

- **The calibration machinery is built and frozen, waiting on data.** `openafr/calibration.py`
  is complete and unit-tested (`tests/test_calibration.py`); its pre-registration
  `work/PREREGISTRATION_calibration.md` is **FROZEN and SHA-pinned**
  (`work/PREREG_calibration.sha256`) and states on its face that the run it governs is **BLOCKED
  on data, not runnable today.** Nothing about the method is in doubt — there is simply no panel
  to run it on.
- **The #132 audit found no qualifying public panel.** `work/RESULTS_diagnostics_panel_spike.md`
  (against frozen bar `1211e973…`) searched for an obtainable public *Candida* collection meeting
  the calibration bar — **N≥150, ≥40 per class, ≥5 variants with ≥5 carriers each, single paired
  genotype+MIC collection** — and found none. Verdict: **FAIL → fallback to concordance-only.**
- **The only adequate ERG11 data is the wrong organism.** The sole variant-diverse, susceptible-
  inclusive public ERG11 data is *C. albicans* (option C — a *pooled* panel that does not yet
  exist assembled). Our re-caller is *C. auris*-specific, so it cannot validate *C. albicans*
  data without new caller/structure work.
- **The azole arm has only a sanity control, not a powered sample.** Just **4 Lockhart clade
  strains** carry paired ERG11 genotype + azole phenotype in hand — enough to sanity-check, not
  to measure accuracy.
- **The novel-variant regime has no continuous score by design.** `openafr/structural.py`
  deliberately emits a *direction*, never a magnitude (`openafr/calibration.py:18-22`), so there
  is no continuous score to fit a probability curve to; fabricating one would claim a precision
  the module refuses to assert.

Separately, azole resistance in *C. auris* is already near-saturated at baseline (azole-marker
event frequency 160/199 = 80.4% [74.3, 85.3]; `work/RESULTS_prevalence.md`), so azole *emergence*
carries little early-warning signal — the dynamic, surveillance-useful axis is echinocandin/FKS1
(10/443 = 2.3%). The scientific center of gravity for this publication is where the measurement
landed.

## Named unlock conditions (what flips each back in scope)

These are recorded so a future rung cannot be waved through on hindsight, and so the paused work
has a concrete restart trigger.

- **Azole/ERG11 arm ←** a **larger paired *C. auris* ERG11 genotype+phenotype collection**
  (the 4-strain set only supports a sanity control). This remains the project's **single most
  valuable next data acquisition** (`docs/DIAGNOSTIC_GO_NO_GO.md` §6).
- **Calibrated probability (#136, Rung B) ←** the assembled **option-C pooled *C. albicans* ERG11
  panel**, plus the three pre-committed costs recorded in the frozen calibration prereg: a
  *separately frozen* prereg amendment for cross-study MIC-method heterogeneity, new *C. albicans*
  caller/structure work, and verified manual transcription — then a held-out reliability PASS
  (ECE ≤0.10, Brier ≤0.20, in-large gap ≤0.10, ≥15/class).

## What remains in scope (unchanged by this decision)

The echinocandin / FKS1 **detection-accuracy** claim — measured and **passing** under the narrowed
high-PPV claim (v3: VME 0/42 = 0.0%, ME 1/46 = 2.2%, abstention 10.2%;
`work/RESULTS_diagnostic_accuracy.md`). The engine remains RUO throughout; no clinical
susceptibility claim is made for any drug class.

## Disposition

- **#160 (A2):** resolved — by formal scope-out, with the stated reason and unlock conditions
  above.
- **#136 (calibration):** the manuscript block is resolved. The track is **paused on a named data
  gap, not killed**; the unlock route is recorded here and in `docs/DIAGNOSTIC_GO_NO_GO.md`
  (Rung B).
