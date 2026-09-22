# Echinocandin coverage-ceiling characterisation — the #137 VME is FKS1-hotspot-shaped, not efflux

**What this resolves.** The #137 echinocandin verdict-accuracy arm FAILed its Rung-A bar with a
very-major-error (VME) rate of **4/40 = 10.0% [4.0%, 23.1%]**
([`RESULTS_diagnostic_accuracy.md`](RESULTS_diagnostic_accuracy.md)). `docs/DIAGNOSTIC_GO_NO_GO.md`
left the next move as an *open* Rung C decision — "document a narrow-panel intended use, **or** extend
coverage to off-panel mechanisms (efflux / non-HS1 FKS1 hotspots)" — because the *mechanism* behind
each of the 4 misses had not been characterised. This memo characterises them, from data already on
disk, and turns that open decision into an evidence-backed one.

**Reproducibility.** Every mechanism label below is copied verbatim from the benchmark's own
`paper_mut=<sub> (hsN)` annotation (`data/earlywarning/recaller_sanity/fks1_benchmark_100.tsv`, truth
= PMC12323592, Misas et al.), joined to the harvested #137 verdicts
(`fks1_accuracy_98.tsv`) by `scripts/characterize_coverage_ceiling.py`. The script reproduces the
measured confusion (TP=36 FN=4 FP=1 TN=46, VME 4/40 = 10.0%) exactly and uses the same
`openafr.backtest.wilson_interval` CI math as the measured run. Nothing here is hand-authored.

## Headline — the 4 misses are all FKS1 HS3, a hotspot the caller never reads

The caller (`openafr/fks1_caller.py`) reads exactly two windows — **HS1 (635–643)** and
**HS2 (1350–1358)**. There is **no HS3 window**. So a substitution around W691 is invisible: the
caller reads HS1/HS2 wild-type, finds nothing, and returns `NO_KNOWN_MARKER` — scored as a VME.

| VME isolate | mechanism (`paper_mut`) | hotspot | caller reads it? |
|---|---|---|---|
| B19617 | W691L | **HS3** | ✗ no HS3 window |
| B19618 | W691L | **HS3** | ✗ no HS3 window |
| B22769 | W691C | **HS3** | ✗ no HS3 window |
| B21978 | none (footnote: M690I, outside the paper's MycoSNP-nf read region) | HS3-adjacent | not read by anyone |

**None of the 4 misses is efflux, ERG3, TR, or any off-*FKS1* mechanism.** Three are FKS1 HS3
(W691L/W691C); the fourth (B21978) is annotated `paper_mut=none` because even the paper's own
pipeline could not resolve its M690I. The C2 ceiling that drove the 10% VME is therefore **tractable
FKS1-hotspot coverage**, not the fundamental off-target wall the go/no-go doc previously assumed.

## The abstentions tell the same story — under-tagging inside a window already read

The 11 `UNCHARACTERIZED_VARIANT` abstentions (excluded from the 2×2, so they do **not** inflate VME)
are almost all HS1 substitutions the caller **already reads but does not panel-tag**:

- **9 × HS1** — 7 × D642Y, 1 × F635C, 1 × F635Y. The HS1 window (635–643) covers residues 635 and
  642, so the caller *sees and emits* these tokens; it just doesn't tag them (the HS1 panel today is
  only S639{F,P,Y}). Widening the HS1 panel is a **pure data edit** — no new window, no new reads.
- **1 × HS2** — B20592 R1354S. The HS2 window is read but panel-tagged with an empty set (a
  documented gap, `openafr/fks1_caller.py:116`); pinning R1354S tags it.
- **1 × undetermined** — B19896 (`paper_mut=Undetermined`); stays an honest abstention.

## Projection (a projection, not a measured re-claim)

Feeding the recovered confusion through the same Wilson math:

| Change | recovers | projected VME | clears bar? |
|---|---|---|---|
| **Measured today** | — | 4/40 = 10.0% [4.0, 23.1] | ✗ |
| **+ HS3 window** tagging W691L/W691C | B19617, B19618, B22769 | **1/40 = 2.5% [0.4, 12.9]** | ✓ point ≤3%, upper ≤15% |
| **+ HS1/HS2 panel widening** (D642Y, F635C/Y, R1354S) | 10 abstentions → detected | 1/49 ≈ 2.0%; abstention ≈ 1/98 | ✓ (also lowers abstention) |

The residual miss is **B21978** (M690I, unresolved even by the source paper) — the honest floor the
extended panel would still not clear, and the reason the conservative projection is 2.5%, not 0%.

**Two integrity guardrails that make this a projection and not a result:**

1. **No overfitting.** The panel mutant set (D642Y, F635C/Y, R1354S, W691L/C, M690I) must be pinned
   from literature that establishes each as an echinocandin-resistance marker **independently of this
   benchmark**. Tagging positions *because they appear in the test set*, then "recovering" them on the
   same set, would be circular. The projection above is only licensed once each mutant is
   independently pinned.
2. **No paper-number without a re-measure.** The projection is arithmetic on the frozen fixture, not a
   re-run. A measured post-extension VME requires the caller change (HS3 window + widened panel) **and**
   a GCP re-download/re-call under a fresh, frozen pre-registration
   (`work/PREREGISTRATION_diagnostic_accuracy_v2.md`). Until then, Rung A stays a **measured NO-GO**;
   `RESULTS_diagnostic_accuracy.md`'s 10% VME is untouched.

## The Rung C decision this licenses

The evidence favours **extend the panel** over **document a permanently narrow panel**: the miss is
FKS1-hotspot-shaped and mostly recoverable inside windows the caller already reads. The concrete,
pre-registered next step (not done in this pass, because it ripples into the prevalence/emergence
resolution semantics — `_FKS1_PANEL_WINDOWS` is derived from `w.panel` — the `format_call` source
strings, stored snapshots, and six test files, and so deserves its own reviewed PR + GCP re-measure):

1. Add an **HS3 `Window`** (688–698, anchor W691 — verified to translate in the pinned reference)
   with a literature-pinned mutant set {W691L, W691C, W691F, M690I}.
2. **Widen the HS1 panel** to add D642{Y} and F635{C,Y}; **pin the HS2 panel** to add R1354{S}.
3. Freeze `work/PREREGISTRATION_diagnostic_accuracy_v2.md` (bar unchanged), then re-run the FKS1
   concordance+accuracy playbook (`work/RUNBOOK_fks1_concordance_and_accuracy.md`) on GCP to produce
   the *measured* post-extension VME and test this projection.

This keeps the caller certified at what it covers (concordance PASS, 100%/100%/100%,
[[fks1-caller-concordance]]) while naming exactly which three windows close the gap — an honest,
tractable Rung C path rather than an off-target dead end. [[diagnostic-accuracy-track]]

## Implementation update (panel extension landed offline; measured re-claim still pending GCP)

The caller change governed by `work/PREREGISTRATION_diagnostic_accuracy_v2.md` is now implemented in
`openafr/fks1_caller.py`: an HS3 `Window` (688–698, anchor W691) plus widened HS1/HS2 panels. Running
the frozen prereg's **anti-circularity rule** against the literature moved two candidate letters:

- **Kept (each independently pinned, NOT from the PMC12323592 benchmark):** HS1 **F635C/Y, D642Y**;
  HS2 **R1354S**; HS3 **M690I, W691L**. Primary independent source: CDC EID 2026;32(5) article 25-0760
  ("Updated Genomic Epidemiologic Description of *Candida (Candidozyma) auris*, United States"), a US
  surveillance collection distinct from the benchmark's SRA reads; W691L additionally CRISPR-confirmed
  (Jacobs et al., AAC 2023, aac.00423-23 / PMC10269051).
- **Dropped:** **W691F** (reported in no source) and **W691C** (its *only* report is the benchmark
  itself, PMC12323592 isolate SRR26666774 → tagging it would be circular). W691C is still emitted
  *verbatim but untagged*, so its benchmark isolate (B22769) leaves the VME numerator as an honest
  `UNCHARACTERIZED_VARIANT` abstention rather than a `NO_KNOWN_MARKER` miss. The projection is therefore
  unchanged (**VME 1/40 = 2.5%**, residual = B21978 M690I, whose reads under-cover HS3).

This also flips the FKS1 resolution contract: `backtest._FKS1_PANEL_WINDOWS` now spans **all three**
windows, so an isolate is "resolved" only when HS1+HS2+HS3 are each called/wild-type (a snapshot whose
source predates HS3 is judged `partial`, never assumed clear). The offline test suite is green; the
**measured** post-extension VME still requires the GCP re-download/re-call under the v2 prereg, so
`RESULTS_diagnostic_accuracy.md`'s 10% and the go/no-go Rung-A verdict are deliberately **untouched**
until that run lands.
