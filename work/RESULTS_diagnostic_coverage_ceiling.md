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
  only S639{F,P,Y}). F635C/F635Y are concordant and safe to tag; **D642Y is not** (see the
  discordance table below — 5 of its 12 carriers are susceptible), so it stays an abstention.
- **1 × HS2** — B20592 R1354S. The HS2 window is read but panel-tagged with an empty set (a
  documented gap, `openafr/fks1_caller.py:116`); pinning R1354S tags it.
- **1 × undetermined** — B19896 (`paper_mut=Undetermined`); stays an honest abstention.

## A marker cannot just be added — it must be phenotype-concordant (the D642Y trap)

Before projecting recovery, each candidate marker's R/S split across all 98 isolates is checked: a
marker that also appears in *susceptible* isolates is phenotype-discordant, and tagging it as
resistance manufactures **major errors** (false-R). The per-marker table
(`scripts/characterize_coverage_ceiling.py`) surfaces one:

| marker | R | S | verdict |
|---|---|---|---|
| S639F / S639P | 15 / 7 | 0 / 0 | concordant (already tagged) |
| S639Y | 13 | 1 | mildly discordant — already tagged today; is the existing FP=1 |
| F635C / F635Y | 1 / 1 | 0 / 0 | concordant |
| R1354S | 1 | 0 | concordant (HS2) |
| **W691L** | 2 | 0 | concordant (HS3, CRISPR-validated) |
| W691C | 1 | 0 | concordant **but not literature-validated** → left untagged (abstain) |
| **D642Y** | **2** | **5** | **DISCORDANT** — 5 susceptible carriers at low MIC (CAS 0.5–1); tagging → 5 ME |

**D642Y is the trap.** It is literature-validated as resistance-implicated, so the naïve "widen the
HS1 panel" move would tag it — but here it is genotype–phenotype discordant (2 R / 5 S), and tagging
it pushes ME to **11.5% [5.4, 23.0]**, blowing the ≤5% bar. Trading the VME problem for an ME problem
is not a fix. The honest move is to **leave D642Y untagged**, so the engine abstains
(`UNCHARACTERIZED_VARIANT`) on it — the truthful "a variant is here but resistance is not reliably
callable," consistent with the caller's existing "uncalled is uncalled" discipline.

## Projection under the honest CLEAN panel (a projection, not a measured re-claim)

Tag only markers that are **both** literature-validated **and** phenotype-concordant here —
S639{F,P,Y} (status quo), F635{C,Y}, R1354{S,H}, and the new HS3 **W691{L}**. Feeding the resulting
confusion through the same Wilson math:

| | confusion | VME | ME | abstention |
|---|---|---|---|---|
| **Measured today** | TP36 FN4 FP1 TN46 | 4/40 = **10.0%** [4.0, 23.1] ✗ | 1/47 = 2.1% ✓ | 11/98 = 11.2% ✓ |
| **HS3 window + CLEAN panel** | TP41 FN1 FP1 TN46 | **1/42 = 2.4%** [0.4, 12.3] ✓ | 1/47 = 2.1% ✓ | 9/98 = 9.2% ✓ |

All three bars project to **PASS**. The recovery is: B19617/B19618 (W691L) → true detections;
B19897 F635Y, B21288 F635C, B20592 R1354S → detections; the 7 R-phenotype D642Y carriers and B22769
(W691C, unvalidated) and B19896 (undetermined) remain honest abstentions; **B21978 (M690I, unresolved
even by the source paper) is the one residual miss** — the floor the extension does not clear, and the
reason the projected VME is 1/42, not 0.

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
   with the CRISPR-validated mutant **W691{L}** only (W691C and M690I are *not* independently
   validated → read but left untagged, so their isolates abstain rather than being asserted).
2. **Widen the HS1 panel** to add **F635{C,Y}** (concordant); **pin the HS2 panel** to add
   **R1354{S,H}**. **Do not tag D642Y** — it is phenotype-discordant here and would fail the ME bar.
3. Freeze `work/PREREGISTRATION_diagnostic_accuracy_v2.md` (bar unchanged), then re-run the FKS1
   concordance+accuracy playbook (`work/RUNBOOK_fks1_concordance_and_accuracy.md`) on GCP to produce
   the *measured* post-extension VME/ME and test this projection.

This keeps the caller certified at what it covers (concordance PASS, 100%/100%/100%,
[[fks1-caller-concordance]]) while naming exactly which three windows close the gap — an honest,
tractable Rung C path rather than an off-target dead end. [[diagnostic-accuracy-track]]
