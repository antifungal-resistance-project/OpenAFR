# Pre-registration v2 — echinocandin verdict accuracy after the FKS1 panel extension (#137)

> **STATUS: FROZEN 2026-09-20, before the panel-extension caller change is written and before any
> post-extension verdict has been scored.** This amends `work/PREREGISTRATION_diagnostic_accuracy.md`
> (the v1 prereg, sha `70d5d2ae…`) for a *specific, named* caller change and re-measure. It fixes, in
> advance: the exact windows/panel the extended caller may add, the anti-circularity rule the panel
> mutant set must obey, the (unchanged) pass bar, and the interpretation of both outcomes. The v1
> mapping, metric definitions, and CLSI denominators are inherited verbatim and NOT re-opened. This
> document's SHA-256 is pinned in `work/PREREG_diagnostic_accuracy_v2.sha256`; the post-extension
> grader will re-check it and the (re-harvested) fixture hash and refuse to certify if either moved.

## Why a v2 is required — the measured FAIL has a characterised, tractable cause

The v1 echinocandin arm FAILed: **VME 4/40 = 10.0% [4.0%, 23.1%]**
(`work/RESULTS_diagnostic_accuracy.md`). The coverage characterisation
(`work/RESULTS_diagnostic_coverage_ceiling.md`, `scripts/characterize_coverage_ceiling.py`) showed the
FAIL is **not** off-target/efflux resistance: all 4 misses are **FKS1 HS3** (B19617/B19618 W691L,
B22769 W691C; B21978 M690I, unresolved even by the source paper), a hotspot the caller has no window
for. The 11 abstentions are mostly **HS1 substitutions the caller already reads but does not tag**
(9 × D642Y/F635, 1 × HS2 R1354S). This v2 governs the caller change that would close that gap and the
re-measure that would test whether the projected **VME 1/42 = 2.4% [0.4, 12.3]** (with ME held at
2.1% by *excluding* the discordant D642Y — see Gate 2 below) materialises.

## The caller change this run governs (fixed now)

The extended `openafr/fks1_caller.py` may add **only** the following, and nothing else:

1. **A new HS3 `Window`** — residues **688–698**, anchor **W691** (verified to translate to `W` in
   the pinned reference `XM_085597048.1`; the load-time numbering check enforces this, and a
   reference whose residue 691 is not `W` must refuse, exactly as HS1/HS2 do today).
2. **Panel mutant sets** — only markers that pass **both** gates below:
   - HS3: **W691 → {L}**  (W691L is CRISPR/Cas9-validated in *C. auris*, AAC 2023
     doi:10.1128/aac.00423-23 / PMC10269051. **W691C and M690I are NOT independently validated** —
     they are read by the new window but left **untagged**, so their isolates abstain rather than
     being asserted resistant.)
   - HS1 (already read): add **F635 → {C, Y}**  (S639 → {F,P,Y} unchanged).
   - HS2 (already read, currently empty): **R1354 → {S, H}**  (R1354H CRISPR-validated, PMC10219442).

No other window, position, or mutant may be added under this prereg. Widening beyond this set is a
new pre-registration, not a silent edit.

### Two gates every tagged mutant must pass (load-bearing)

**Gate 1 — independent literature validation.** Each tagged mutant must be established as an
echinocandin-resistance marker in literature **independent of PMC12323592** (AAC 2022
doi:10.1128/aac.01243-22, AAC 2023 doi:10.1128/aac.00423-23, PMC10219442 for R1354H). Recovering a
marker on the same set that motivated tagging it is circular; a position that cannot be independently
pinned is left untagged (its isolates abstain), even if that leaves the VME above bar. **W691C and
M690I fail this gate** and are therefore untagged.

**Gate 2 — phenotype concordance (the D642Y exclusion).** A literature-validated marker that also
appears in **susceptible** isolates in the benchmark is phenotype-discordant; tagging it manufactures
major errors. `scripts/characterize_coverage_ceiling.py` reports the per-marker R/S split, and
**D642Y is discordant (2 R / 5 S, the S carriers at low MIC CAS 0.5–1)** — tagging it drives ME to
11.5% [5.4, 23.0], failing the ≤5% bar. **D642Y is therefore NOT tagged**; the engine abstains on it
(`UNCHARACTERIZED_VARIANT`), the honest "variant seen, resistance not reliably callable." Excluding a
discordant marker is *conservative* (it removes a false-R source, it does not inflate recovery), and
the discordance itself is a reported finding — not a silent bar-tune. If the eventual GCP re-measure
were to *include* D642Y, this prereg is violated.

### Scope guard — the ripple is acknowledged, not hidden

Adding a panel-bearing HS3 window changes `backtest._FKS1_PANEL_WINDOWS` (derived from `w.panel`), the
`format_call` per-window source strings (`sra-fks1-recaller:HS1=…,HS2=…,HS3=…`), the prevalence
resolution semantics, and stored snapshots. The implementing PR must update the affected tests
(`tests/test_fks1_caller.py`, `test_fks1_prevalence.py`, `test_recall_fks1_orchestration.py`,
`test_earlywarning_snapshot.py`, `test_interpret.py`, `test_accuracy_harvest.py`) and re-freeze any
snapshot fixtures — this is why the change is a separate reviewed PR, not folded into the
characterisation.

## The re-measure protocol (fixed now)

The post-extension verdicts must be **re-harvested from a live caller pass**, not edited into the
existing fixture. Following `work/RUNBOOK_fks1_concordance_and_accuracy.md`, a GCP session
re-downloads the PMC12323592 reads, re-runs the concordance grader (the extended caller must still
PASS token concordance — a new window must not break existing HS1/HS2 calls), and re-harvests
`fks1_accuracy_98.tsv`'s `called_verdict` column from the extended caller. The v1 fixture is
preserved (renamed/versioned), never overwritten in place, so the FAIL remains auditable.

## Pass bar (inherited from v1, unchanged)

    PASS         if  VME point <= 0.03  AND  VME Wilson UPPER <= 0.15
                 AND ME point <= 0.05   AND  abstention <= 0.30
    FAIL         if any gated metric's point is above its bar
    UNDERPOWERED if the resistant arm < 15 scored OR the susceptible arm < 15 scored

The bar is deliberately identical to v1: the panel extension must clear the *same* bar the narrow
panel failed, or the extension does not license a Rung-A accuracy claim.

## Pre-committed projection (the falsifiable prediction)

On the frozen fixture, under the CLEAN panel above, the projection is **VME 1/42 = 2.4% [0.4%,
12.3%]**, **ME 1/47 = 2.1% [0.4%, 11.1%]**, **abstention 9/98 = 9.2% [4.9%, 16.5%]** — all three bars
project to PASS (`scripts/characterize_coverage_ceiling.py`). The single residual VME is **B21978**
(M690I, unresolved even by the source paper). Recording these now makes the re-measure falsifiable in
both directions: a measured post-extension **VME** materially worse than 2.4% (e.g. if reads do not
cover HS3 at depth, turning expected W691L detections into `partial`/abstentions) **or** a measured
**ME** above 5% (e.g. if a tagged marker turns out discordant on the live calls) is a reportable miss
of the prediction, not a silently-revised expectation.

## Interpretation, committed now (both directions)

- **PASS.** The extended panel clears the same bar the narrow panel failed; the echinocandin arm earns
  a measured RUO accuracy claim, and the Rung-A no-go is lifted for echinocandin. Folded into
  `RESULTS_diagnostic_accuracy.md` (v2) and the go/no-go register.
- **FAIL (extension insufficient or reads under-cover HS3).** A genuine negative: the projected
  recovery did not materialise. Reported with which isolates stayed missed and why (coverage vs.
  mechanism), and the narrow-panel intended-use decision (Rung C, the other fork) is revisited.
- **UNDERPOWERED** is not expected here (the echinocandin arm already scored ≥15/arm at v1) but is
  retained for completeness.

## Not in scope (declared)

- Efflux (TAC1/MRR1/CDR1), ERG3, and promoter/TR mechanisms — still uncovered; a resistant isolate
  carrying only those remains an honest miss, and this extension makes no claim over them.
- *C. albicans* / cross-species breadth (Rung C species scope) and the probability field (#136).
- The azole/ERG11 arm, which stays separately UNDERPOWERED under the v1 prereg.
