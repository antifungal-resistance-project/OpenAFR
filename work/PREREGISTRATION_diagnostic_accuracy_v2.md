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
re-measure that would test whether the projected **VME 1/40 = 2.5% [0.4, 12.9]** materialises.

## The caller change this run governs (fixed now)

The extended `openafr/fks1_caller.py` may add **only** the following, and nothing else:

1. **A new HS3 `Window`** — residues **688–698**, anchor **W691** (verified to translate to `W` in
   the pinned reference `XM_085597048.1`; the load-time numbering check enforces this, and a
   reference whose residue 691 is not `W` must refuse, exactly as HS1/HS2 do today).
2. **Panel mutant sets** for the resistance positions, added to the relevant windows' `panel`:
   - HS3: **W691 → {L, C, F}**, **M690 → {I}**
   - HS1 (already read): add **D642 → {Y}**, **F635 → {C, Y}** (S639 → {F,P,Y} unchanged)
   - HS2 (already read, currently empty): **R1354 → {S}**

No other window, position, or mutant may be added under this prereg. Widening beyond this set is a
new pre-registration, not a silent edit.

### Anti-circularity rule (load-bearing)

Each mutant letter above **must be independently established as an echinocandin-resistance marker in
the FKS1-resistance literature — NOT taken from this benchmark's `paper_mut` column.** The
pre-registration is void for any position whose resistance status rests only on the isolates being
scored: recovering a marker on the same set that motivated tagging it is circular and does not count
as validation. The PR that implements the change must cite, per position, a source independent of
PMC12323592. Positions that cannot be independently pinned are dropped from the panel (and the
affected isolates remain honest abstentions), even if that leaves the VME above bar.

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

On the frozen fixture, the conservative projection is **VME 1/40 = 2.5% [0.4%, 12.9%]** (only B21978,
M690I unresolved-by-source, remains a miss). Recording it now makes the re-measure falsifiable: a
measured post-extension VME materially worse than this — e.g. if the reads do not actually cover HS3
at depth, converting expected detections into `partial`/abstentions rather than calls — is a
reportable miss of the prediction, not a silently-revised expectation.

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
