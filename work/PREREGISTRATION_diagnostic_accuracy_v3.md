# Pre-registration v3 — echinocandin verdict accuracy under a high-PPV panel tier (#137)

> **STATUS: FROZEN 2026-09-27, before the PPV-tier code change is written and before any re-tiered
> verdict has been scored.** This amends `work/PREREGISTRATION_diagnostic_accuracy.md` (v1, sha
> `70d5d2ae…`) and `work/PREREGISTRATION_diagnostic_accuracy_v2.md` (v2). It fixes, in advance: the
> exact PPV tiering the caller/verdict may apply, the anti-circularity rule that tiering must obey, the
> (unchanged) pass bar, the re-derivation protocol, and the interpretation of both outcomes. The v1
> mapping, metric definitions, and CLSI denominators are inherited verbatim and NOT re-opened. This
> document's SHA-256 is pinned in `work/PREREG_diagnostic_accuracy_v3.sha256`; the re-tiered grader
> re-checks it and the derived fixture hash and refuses to certify if either moved.

## Why a v3 is required — the v2 FAIL is a marker-PPV ceiling, not a coverage gap

The v2 echinocandin arm (HS3 window + widened HS1/HS2, #154) did exactly what its projection predicted
for **detection**: **VME collapsed 10% → 0.0% [0.0, 7.9], sensitivity 100%**. But the failure mode
**flipped** — **major error rose 2.1% → 13.5% [6.7, 25.3]** (bar ≤5%), a FAIL on *over-calling*
(`work/RESULTS_diagnostic_accuracy.md`). The cause is **marker PPV, not coverage**: two panel positions
are known FKS1 changes but **do not predict the phenotype**:

- **D642Y** — carriers split **2 R / 5 S** in the benchmark (~29% PPV). Independently, the benchmark's
  own expert `expected_panel` scores every D642Y carrier as `-` (wild-type), and the CDC EID 25-0760
  collection lists D642Y among HS1 changes without treating it as resistance-defining.
- **M690I** — splits **1 R / 1 S** here and is phenotype-uncertain/unresolved in the source literature.

More markers cannot fix this: dropping D642Y silently re-opens VME (its 2 R-only carriers become
misses); tagging it over-calls. The honest move is a **two-tier panel** that DETECTS high-PPV markers
and ABSTAINS on low-PPV positions — reported as `UNCHARACTERIZED_VARIANT`, never as susceptible and
never as resistance. This is the "re-scoped, high-PPV intended-use claim, pre-registered fresh" the v2
result and `docs/DIAGNOSTIC_GO_NO_GO.md` (Rung A/C) committed to.

## The code change this run governs (fixed now)

The change may add **only** a PPV tier over the *existing* v2 panel — no new window, position, or
mutant letter, and no change to what `is_panel_token`/prevalence counts as a panel token:

1. **`openafr/fks1_caller.py`** defines a `LOW_PPV_PANEL = {642: {Y}, 690: {I}}` and exposes
   `panel_tier(token) -> 'core' | 'low_ppv' | None` and `is_core_panel_token`. `call_window` /
   `call_windows` additionally emit `core_panel_hits` (the high-PPV subset of `panel_hits`). Every
   token is still emitted and panel-tagged exactly as in v2.
2. **`openafr/verdict.py`** `_classify` drives `RESISTANCE_MARKER_DETECTED` off the caller's
   `core_panel_hits` when present (FKS1); a panel hit that is present but **not** core falls through to
   `UNCHARACTERIZED_VARIANT`. ERG11 exposes no such key and is unchanged (every panel hit detects).

The **CORE (high-PPV) tier** is exactly: **F635{C,Y}, S639{F,P,Y}** (HS1); **R1354{S}** (HS2);
**W691{L}** (HS3). The **LOW-PPV tier** is exactly **D642Y** and **M690I**. No other position moves
tiers under this prereg.

### Anti-circularity rule (load-bearing)

A position's tier **must rest on evidence independent of the PMC12323592 isolates being scored.** The
demotions are pinned so:

- **D642Y → low-PPV:** the PMC12323592 authors' *own expert panel* scores D642Y carriers as `-`
  (wild-type) — an independent expert judgement, distinct from the per-isolate R/S phenotype this run
  grades against — corroborated by CDC EID 25-0760 not treating it as resistance-defining.
- **M690I → low-PPV:** reported as a single / phenotype-unresolved association in the source
  literature, i.e. not an established resistance marker.

The core tier is the set of canonical Perlin FKS1 hotspot substitutions with consistent R association
(HS1 S639/F635; HS2 R1354S; HS3 W691L, additionally CRISPR-confirmed, AAC 2023 aac.00423-23). The rule
is **"tier by independent literature / expert-panel status,"** stated as a principle — not a
post-hoc pick of this cohort's discordant rows. A tier assignment that could only be justified by the
scored isolates' own phenotypes is void.

## The re-derivation protocol (fixed now)

No SRA re-read is required or permitted for v3: the panel re-tiering only re-disposes calls the caller
already made. `scripts/harvest_diagnostic_accuracy_v3.py` re-derives each isolate's `called_verdict` by
running the **production** `echinocandin_verdict()` over the caller output **already frozen in the v2
fixture** (the per-isolate `panel_hits` recorded verbatim in its `note` column). The v2 fixture is
preserved un-overwritten; the derived `data/earlywarning/recaller_sanity/fks1_accuracy_98_v3.tsv` is a
new file whose hash is pinned. The harvester asserts the invariant that a verdict may change **only**
from `RESISTANCE_MARKER_DETECTED` → `UNCHARACTERIZED_VARIANT`, and only for a low-PPV-only call.

## Pass bar (inherited from v1/v2, unchanged)

    PASS         if  VME point <= 0.03  AND  VME Wilson UPPER <= 0.15
                 AND ME point <= 0.05   AND  abstention <= 0.30
    FAIL         if any gated metric's point is above its bar
    UNDERPOWERED if the resistant arm < 15 scored OR the susceptible arm < 15 scored

The bar is identical to v1/v2: the re-tiered panel must clear the *same* bar, or the tiering does not
license a Rung-A accuracy claim.

## Pre-committed projection (the falsifiable prediction)

Re-tiering moves the 9 low-PPV-only isolates (7 × D642Y, 2 × M690I) from the confusion matrix to
abstention: **3 resistant** (D642Y-only B20717/B21114; M690I B21978) leave TP → abstain (NOT counted
as very-major misses — they are honest `UNCHARACTERIZED_VARIANT`), and **6 susceptible** (5 × D642Y,
M690I B20326) leave FP → abstain. The projected v2→v3 confusion is **TP=42 FN=0 FP=1 TN=45**:

- **VME 0/42 = 0.0%** (Wilson upper ≈ 8.4% ≤ 15%) — ✓
- **ME 1/46 = 2.2%** (only the single genuinely-discordant S639Y-S carrier B20673 remains) ≤ 5% — ✓
- **abstention 10/98 = 10.2%** ≤ 30% — ✓; **n_res 42 ≥ 15, n_sus 46 ≥ 15**

→ **projected PASS.** Recording it now makes the re-derivation falsifiable: a measured result that
differs (e.g. a core marker unexpectedly demoted, or a residual FP/FN) is a reportable miss, not a
silently-revised expectation.

## Interpretation, committed now (both directions)

- **PASS.** The high-PPV core panel clears the bar the v2 panel failed. The echinocandin arm earns a
  measured RUO accuracy claim **for the narrowed intended use** — *detection over validated high-PPV
  FKS1 markers (S639/F635/W691/R1354S core), abstaining on low-PPV positions (D642Y/M690I)* — and
  Rung A is lifted **for that narrowed claim only**. Folded into `RESULTS_diagnostic_accuracy.md` (v3)
  and the go/no-go register. The abstention rate is reported as the honest price of the narrowing.
- **FAIL.** The tiering did not clear the bar (e.g. a remaining FP/FN the projection missed). A genuine
  negative: reported with which isolates and why, and the intended-use claim is revisited rather than
  re-tiered again to fit.
- **UNDERPOWERED** is not expected (≥15/arm survives the re-tiering) but is retained for completeness.

## Not in scope (declared)

- Any change to the *panel membership* or windows (frozen at v2) — only the PPV tier is added.
- Efflux (TAC1/MRR1/CDR1), ERG3, promoter/TR mechanisms; *C. albicans* / cross-species breadth; and
  the probability field (#136) — all still out of scope.
- The azole/ERG11 arm, which stays separately UNDERPOWERED.
