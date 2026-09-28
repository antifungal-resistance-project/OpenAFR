# Diagnostic verdict-accuracy vs phenotypic AST (#137) — echinocandin arm, measured

**The engine's verdict accuracy, measured against phenotype.** When the resistance-interpretation
engine emits a categorical verdict for a *C. auris* isolate, how often does it agree with that
isolate's echinocandin R/S phenotype, in CLSI-M23 error terms (very-major / major error)? This is
Rung A of the clinical gate ladder (`docs/DIAGNOSTIC_GO_NO_GO.md`) — until now **blocked on paired
genotype+phenotype data** and expected to read underpowered. The FKS1 concordance run
([`RESULTS_fks1_concordance.md`](RESULTS_fks1_concordance.md)) unblocked the echinocandin arm: its
reads→caller pass over the 98-isolate PMC12323592 benchmark **harvested the paired fixture in the
same session**, at zero extra download.

**Pre-registration:** `work/PREREGISTRATION_diagnostic_accuracy.md` (sha `70d5d2ae…`), bar frozen
before any isolate was scored. **Harvested fixture:**
`data/earlywarning/recaller_sanity/fks1_accuracy_98.tsv` (sha `1aeffea5…`, pinned in
`work/PREREG_diagnostic_accuracy.sha256`) — its `called_verdict` column is
`echinocandin_verdict(call=...)` emitted from the live caller output per isolate (**filled by the
orchestration, not hand-authored**, per the prereg rule); `susceptibility` copied verbatim from the
frozen concordance truth. **Code:** `openafr/accuracy.py`, `scripts/validate_diagnostic_accuracy.py`.
Both integrity hashes matched at run time.

> **UPDATE 2026-09-27 — a v3 re-tier (below) supersedes the interpretation of both v1 and v2.** The
> v2 and v1 sections are preserved verbatim as the auditable priors. Read v3 first.

---

# v3 re-tier — high-PPV FKS1 panel tier (#137), measured → PASS

**What changed.** The v2 FAIL was a **marker-PPV ceiling**, not a coverage gap: the panel detected
resistance (VME → 0%) but over-called it (ME 13.5%) because two panel positions — **D642Y** and
**M690I** — are known FKS1 changes that *do not predict the phenotype* (D642Y splits 2 R / 5 S; M690I
1 R / 1 S). v3 adds a **PPV tier** over the *same* v2 panel (no new markers): only high-PPV **core**
markers (F635/S639 HS1, R1354S HS2, W691L HS3) emit `RESISTANCE_MARKER_DETECTED`; a low-PPV-only call
(D642Y/M690I) emits `UNCHARACTERIZED_VARIANT` — an honest abstain, never resistance and never
susceptible. Frozen in `work/PREREGISTRATION_diagnostic_accuracy_v3.md` (sha pinned in
`work/PREREG_diagnostic_accuracy_v3.sha256`), with a pre-committed falsifiable projection of the exact
re-tiered confusion. **Code:** the PPV tier (`openafr.fks1_caller.panel_tier` / `core_panel_hits`) and
the tier-aware verdict (`openafr.verdict._classify`). **Fixture:**
`data/earlywarning/recaller_sanity/fks1_accuracy_98_v3.tsv` (sha `02fc37e1…`) — **derived, not
re-read**: `scripts/harvest_diagnostic_accuracy_v3.py` re-emits each `called_verdict` by running the
production `echinocandin_verdict()` over the caller output *already frozen in the v2 fixture*, so no SRA
re-read/GCP was needed. The v1 and v2 fixtures are preserved un-overwritten. All prereg + fixture
hashes matched at grade time.

## Headline — the projection held exactly: Rung A PASS for the narrowed high-PPV claim

> Scored **88** isolates (**42 R / 46 S**); **10 abstained**, 0 unresolved, 0 no-phenotype — all excluded.
> Confusion (R/S × call): **TP=42 FN=0 FP=1 TN=45**.

| Metric | Value | Wilson 95% | Bar | Verdict | v2 → v3 |
|---|---|---|---|---|---|
| **Very-major error** (missed R) | **0/42 = 0.0%** | [0.0%, 8.4%] | point ≤3% **and** upper ≤15% | ✓ | 0.0% → **0.0%** |
| **Major error** (false R) | **1/46 = 2.2%** | [0.4%, 11.3%] | ≤5% | ✓ | 13.5% → **2.2%** |
| Sensitivity (= 1 − VME) | 42/42 = 100.0% | [91.6%, 100%] | — | — | 100% → 100% |
| Specificity (= 1 − ME) | 45/46 = 97.8% | [88.7%, 99.6%] | — | — | 86.5% → 97.8% |
| Categorical agreement | 87/88 = 98.9% | [93.8%, 99.8%] | — | — | 92.8% → 98.9% |
| Abstention (uncharacterised) | 10/98 = 10.2% | [5.6%, 17.8%] | ≤30% | ✓ | 1.0% → 10.2% |

**Pre-registered verdict: PASS** — every gated metric clears its bar, matching the frozen projection
(VME 0/42, ME 1/46, abstention 10/98) exactly.

## What it means — a defensible, honestly-narrowed intended use

The 9 low-PPV-only carriers (7 D642Y, 2 M690I) move from the confusion matrix into **abstention**:

- The **3 resistant** low-PPV-only carriers (D642Y-only B20717/B21114; M690I B21978) are now honest
  `UNCHARACTERIZED_VARIANT` — **not very-major misses**. Abstaining is the contract-correct outcome
  (rule 3): the engine says "a known-but-low-PPV change is present; I cannot call resistance from it,"
  which is neither a false-negative nor a susceptibility claim.
- The **6 susceptible** low-PPV-only carriers (5 D642Y, M690I B20326) likewise abstain instead of being
  over-called R — this is what fixes the v2 major-error failure.
- Only **B20673** (S639Y, phenotypically S) remains a single major error — irreducible genotype↔phenotype
  noise at an otherwise-reliable core position, well within the ≤5% bar.

**The claim this licenses (and no more):** *detection over validated high-PPV FKS1 markers
(F635/S639/R1354S/W691L core), abstaining on low-PPV positions (D642Y/M690I).* The 10.2% abstention is
the honest price of that narrowing, reported openly. This is **not** a claim to classify every isolate
R/S — the v2 result already showed no FKS1-genotype panel can do that on this cohort — it is a claim
that *when the engine detects a core marker, that detection is accurate*, and that it abstains rather
than guess elsewhere.

**Anti-circularity.** The tier is pinned to evidence independent of the scored isolates' own phenotype:
D642Y is demoted because the PMC12323592 authors' *own expert panel* scores it wild-type (`-`) and CDC
EID 25-0760 does not treat it as resistance-defining; M690I because it is phenotype-unresolved in the
source literature. The core tier is the canonical Perlin FKS1 hotspot set (W691L additionally
CRISPR-confirmed). The rule is "tier by independent literature/expert-panel status," frozen before the
re-grade.

## Bottom line for the clinical gate

Rung A is now a **measured GO for echinocandin under the narrowed high-PPV claim** — the first arm to
clear the frozen bar. It does **not** lift the broader limits: the engine still makes no susceptibility
call, no probability (#136), and no coverage of efflux/ERG3/TR or *C. albicans*; the azole/ERG11 arm
stays UNDERPOWERED. What v3 establishes is that the *detection* claim, honestly scoped to high-PPV
markers, is accurate and defensible — the over-call ceiling was a claim-scope problem, and scoping the
claim correctly resolves it.

---

# v2 re-measure — FKS1 HS3 + widened HS1/HS2 panel (#154), measured

**What changed.** The v1 FAIL was traced to marker coverage (all 4 very-major misses were FKS1 **HS3**,
a hotspot the caller had no window for). PR #154 added an HS3 window (anchor W691) and widened the
HS1/HS2 panel tags, exactly as frozen in `work/PREREGISTRATION_diagnostic_accuracy_v2.md`
(sha `014094c6…`), with a **pre-committed falsifiable projection of VME 1/40 = 2.5% [0.4, 12.9]**. A
fresh GCP reads→caller pass over the same 98-isolate PMC12323592 benchmark re-harvested the paired
fixture (`data/earlywarning/recaller_sanity/fks1_accuracy_98_v2.tsv`, sha `a7c5056f…`, pinned in
`work/PREREG_diagnostic_accuracy_v2.sha256`; the v1 fixture is preserved un-overwritten). Both prereg
hashes and the re-harvested fixture hash matched at grade time.

## Headline — the projection held, but the failure mode flipped: FAIL on major error

> Scored **97** isolates (**45 R / 52 S**); 1 abstained, 0 unresolved, 0 no-phenotype — all excluded.
> Confusion (R/S × call): **TP=45 FN=0 FP=7 TN=45**.

| Metric | Value | Wilson 95% | Bar | Verdict | v1 → v2 |
|---|---|---|---|---|---|
| **Very-major error** (missed R) | **0/45 = 0.0%** | [0.0%, 7.9%] | point ≤3% **and** upper ≤15% | ✓ | 10.0% → **0.0%** |
| **Major error** (false R) | **7/52 = 13.5%** | [6.7%, 25.3%] | ≤5% | ✗ **BELOW BAR** | 2.1% → **13.5%** |
| Sensitivity (= 1 − VME) | 45/45 = 100.0% | [92.1%, 100%] | — | — | 90.0% → 100% |
| Specificity (= 1 − ME) | 45/52 = 86.5% | [74.7%, 93.3%] | — | — | 97.9% → 86.5% |
| Categorical agreement | 90/97 = 92.8% | [85.8%, 96.5%] | — | — | 94.3% → 92.8% |
| Abstention (uncharacterised) | 1/98 = 1.0% | [0.2%, 5.6%] | ≤30% | ✓ | 11.2% → 1.0% |

**Pre-registered verdict: FAIL** — but on the *opposite* axis from v1. The HS3 extension did exactly
what the projection predicted: **VME collapsed 10% → 0%** (the two W691L misses B19617/B19618 and the
M690I miss B21978 are now correctly detected; the pre-committed VME projection of 2.5% materialised,
in fact better). But the widened HS1 panel over-calls: **ME rose 2.1% → 13.5%**, above the ≤5% bar.

## What it means — D642Y is not a reliable resistance marker in this cohort (a biological ceiling)

The 7 major errors are **genotype-correct but phenotype-discordant** — the caller read each genome
right; the tagged position just does not predict the phenotype:

- **5 of 7 are D642Y** (B20464, B20681, B20689, B20702, B20704), all phenotypically **susceptible**.
  Across the whole benchmark, D642Y carriers split **2 R / 5 S** — a ~29% positive predictive value.
  There is **no genotype rule that wins**: tag D642Y and you get 5 major errors; drop it and its 2
  resistant carriers (B20717, B21114, which carry *only* D642Y) become very-major misses (VME → 4.4%,
  which also fails). Genotype alone cannot separate D642Y's R from S carriers.
- **1 M690I** (B20326, S) — the same HS3 position that recovered a true VME miss (B21978, R) also
  over-calls one susceptible carrier: M690I splits 1 R / 1 S here.
- **1 S639Y** (B20673, S) — a single discordant carrier of an otherwise-reliable core marker (every
  other S639F/P/Y isolate is R); irreducible genotype↔phenotype noise, not a panel defect.

**Convergent, independent falsification of D642Y.** The v2 anti-circularity rule pinned D642Y to
literature *independent* of this benchmark (CDC EID 25-0760). The benchmark then independently rejects
it on **two** axes: (a) its expert `expected_panel` column lists **`-` (wild-type)** for every D642Y
carrier — the authors do not count D642Y as a resistance panel marker — and (b) 5/7 carriers are
phenotypically susceptible. This is the anti-circularity discipline working as designed: an
independently-pinned marker, tested against held-out truth, failed.

**Concordance footnote (not a regression of existing calls).** Graded against the *narrow* S639-only
`expected_panel`, the v2 caller's token specificity drops to 75.4% (15 FP) — but every one of those
"FP" is a v2-*new* marker (D642Y/F635/W691L/R1354S/M690I) that the benchmark's narrow panel scores as
`-`; **all 36 original S639 calls remain 100% concordant**. The v2 caller is a broader panel than the
concordance truth set, so that specificity number reflects truth-set scope, not a broken call. The
clinically meaningful judgment is the **phenotype** arm above.

## Bottom line for the clinical gate

Rung A stays a **measured NO-GO** for echinocandin — but the reason has moved from *under-detection*
(v1 VME) to *over-calling* (v2 ME), and the over-call is a **biological ceiling, not a coverage gap**:
on this cohort no FKS1-genotype panel clears both the VME and ME bars simultaneously, because D642Y
(and, weakly, M690I) carriers are phenotypically split. More markers cannot fix this; it is the limit
of a genotype→verdict engine against phenotype for these positions. The honest intended-use claim
narrows accordingly (detection-only over *validated, high-PPV* markers — i.e. the S639/W691 core —
explicitly excluding low-PPV positions like D642Y), which is a **new pre-registration**, not a silent
panel edit.

---

# v1 (superseded) — narrow S639-only panel

## Headline — echinocandin arm FAILs the clinical VME bar

> Scored **87** isolates (**40 R / 47 S**); 11 abstained (`UNCHARACTERIZED_VARIANT`), 0 unresolved,
> 0 no-phenotype — all excluded, never scored as a negative.
> Confusion (R/S × call): **TP=36 FN=4 FP=1 TN=46**.

| Metric | Value | Wilson 95% | Bar | Verdict |
|---|---|---|---|---|
| **Very-major error** (missed R) | 4/40 = 10.0% | [4.0%, 23.1%] | point ≤3% **and** upper ≤15% | ✗ **BELOW BAR** |
| **Major error** (false R) | 1/47 = 2.1% | [0.4%, 11.1%] | ≤5% | ✓ |
| Sensitivity (= 1 − VME) | 36/40 = 90.0% | [76.9%, 96.0%] | — | — |
| Specificity (= 1 − ME) | 46/47 = 97.9% | [88.9%, 99.6%] | — | — |
| Categorical agreement | 82/87 = 94.3% | [87.2%, 97.5%] | — | — |
| Abstention (uncharacterised) | 11/98 = 11.2% | [6.4%, 19.0%] | ≤30% | ✓ |

**Pre-registered verdict: FAIL** — VME point 10.0% (bar ≤3%) and VME Wilson upper 23.1%
(bar ≤15%) are both above bar. ME and abstention pass.

## What it means

1. **An honest FAIL, not a caller defect.** The concordance run certified the caller's *tokens* as
   perfect (sensitivity/specificity/identity all 100%). The accuracy FAIL is therefore **not** the
   caller misreading the genome — it is the **narrow-panel mechanism-coverage ceiling** (go/no-go
   C2) made quantitative: 4 of 40 phenotypically-resistant isolates carry resistance that is not an
   HS1 S639 substitution, so the engine correctly returns `NO_KNOWN_MARKER`/`UNCHARACTERIZED_VARIANT`
   and, scored against phenotype, misses them. This is the same signal as the concordance
   resistance-detection ceiling (36/46 = 78.3%), viewed through the CLSI error frame.
2. **Detection-only ≠ diagnostic — now measured, not asserted.** The engine's day-one RUO claim
   never promised susceptibility calls; C2 always stated that `NO_KNOWN_MARKER` is not "susceptible."
   This run replaces that assertion with a measured VME of 10% [4.0, 23.1] — the exact number the
   go/no-go doc pre-committed to reporting as the mechanism-coverage ceiling.
3. **Rung A is a measured NO-GO for echinocandin, not a blocked one.** The gate ran and the arm did
   not clear it. Lifting it needs *broader marker coverage* (Rung C: efflux / non-HS1 FKS1 hotspots),
   not more compute — the caller is already perfect at what it covers.
4. **The azole/ERG11 arm remains separately underpowered**, and **calibration (#136) remains
   blocked** on the non-public option-C paired *C. albicans* panel. This session unblocked only the
   echinocandin arm; those ceilings are unchanged. [[calibration-track-blocked]]
