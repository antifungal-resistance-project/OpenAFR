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
