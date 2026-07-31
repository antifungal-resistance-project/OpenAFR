# Run 2 — pre-registered held-out validation

Date: 2026-07-30
Pre-registration: `work/PREREGISTRATION_run2.md` (sha256 7c89d78c..., verified unmodified
by the grading script at run time)
Test set: 7 held-out azole antifungals + 348 property-matched decoys (355 ranked, 2 docking
failures, both decoys)
Protocol: fixed in advance — box 26 A @ 70.61/66.28/4.18, exhaustiveness 32, num_modes 20,
seed 42, rigid 5TZ1 chain A + heme. No per-molecule tuning.

## RESULT: GATE PASSED

Pre-registered bar (unchanged from run 1, not lowered): **AUC >= 0.70 AND EF@1% >= 5.0x**,
evaluated on mode C only.

| mode | AUC | EF@1% | EF@5% | BEDROC | role |
|---|---|---|---|---|---|
| **C — iron-approach distance** | **0.794** | **12.68x** | 5.63x | 0.325 | PRIMARY — **PASS** |
| A — Vina score alone | 0.471 | 0.00x | 0.00x | 0.002 | secondary — worse than random |
| D — combined rank | 0.691 | 0.00x | 2.82x | 0.126 | secondary |

Hypothesis H1, formed on run-1 data and tested here on molecules it had never seen, is
**supported**.

## Statistical support

    permutation test (20,000 label shuffles): p = 0.0028
    bootstrap 95% CI over actives:            AUC 0.590 - 0.946
    active rank positions (of 355):           3, 8, 29, 43, 49, 124, 274

Five of seven held-out drugs land in the top 50 of 355 (top 14%). One lands at 274.

## The headline contrast

On held-out data, the docking program's own affinity score ranks known antifungals
**worse than random** (AUC 0.471). The geometric criterion — how close the molecule can
bring a nitrogen to the heme iron — ranks them at AUC 0.794. Same poses, same run; the
only difference is which number you sort by.

For azoles against CYP51, geometry beats the scoring function. This is now a
pre-registered, held-out result rather than a post-hoc observation.

## Honest limitations

1. **n = 7 actives.** The bootstrap 95% CI lower bound (0.590) sits **below** the 0.70 pass
   bar. The point estimate passes and the permutation test rules out chance, but a true AUC
   under the bar cannot be excluded. More held-out actives would tighten this.
2. **Decoys are presumed inactive, not verified.** Any real binder among the 348 decoys
   depresses the measured performance; any systematic bias in decoy selection inflates it.
3. **One rigid receptor.** 5TZ1 is *C. albicans*, not *C. auris*, and rigid docking cannot
   model induced fit.
4. **Applicability domain is <= 45 heavy atoms**, declared in advance. The
   posaconazole/itraconazole class remains a known blind spot.
5. **This validates a ranking criterion, not a drug.** A high rank is a hypothesis worth
   testing in a wet lab. Nothing here predicts that any molecule will kill a fungus.

## What this earns

Under the pre-registered protocol and within the declared applicability domain, the
pipeline demonstrably recovers known antifungals from property-matched decoys better than
chance and far better than the docking score alone. Screening new molecules is now
justified — with results reported as ranked hypotheses carrying these limitations, never
as hits.

## Reproduce

    ./scripts/run_screen2.sh
    conda run -n openafr --no-capture-output python scripts/validate_gate2.py \
      work/screen2 work/receptor_A.pdb        # exits 0 on pass, 1 on fail
