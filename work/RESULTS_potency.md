# Results — does the N–Fe geometry track CONTINUOUS potency? (#81)

**Verdict: the pre-registered primary hypothesis (H1) FAILS.** On the direct on-target
*C. albicans* CYP51 **enzyme** set, the geometry criterion does **not** rank continuous
potency (Spearman ρ = **+0.103**, one-sided permutation p = 0.72 — not even the predicted
sign). The well-powered but confounded **whole-cell** secondary shows a weak, robustly
significant negative correlation (ρ = **−0.223**, p = 1e-05) that, per the pre-registration,
is *consistent-with* but not *evidence-for* on-target potency tracking. **The honest,
settled conclusion: the geometry criterion enriches actives (binary, AUC ~0.79–0.85) but does
not demonstrably rank how potent a molecule is.** This upgrades nothing; it closes the #81
question with a clean null and one honest caveat-bounded secondary signal.

Frozen inputs (all re-verified unmodified by the grader at run time):

    pre-registration   work/PREREGISTRATION_potency.md   sha256 5d3416d9907ffa9c...
    protocol           protocol.yaml                     sha256 6eca16b1be2e5e33...
    receptor           work/receptor_A.pdb (5TZ1 + heme) sha256 a1a8410868d1...
    datasets           data/ligands/potency_{enzyme,wholecell}.{tsv,smi}  (4 sha256 pinned)

## How to reproduce

```bash
# 1. datasets (network; deterministic build). Already frozen — rebuild only to re-verify.
conda run -n openafr python scripts/make_potency_set.py fetch
conda run -n openafr python scripts/make_potency_set.py build

# 2. dock both sets under the frozen protocol (5TZ1, box 26 Å, exhaustiveness 32, seed 42)
conda run -n openafr python scripts/prep_ligands.py data/ligands/potency_enzyme.smi   work/potency_enzyme_ligands
conda run -n openafr python scripts/prep_ligands.py data/ligands/potency_wholecell.smi work/potency_wholecell_ligands
conda run -n openafr scripts/screen.sh work/potency_enzyme_ligands   work/potency_screen_enzyme
conda run -n openafr scripts/screen.sh work/potency_wholecell_ligands work/potency_screen_wholecell

# 3. grade (re-checks every frozen hash, then correlates)
conda run -n openafr python scripts/validate_gate_potency.py
```

## Primary — enzyme (CHEMBL1780), on-target binding potency

pChEMBL 5.94–8.00 (a narrow ~2-log range; see Limitation 1), n = 32 in-domain compounds,
14/32 of which actually coordinate the iron (closest N–Fe ≤ 3.0 Å).

| subset | n | Spearman ρ | 95% bootstrap CI | one-sided perm p | vs bar (ρ ≤ −0.30 & p < 0.05) |
|---|---:|---:|---|---:|---|
| **H1: full in-domain** | 32 | **+0.103** | [−0.311, +0.477] | 0.715 | **FAIL** |
| leakage-guarded (near-known azoles excluded) | 28 | +0.177 | [−0.267, +0.565] | 0.816 | (reported) |
| iron-bound only (N–Fe ≤ 3.0 Å) | 14 | −0.440 | [−0.824, +0.156] | 0.057 | (reported) |

**Reading it honestly.** The full-set ρ is small and *positive* — the opposite of the
hypothesis — so H1 fails outright, not marginally. Removing the four near-known azoles moves
it further positive, so the failure is not a leakage artifact being masked. The one subset
that trends the predicted way is *iron-bound-only* (ρ = −0.440): **among the 14 molecules that
reach the iron, a tighter reach weakly tracks higher potency**, but its CI crosses zero
(p = 0.057) and — decisively — the pre-registered stopping rule forbids elevating a non-gating
subset to a claim after a failed primary. It is recorded as a hypothesis for a future,
properly-powered, separately-registered test, nothing more.

## Secondary — whole-cell *C. albicans* potency (confounded; does NOT gate)

pChEMBL 4.00–10.97 (a wide range), n = 435 (6 of the 441 in-domain compounds — CHEMBL1200883,
CHEMBL1382917, CHEMBL1450615, CHEMBL1484167, CHEMBL4469221, CHEMBL4592498 — failed PDBQT
conversion in prep with a vina-rejected atom tag; logged, not silently dropped).

| subset | n | Spearman ρ | 95% bootstrap CI | one-sided perm p |
|---|---:|---:|---|---:|
| whole-cell in-domain | 435 | **−0.223** | [−0.313, −0.128] | 1e-05 |

**Reading it honestly.** This *is* a real, negative, highly significant correlation in the
predicted direction — the tight CI (well clear of zero at n=435) rules out chance. But three
things keep it from rescuing H1, exactly as pre-registered:

1. **It is confounded.** Whole-cell potency folds in permeability, efflux and target abundance;
   a whole-cell pMIC is not a CYP51 measurement. The correlation is *consistent with* geometry
   tracking potency, not *evidence for* it.
2. **It is the binary signal in continuous clothing.** The whole-cell set spans chemotypes
   (azole and non-azole); the same "reaches-the-iron azoles are the potent ones" effect that
   drives the validated binary AUC ~0.79 will produce a modest negative ρ here. ρ² ≈ 0.05 — the
   geometry explains ~5% of the rank variance — so this is a weak monotone shadow of the
   chemotype enrichment, not a within-class potency ruler.
3. **The clean on-target test (the primary) says no.** Where the confound is removed (enzyme
   data), the relationship vanishes. That is the more informative of the two results.

## What this means for the tool (the durable takeaway)

- The validated claim stays exactly where it was: **novel-chemotype triage / active-enrichment
  at AUC ~0.79–0.85** (`work/RESULTS_run2_holdout.md`, `work/RESULTS_verified_inactives.md`).
  #81 asked whether it could be upgraded to "predicts potency." It cannot, on the evidence.
- This corroborates the geometry-ceiling story from a new direction: the N–Fe distance reads
  *mechanism-class membership* (does this molecule coordinate the heme iron like an azole?),
  not *fine-grained potency*. The rigid single conformer (Limitation 4) is the most likely
  reason a continuous within-target ordering does not emerge — the same ceiling that closed the
  within-azole ranking line (looks #7–#9, `work/RESULTS_ensemble_confirm.md`).
- A score readout must therefore keep saying what it already says: *a tight N–Fe is a
  hypothesis for a wet lab that this molecule engages CYP51, never a potency estimate.*

## Limitations (as pre-registered, and how each bit here)

1. **Narrow on-target range.** Enzyme pChEMBL spans only 5.94–8.00; a truncated dependent
   variable mechanically shrinks ρ. Real, but not the whole story — the sign is wrong, not just
   the magnitude.
2. **The ≤45-HA domain truncates the potent tail.** Posaconazole/itraconazole (the most potent
   azoles) are out-of-domain and off the plot, capping the achievable range further.
3. **Whole-cell ≠ enzyme** — the secondary's confound, above.
4. **Single rigid conformer, single seed** — the most probable mechanistic reason no continuous
   ordering emerges; a flexible/ensemble receptor is the pre-registered path if this is ever
   re-opened (separate registration required).
5. **pChEMBL heterogeneity** — pooling IC50/Ki/Kd/EC50 across labs adds noise the median only
   partly absorbs.

## Stopping rule (honored)

H1 failed; no subset of these poses is promoted to a potency claim. The iron-bound-only trend
and any flexible-receptor follow-up are future, separately-registered experiments. This file,
`work/PREREGISTRATION_potency.md`, the frozen datasets and `scripts/validate_gate_potency.py`
are the complete record.
