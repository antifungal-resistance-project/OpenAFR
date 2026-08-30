# Pre-registration — does the N–Fe geometry track CONTINUOUS potency? (#81)

**Status: FROZEN before any geometry (closest N–Fe distance) was computed for these
compounds.** The two datasets below were built first, from ChEMBL, by a deterministic
label-only pipeline (`scripts/make_potency_set.py build`) that never touches a docked pose —
it reads `pchembl_value` (the dependent variable) and the applicability triage (heavy-atom
count + Tanimoto, both pose-free). Their SHA-256 are pinned here; this document is then frozen
by its own SHA-256 (`work/PREREG_potency.sha256`). The grader
(`scripts/validate_gate_potency.py`) re-checks this document's hash, the `protocol.yaml` hash,
the receptor-anchor hash and the two dataset hashes before it will report a correlation, and
refuses to certify a result if any moved. **No N–Fe distance for any potency-set compound had
been read at freeze time.**

## The question (#81)

Every validated OpenAFR headline is *binary*: the N–Fe geometry criterion **enriches actives**
over inactives (AUC ~0.79–0.85; `work/RESULTS_run2_holdout.md`, `work/RESULTS_verified_inactives.md`).
A wet lab would bet harder on a *monotonic* claim — that a **tighter** nitrogen-to-heme-iron
approach goes with a **higher measured potency**, not merely with the "active" label. This is
the one pre-registered test of that stronger claim.

## The geometry (independent variable) — UNCHANGED, the validated criterion

The criterion is the project's primary metric, computed exactly as `scripts/validate_gate2.py`
computes it: **mode C = the minimum ligand-nitrogen-to-heme-iron distance** over the docked
poses, under the hash-frozen protocol and the validated 5TZ1 receptor. Nothing new is invented.

    receptor     work/receptor_A.pdb  (5TZ1 chain A + heme)
                 sha256 a1a8410868d1b9ee2835117cf6327ba4721ca3830a4d5689b7e9ea0ee1eee22f
    protocol     protocol.yaml  sha256 6eca16b1be2e5e33ae42d4b1b7a63a0e9b891a5dd6cbfb19642dd3e9a149aba4
                 box 26 Å @ (70.61, 66.28, 4.18), exhaustiveness 32, 20 modes, seed 42
    prep         scripts/prep_ligands.py → screen.sh → analyze_poses, the score_molecule.py path

Smaller closest N–Fe = tighter iron approach = the geometry's "more like a real azole" end.

## The potency (dependent variable) — ChEMBL pChEMBL, one comparable −log scale

`pchembl_value` = −log10(molar potency), folding IC50/Ki/Kd/EC50/Potency onto one axis
(pChEMBL 7 = 100 nM, 8 = 10 nM; higher = more potent). It exists only for exact ('=')
molar-standardised records, so censored ('>'/'<') and raw-MIC records carry none and are
dropped by construction (`openafr.potency.median_pchembl`, `POTENCY_TYPES`). Each compound is
collapsed to **one median pChEMBL** across its qualifying records, so no compound is weighted
by how often it was assayed. MIC is deliberately **not** on this axis: ChEMBL emits no
pchembl_value for raw MIC, so it never enters the aggregation.

## The two datasets (frozen)

Built by `scripts/make_potency_set.py`; only compounds inside the validated applicability
envelope (≥1 aromatic N, ≤45 heavy atoms — `scripts.applicability_triage`) are docked, because
the method makes **no claim** outside it. Out-of-domain compounds are excluded, not scored.

| set | ChEMBL selector | in-domain n | pChEMBL range | role |
|---|---|---:|---|---|
| **enzyme** | target **CHEMBL1780** (*C. albicans* CYP51 / sterol 14α-demethylase) | **32** | 5.94 – 8.00 | **PRIMARY** — direct on-target binding potency, no whole-cell confound |
| **wholecell** | target_organism *Candida albicans* | **441** | 4.00 – 10.97 | SECONDARY — well-powered but confounded (uptake / efflux / target abundance) |

    data/ligands/potency_enzyme.tsv     sha256 57de757b8c6644424cfea7b99d611998ff4bdd58b983d0555c0b5c1384cc96bc
    data/ligands/potency_enzyme.smi     sha256 bf38daaa98d997db0f0b62d394443e5540ebb7f9fc071502b3912386d61e0d5b
    data/ligands/potency_wholecell.tsv  sha256 03e2d4fa2dc92f0600d7cb5cd872eb1bad08280a437185b664d86e25d87b2217
    data/ligands/potency_wholecell.smi  sha256 94cbef9c27e729666f8c799ef75f6392fe61fa8db201d13778b76534f8694ba8
    ChEMBL cache (raw pull): work/chembl_potency/{enzyme,wholecell}.jsonl + MANIFEST.json

The `.smi` (SMILES + ChEMBL id) is what gets docked; the `.tsv` carries every compound's
median pChEMBL, record count, standard_types, heavy-atom count, triage verdict and (for
near-known) its nearest validated active. The four **near-known** compounds in each set are
near-duplicates (Tanimoto ≥ 0.70) of the 15 validated azoles — in enzyme: CHEMBL106
(fluconazole), CHEMBL638 (voriconazole), CHEMBL157101 + CHEMBL5555392 (ketoconazole).

## Primary hypothesis (H1) — the only thing that gates

> **H1.** Over the **enzyme** in-domain set (n=32), the Spearman rank correlation between
> closest N–Fe distance and median pChEMBL is **negative**: tighter iron approach ↔ higher
> potency.

**Statistic (fixed):** Spearman ρ (`openafr.potency.spearman`, Pearson of average ranks).
**Significance (fixed):** one-sided permutation p-value in the predicted direction
(`alternative="less"`, 100,000 permutations, seed 42). **Uncertainty (fixed):** 95%
percentile bootstrap CI (10,000 resamples of the (N–Fe, pChEMBL) pairs, seed 42).

> **PASS bar (pre-committed):** ρ ≤ **−0.30** AND one-sided permutation **p < 0.05**.

ρ ≤ −0.30 is a modest-but-real monotonic effect (the geometry explains ≈ a tenth of the rank
variance); the bar is set here, before any distance is read, so it can neither be lowered to
manufacture a pass nor raised to dismiss one. The exact ρ, its CI and p are reported **whatever
they are** — a pre-registered honest null is a valid result and is expected to be plausible
here (the geometry is known to read chemotype, not within-class potency:
`work/RESULTS_verified_inactives.md`, `work/TODOS.md`).

## Secondary analyses (reported, do NOT gate, cannot rescue a failed H1)

1. **Whole-cell** ρ/p/CI on the 441-compound set — well-powered, but confounded, so it is
   context, not proof of an on-target potency relationship. A passing whole-cell result with a
   failed enzyme H1 is reported as exactly that: consistent-with, not evidence-for, on-target
   potency tracking (the confound is un-removed).
2. **Leakage-guarded enzyme** ρ with the 4 near-known azoles **excluded** (n=28) — shows the
   effect is not carried by re-scoring near-duplicates of the training actives.
3. **Iron-bound-only enzyme** ρ restricted to compounds that actually coordinate the iron
   (closest N–Fe ≤ 3.0 Å) — tests whether, *among molecules that reach the iron*, how closely
   they reach tracks potency (the within-mechanism question).

## Stopping rule (pre-committed)

- If **H1 fails**, the settled conclusion is: *the geometry criterion enriches actives but does
  not demonstrably rank continuous on-target potency* — recorded as such, no further subset
  hunting on these poses. A passing secondary does **not** convert a failed H1 into a potency
  claim.
- If **H1 passes**, the claim upgrades to "geometry tracks on-target potency (weak-moderate,
  ρ≈X, CI …)" — never "predicts potency", and always carrying the limitations below.
- Either way this document, the datasets and the grader are the whole record; the result file
  is `work/RESULTS_potency.md`.

## Limitations (fixed in advance, carried into the result)

1. **Narrow on-target dynamic range.** The enzyme pChEMBL spans only 5.94–8.00 (~2 log units);
   a truncated range mechanically shrinks any correlation and widens its CI at n=32. This is a
   property of what CYP51-enzyme data ChEMBL has, not a choice.
2. **The applicability domain truncates the potent tail.** The ≤45-heavy-atom envelope excludes
   the largest, most potent azoles (posaconazole, itraconazole are themselves out-of-domain),
   so the very top of the potency axis is not on the plot. Honest, and pre-registered
   (`PREREGISTRATION_run2.md`), but it caps the achievable ρ.
3. **Whole-cell ≠ enzyme.** The secondary's potency is growth inhibition, confounded by
   permeability, efflux and target abundance — a compound can bind CYP51 and still shift its
   whole-cell pMIC. The confound depresses OR inflates; it cannot be signed away.
4. **Single rigid conformer, single seed.** The N–Fe distance is one 5TZ1 pose set at seed 42;
   the ensemble looks (#8/#9) show the rigid conformer is a hard ceiling for within-class
   ranking, and continuous potency is a within-class question. Expect this to bound the result.
5. **pChEMBL heterogeneity.** Pooling IC50/Ki/Kd/EC50 across assays and labs adds measurement
   noise to the dependent variable (median mitigates, does not remove it).
