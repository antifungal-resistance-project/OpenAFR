# Pre-registration — WIDEN the active set against MEASURED inactives to retire the n=7 tip (Limitation #2)

> **STATUS: FROZEN 2026-09-08, before any AUC of the widened active set against the measured
> inactives was scored.** The active sample (`data/ligands/active_power_sample.smi`) and the
> measured-inactive set (`data/ligands/verified_inactives.smi`) were both frozen in earlier looks;
> their sha256s are pinned below and re-checked by the grader. This document is frozen by its own
> SHA-256. The grader (`scripts/validate_gate_verified_power.py`) re-checks this document's hash,
> `protocol.yaml`, and both input hashes, and refuses to certify a pass if any moved.

## Why this run is required — the "single largest caveat," still at n=7

Two limitations of the preprint (`work/PREPRINT_geometry_ceiling.md` §6) were each an *n=7*
tip. Limitation #1 (actives vs **property-matched, presumed-inactive** decoys) was retired at
power by look #10 (`work/RESULTS_active_power.md`): the identical gate with the active side
widened 7 → N=300 gave **AUC 0.750, 95% CI 0.717–0.782**, lower bound clearing the bar.

Limitation #2 — *"decoys are presumed inactive, not verified … the single largest caveat"* — was
attacked by look #6 (`work/RESULTS_verified_inactives.md`), which replaced the presumed decoys
with **279 compounds MEASURED not to inhibit *C. albicans***. It passed on the point estimate
(**AUC 0.716**) but on the **same 7 held-out actives**, so its bootstrap 95% CI was **0.506–0.885**
— the lower bound again sits **below** the 0.70 bar. So the strongest version of the caveat (real
azoles vs *measured* non-inhibitors) has never been powered.

This run does to Limitation #2 exactly what look #10 did to Limitation #1: it re-runs look #6's
contrast — same rigid 5TZ1, same mode C, same 279 measured inactives, same frozen protocol and
seed — with **the active side widened from 7 to the same N=300 sample look #10 used**, so the CI
becomes tight and its lower bound can be compared against the bar. One question: **once the
sampling noise of 7 actives is removed, does the geometry criterion separate real azoles from
compounds MEASURED not to work, with its AUC lower bound above the usefulness bar?**

This is a *power* run, not a new criterion and not a new dataset. Both the actives and the
inactives are pre-existing frozen files; only which actives are graded against the measured
inactives changes (7 → 300).

## This run is NOT blind — and says so

Look #6's subgroup AUCs (non-azole 0.810, azole-bearing 0.650) and look #10's secondary S1
(N=300 actives vs the **non-azole** subset, 0.774) are already known. So, exactly as look #10
framed itself, **the value here is the tight interval, not blindness.** The two numbers this run
measures for the first time are (a) the **primary**: N=300 actives vs **all 279** measured
inactives, and (b) the **crux secondary**: N=300 actives vs the **azole-bearing** measured
inactives (look #6 saw 0.650 there on only 7 actives). The pre-committed interpretation below is
fixed before either is read.

## Single-variable design — everything is held fixed from look #6 except the actives

| | |
|---|---|
| receptor | rigid **5TZ1** chain A + heme, single conformer — the product receptor, **unchanged** |
| criterion | **mode C** — `openafr.pdbqt.min_nitrogen_iron_distance` — **unchanged** |
| protocol | `protocol.yaml` unchanged — box 26 Å @ 70.61/66.28/4.18, exhaustiveness 32, 20 modes, seed 42 |
| **measured inactives** | the **same 279** compounds MEASURED not to inhibit *C. albicans* (`data/ligands/verified_inactives.smi`), **unchanged** from look #6 |
| **actives** | **CHANGED** — the N=300 blind measured-active azole sample (look #10's `data/ligands/active_power_sample.smi`), replacing look #6's 7 held-out azoles |

Only the actives change, so any difference from look #6's 0.716 is attributable to the active
sample size, not to method drift. The N=300 sample was built by the same frozen rule described in
`work/PREREGISTRATION_active_power.md` (seed-7 deterministic draw from the 3,609-molecule
`verified_actives.smi` population; ≤45 heavy atoms, ≥1 aromatic N, max Tanimoto <0.70 vs all 15
project actives + the 5TZ1 co-crystal), so every active is genuinely distinct chemistry (no
leakage) and in the declared applicability domain.

> **PINNED (grader re-checks all three before scoring):**
> - `data/ligands/active_power_sample.smi` sha256 `ba8d1f3cb05717017ababa321c8a74b4e0122c5eae94c91ceab92915dd3c37ae` (300 actives, seed 7)
> - `data/ligands/verified_inactives.smi` sha256 `37a20eec1c026c42d3c8efb53c1dd27cda2f9057ca541cfc1afed21a35faea6b` (279 measured inactives)
> - `protocol.yaml` sha256 `6eca16b1be2e5e33ae42d4b1b7a63a0e9b891a5dd6cbfb19642dd3e9a149aba4`

## Primary hypothesis (H1) — the only thing that gates

Ranking the pooled molecules (the N=300 measured-active azole sample + the 279 measured
inactives) by **mode C on rigid 5TZ1** separates measured-active azoles from measured
non-inhibitors at **AUC ≥ 0.70** (`protocol.yaml` `min_auc`, unchanged), AND the **bootstrap 95%
CI lower bound is reported**. The retirement of Limitation #2 is graded on the **lower bound**,
not merely the point estimate — the whole reason this run exists is that at n=7 the lower bound
(0.506) sat below the bar. Unranked actives (docking failures) are ranked **last**, never dropped
(the anti-inflation rule from Run 2 and look #6).

EF is **not** part of the gate (seed-fragile, preprint §3.2); reported for continuity only.

## Secondary — all pre-specified, none gating

1. **S1 — non-azole subgroup (novel-chemotype triage, the product claim).** N=300 actives vs the
   **non-azole** measured inactives. Comparators: look #6 0.810 (n=7), look #10 S1 0.774 (N=300,
   GCP). This run re-measures it on **freshly, locally docked** poses, so it also serves as a
   GCP→local **reproducibility check** of look #10 S1. Reported with bootstrap 95% CI.
2. **S2 — azole-bearing subgroup (the crux, the sharpest form of the caveat).** N=300 actives vs
   the **azole-bearing** measured inactives — can pose geometry tell a *working* azole from a
   *failed* azole analogue? Look #6 saw **0.650** here on 7 actives. Reported with bootstrap 95% CI.
3. **S3 — the tip + statistical support.** Best true active's percentile rank; how many measured
   inactives outrank it; composition of the validated iron-bound band (2.47–2.88 Å);
   permutation p (20,000 shuffles); bootstrap 95% AUC CI.
4. **S4 — domain-edge check.** AUC in heavy-atom halves of the actives, to confirm the widened
   sample is not carried by a molecular-size artifact.

## Pre-committed interpretation

**Primary (gates Limitation #2's retirement):**

- **RETIRED — H1 AUC ≥ 0.70 AND CI lower bound ≥ 0.70.** Limitation #2 is discharged at power:
  with the sampling noise of n=7 removed, the criterion separates real azoles from compounds
  *measured* not to inhibit, lower bound above the bar. The decoy ceiling is confirmed **not** an
  artifact of decoy construction (it was already measured; now it is powered). Preprint §6
  Limitation #2 is updated from "the single largest caveat" to the measured tight interval, and
  §4.1 gains the powered number.
- **MARGINAL — point ≥ 0.70 but CI lower bound < 0.70.** The gate passes on the point estimate but
  a sub-bar true AUC against measured inactives cannot be excluded even at N=300. Reported honestly
  as a narrowed-but-not-cleared interval; Limitation #2 is *quantified*, not retired.
- **FAIL — point estimate < 0.70.** Look #6's 0.716 was itself an n=7 fluctuation; broad
  enrichment against *measured* inactives does not hold at power. This is a material negative about
  the product's headline and is reported as prominently as the pass branches — the novel-chemotype
  triage scope would need re-statement, not a footnote.

**Secondary S2 (the within-warhead crux — reported, does NOT gate, interpretation fixed now):**

- **S2 AUC < 0.70 at N=300** → the within-azole-warhead ceiling is **confirmed real**, not an n=7
  fluctuation (consistent with look #9's independent closure of within-class ranking). This is the
  expected result and it *strengthens* the honest scope: the tool ranks chemotype, not potency.
- **S2 AUC ≥ 0.70 at N=300** → a surprising, well-powered reopening of the within-warhead question.
  It does **not** by itself change the product scope (novel-chemotype triage) and does **not**
  retire anything; per the stopping rule it triggers a *separate* pre-registered confirmation on an
  independent active draw before any within-class claim, exactly as look #8→#9 required.

## Stopping rule (binding)

**One sample, one look.** The actives are the frozen N=300 seed-7 sample; the inactives are the
frozen 279. Neither will be re-sampled/re-drawn and re-graded if a number disappoints. Receptor,
criterion, and protocol are unchanged. Any further work is a new, separately pre-registered dataset.

## Limitations (fixed in advance)

1. **Not blind** (see above) — the value is the interval, not concealment of the point estimate.
2. **Whole-cell "active" ≠ CYP51 "active"**, and **"inactive" = measured non-inhibitor of the
   whole cell** — both inherited from looks #6/#9/#10; each mislabel is conservative for the gate
   (a mislabelled active ranks low; a true binder among the inactives depresses AUC).
3. **Measured inactives are enriched for the warhead by construction** (58.4% azole-bearing) — the
   only compounds anyone measures against *C. albicans* are antifungal programs. This makes the
   primary a **harder** benchmark than look #10's property-matched decoys, by design; the drop from
   0.750 is expected and is the point of the S1/S2 split.
4. **Rigid single conformer / ≤45 heavy-atom domain** — unchanged; posaconazole/itraconazole
   long-tail remains the declared blind spot.
5. **A rigid-docking false negative ranks LAST, never dropped** — conservative for the gate.

## Reproduce (all local — Vina 1.2.7 / obabel 3.1.0 on Apple Silicon, the look #6 host class)

```
conda activate openafr
python scripts/prep_receptor.py                                             # -> work/receptor.pdbqt (from 5TZ1)
python scripts/prep_ligands.py data/ligands/active_power_sample.smi   work/verified_power_ligands
python scripts/prep_ligands.py data/ligands/verified_inactives.smi    work/verified_power_ligands
RECEPTOR=work/receptor.pdbqt scripts/screen.sh work/verified_power_ligands work/screen_verified_power
python scripts/validate_gate_verified_power.py                              # exits 0 pass / 1 fail
```
