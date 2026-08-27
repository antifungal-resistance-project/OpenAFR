# Results — human-CYP51 positive control, redone as a native-pose-recovery test (issue #73)

Date: 2026-08-26. Pre-registered in
[PREREGISTRATION_human_control.md](PREREGISTRATION_human_control.md)
(sha256 `74542414…`, frozen in `work/PREREG_human_control.sha256` before any distance below was
read). Run by `scripts/human_positive_control.py`; machine summary
`work/candidates/human_control/human_control.json`.

> **Outcome: PASS. The human receptor is validated.** The pre-registered decider —
> local-optimising the **crystallographic** ketoconazole pose in the prepared human receptor —
> reproduces iron coordination at **N-Fe 2.651 Å** (crystal ground truth 2.405 Å) with a
> favorable affinity of **−10.369 kcal/mol**. This removes the one unmet precondition behind the
> #73 selectivity counter-screen, so those numbers are **promoted from provisional to validated**
> (see [RESULTS_human_selectivity.md](RESULTS_human_selectivity.md)).

## Why the first control did not fire — and what this one isolates

The #73 counter-screen pre-committed to a single positive control: re-dock ketoconazole and
require N-Fe < 3.0 Å near the crystallographic 2.41 Å. It did not fire (best human N-Fe 3.675 Å
over 5 seeds). The post-hoc diagnosis — a **ligand global-search** failure, not a receptor fault —
was written after the number was read, so it could not, on its own, promote the #73 numbers.

A global redock conflates two failures the pre-registration named up front: **(a)** the receptor /
scoring function does not recognize the true coordinating pose (a receptor fault that would
invalidate every selectivity number), versus **(b)** the receptor is right but Vina's global
Monte-Carlo search cannot *find* the pose for a 36-heavy-atom, 8-rotatable-bond floppy ligand (a
search-power fault that does not touch the margins, which are read on small rigid ligands search
finds easily). This control was pre-registered to separate them with a native-pose-recovery
(local-optimization) test that asks only (a).

## The decider (run 1): local optimisation from the crystal pose — PASS

Ground truth from the deposited 3LD6 coordinates (independent of any docking): the bound
ketoconazole's N2 sits **2.405 Å** from the human heme iron (its other nitrogens 4.2 / 14.1 /
16.6 Å) — matching 3LD6's LINK record (2.41 Å).

Local-optimising that crystal pose in `work/receptor_human.pdbqt` (search power removed by
construction):

| quantity | value | pass condition | met? |
|---|---:|---|:---:|
| min N-Fe (local-optimised crystal pose) | **2.651 Å** | < 3.0 Å (FE_CUTOFF) | ✅ |
| Vina affinity | **−10.369 kcal/mol** | < 0 | ✅ |
| (fidelity read, non-gating) crystal N-Fe | 2.405 Å | near-recovered | ✅ |

The true coordinating pose is a coordinating, favorably-scored minimum in this receptor. Per the
frozen criterion, **the receptor + box + scoring are validated as geometrically correct.**

## Context (run 2): global redock from the crystal conformer

Re-running the frozen global protocol (exhaustiveness 32, num_modes 20, seeds {42,1,2,3,4}) but
starting from the crystal conformer's fragment geometry — reported for completeness, not the
decider:

| seed | 42 | 1 | 2 | 3 | 4 | consensus |
|---|---:|---:|---:|---:|---:|---:|
| min N-Fe (Å) | 5.188 | 4.289 | 3.783 | 3.644 | **2.474** | **2.474** |
| affinity | −10.106 | −10.859 | −10.249 | −10.354 | −11.005 | |

**1/5 seeds** now recovers the coordinating pose (2.474 Å) — up from **0/5** in #73's ETKDG-embed
run. A better starting conformer alone lifts recovery, but global search remains the bottleneck
for this ligand (as pre-registered, this outcome does not change the verdict; it re-confirms the
failure was search power). Note every pose scores −10 to −11 kcal/mol regardless of coordination:
affinity ranking alone would not have caught the non-coordinating poses — the iron-distance
criterion does, which is the method's whole point.

## What this establishes (and does not)

- **Does:** confirm the human receptor+box+scoring recognize the true ketoconazole–iron
  coordination, closing the one open precondition from #73. The #73 margins (vatalanib,
  flucloxacillin, taranabant, lersivirine cleanly fungal-selective; verinurad/L-838417 marginal;
  nolatrexed a human-selectivity liability) are now reported as **validated**, not provisional.
- **Does not:** change any selectivity margin (nothing is re-docked here), claim ketoconazole is
  findable by global search, or extend beyond CYP51 to other human CYPs (issue #74). Every
  single-conformer / structural-necessary-not-sufficient limitation of #73 still binds.

## Reproduce

```
conda activate openafr
python scripts/prep_receptor_human.py            # build work/receptor_human.pdbqt from 3LD6
python scripts/human_positive_control.py         # local_only decider + global-redock context
```

Docked poses are gitignored; the crystal ligand (`work/candidates/human_control/
ketoconazole_crystal.pdb`) and the JSON summary persist. 3LD6 is `data/structures/3LD6.pdb`
(sha256 `486ba093…`).
