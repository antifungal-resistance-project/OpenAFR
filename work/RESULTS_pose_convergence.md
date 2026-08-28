# Per-candidate pose-convergence diagnostics (#71)

Date: 2026-08-27. Pre-registered in
[PREREGISTRATION_pose_convergence.md](PREREGISTRATION_pose_convergence.md)
(sha `f65adbdb46ec…`, frozen BEFORE any pose-to-pose RMSD was read).
Grader: [`scripts/pose_convergence.py`](../scripts/pose_convergence.py); pure clustering in
[`openafr/pose_convergence.py`](../openafr/pose_convergence.py) (10 known-answer tests).
Data: [`work/candidates/pose_convergence.json`](candidates/pose_convergence.json).

## The question

The scorecard reports each candidate's **single closest** nitrogen-to-iron distance — the best
of the 20 poses one Vina search returns. Seed-stability (#68) asks whether that number
reproduces *across* independent searches. This axis asks the complementary *within*-search
question: of the 20 poses one search returns, do they **agree** on where the ligand sits, or
does the winning pose stand almost alone while the rest scatter across the pocket? A tight
N–Fe from one pose is weaker evidence if the other 19 land somewhere else.

**This is a per-candidate confidence axis, not a gate.** It demotes/annotates, never promotes.

## Method (frozen)

The one **frozen seed-42 dock** per candidate (protocol.yaml unchanged: 26 Å box,
exhaustiveness 32, num_modes 20), receptor `work/receptor_A.pdb` (WT *C. albicans* 5TZ1
stand-in). The 20 returned poses are clustered by **unaligned, index-matched heavy-atom RMSD**
at 2.0 Å (AutoDock convention), energy-ordered leader clustering; the top cluster is the
best-affinity pose's. **CONVERGED** iff ≥ 50 % of the poses fall in that cluster, else
**DISPERSED**. Reported alongside: the affinity rank of the coordinating (min-N–Fe) pose and
whether it lies inside the dominant cluster.

**Reproduction check passed:** the min N–Fe over the 20 poses equals the seed-stability run's
stored `per_seed_fungal["42"]` for **all 16/16** ligands (≤ 0.01 Å) — the fresh dock
reproduced the frozen protocol exactly.

## Result — the search does not converge under this protocol (every candidate DISPERSED)

| candidate | N–Fe (Å) | poses | top-cluster | n_clusters | coord rank | coord in top | verdict |
|---|---:|:--:|:--:|:--:|:--:|:--:|:--|
| flutrimazole* | 2.371 | 19 | 5 % | 15 | 8 | no | DISPERSED |
| verinurad | 2.417 | 20 | 5 % | 19 | 11 | no | DISPERSED |
| vatalanib | 2.479 | 19 | 10 % | 17 | 13 | no | DISPERSED |
| flucloxacillin | 2.558 | 19 | 5 % | 17 | 15 | no | DISPERSED |
| taranabant | 2.586 | 20 | 5 % | 17 | 12 | no | DISPERSED |
| **lersivirine** | 2.599 | 20 | 5 % | 20 | **1** | **yes** | DISPERSED |
| LDN-27219 | 2.603 | 20 | 10 % | 15 | 19 | no | DISPERSED |
| R-1479 | 2.621 | 20 | 10 % | 16 | 9 | no | DISPERSED |
| ravuconazole* | 2.653 | 20 | 10 % | 17 | 6 | no | DISPERSED |
| L-838417 | 2.667 | 20 | 10 % | 16 | 11 | no | DISPERSED |
| olprinone | 2.680 | 19 | 5 % | 14 | 15 | no | DISPERSED |
| nolatrexed | 2.692 | 20 | 5 % | 19 | 8 | no | DISPERSED |
| letrozole* | 2.719 | 20 | 5 % | 19 | 4 | no | DISPERSED |
| fluconazole* | 2.740 | 20 | 10 % | 14 | 11 | no | DISPERSED |
| voriconazole* | 2.799 | 20 | 5 % | 15 | 6 | no | DISPERSED |
| ketoconazole* | 2.921 | 17 | 18 % | 14 | 2 | yes | DISPERSED |

\* known-azole control.

**The headline is protocol-level, not per-candidate: no candidate's returned poses converge.**
The 20 modes fragment into 14–20 distinct 2-Å clusters — nearly one cluster per pose — and for
most candidates a large share of the modes dock **> 6 Å away** from the best pose, i.e. in a
spatially distinct region of the deliberately large 26 Å box (11/20 for verinurad, 12/19 for
vatalanib, 13/20 even for lersivirine). The reported tight N–Fe is therefore **one specific
mode among ~20 scattered ones**, not a basin the single search settled into. This is the
within-search face of the same reliability caveat that made the 5-seed consensus (#68)
necessary, and it applies to essentially the whole protocol — the 26 Å box (sized to "fit the
long azoles") plus num_modes 20 all but guarantees the returned modes spread out.

## What still discriminates between candidates

The CONVERGED/DISPERSED flag does not (all DISPERSED). The **coordinating-pose rank** does —
the affinity rank of the pose that gives the reported N–Fe:

- **lersivirine is the standout: coord rank 1, in the top cluster.** It is the **only novel
  candidate** whose reported coordinating geometry is *also* Vina's single best-energy pose —
  the two independent signals point at the same pose. (ketoconazole, a control, is the only
  other in-top case, at rank 2.) This is consistent with lersivirine being one of the two
  candidates clean on every other measured axis in the scorecard.
- For **most novel candidates the coordinating pose is a mid-to-high-energy mode**
  (flucloxacillin rank 15, LDN-27219 19, vatalanib 13, taranabant 12, verinurad/L-838417 11):
  the geometry that shortlisted them is *not* the pose Vina scores best.

**Read the rank honestly — a high coordinating rank is NOT a demotion.** The project's central
result is that Vina's affinity score is *untrustworthy* on this system (geometry AUC ~0.79 vs
Vina's own AUC 0.47). So a coordinating pose that is *not* Vina's best-energy pose is close to
the **expected** state — it is the very disagreement that motivates geometry rescoring in the
first place. What the rank measures is **agreement between the two signals**: when the
coordinating pose is *also* low-energy (lersivirine, rank 1), geometry and the docking score
concur, which is extra reassurance; when it is high-energy, they disagree, which is the norm
here and disqualifies nothing.

## What this establishes / does not

- **Does:** add one honest per-candidate column and surface a real protocol-level caveat — the
  single-dock poses do not converge, so every reported N–Fe is one scattered mode among ~20;
  and identify **lersivirine** as the one novel candidate where the coordinating geometry and
  the docking score agree on the same, best-energy pose.
- **Does not:** call any molecule a hit or promote one. Convergence is pose *reproducibility
  within one search*, **not** thermodynamic stability and **not** activity — coordination
  stays necessary-not-sufficient ([RESULTS_axial.md](RESULTS_axial.md)). DISPERSED is not
  disqualifying (it applies to the whole protocol); a high coordinating rank is not a demotion
  (the method distrusts the Vina score by design). The unaligned RMSD is conservative — it can
  call orientation-scattered poses DISPERSED even when their centroids nearly coincide.

## Reproduce

    conda run -n openafr python scripts/prep_receptor.py                     # -> work/receptor.pdbqt
    conda run -n openafr python scripts/prep_ligands.py \
        work/candidates/shortlist_focus.smi work/candidates/cand_ligands
    conda run -n openafr python scripts/pose_convergence.py \
        --ligands work/candidates/cand_ligands \
        --workdir work/candidates/dock_convergence --json > work/candidates/pose_convergence.json

Deterministic (fixed protocol seed). Docking needs vina + obabel from the pinned `openafr`
env; the clustering/grading is pure and offline.
