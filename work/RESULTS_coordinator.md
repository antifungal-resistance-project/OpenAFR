# Results — which nitrogen coordinates the iron (coordinator identity)

Date: 2026-09-02. Pre-registered in
[PREREGISTRATION_coordinator.md](PREREGISTRATION_coordinator.md) (sha256 `b595de48…`,
`work/PREREG_coordinator.sha256`), frozen before the full-shortlist run. Protocol hash
unchanged. Pure classifier in [`openafr/coordinator.py`](../openafr/coordinator.py)
(7 known-answer tests, [`tests/test_coordinator.py`](../tests/test_coordinator.py));
runner [`scripts/coordinator_identity.py`](../scripts/coordinator_identity.py).

> **A per-candidate confidence axis, not a gate. It annotates and demotes; it never
> promotes, never re-runs the gate, never claims activity.** It answers one question the
> validated criterion leaves open: the rank measures the closest approach of **any**
> nitrogen to the heme iron — but *which* nitrogen? The method was validated on the azole
> mechanism (an aromatic ring N donating to the iron); a nitrile / amine / azide nitrogen
> that posts the same distance is coordinating through chemistry the method never validated.

## Why this run exists — and why verinurad triggered it

This axis was built during the "harden the single Priority-A candidate" pass. The dossier
(#119) put exactly one molecule at Priority A — **verinurad** — flagged only
"coordinator-fragile" (#87: a single aromatic N). Reading the criterion
(`openafr.pdbqt.min_nitrogen_iron_distance`) exposed a sharper question than the flag posed:
it takes the nearest nitrogen of *any* element-N type, and verinurad has **two** nitrogens —
a pyridine ring N (aromatic, the assumed azole-mimetic coordinator) and a **nitrile** N
(`C#N`, terminal, sp). Which one sets its rank?

## The classifier validates on every control

The 6 known azoles coordinate through a triazole/imidazole ring nitrogen by construction, so
the classifier **must** call them `aromatic-ring`. It does — **6/6 controls in-mechanism**
(fluconazole, voriconazole, ravuconazole, letrozole, flutrimazole, ketoconazole). The purely
geometric rule (a nitrile is the only nitrogen with a single sub-1.25 Å triple bond; a ring N
carries two ~1.34 Å bonds) recovers the known azole mechanism on every positive control. That
is the licence to trust its calls on the novel candidates.

## The result (5-seed consensus, WT *C. albicans* 5TZ1)

`report N-Fe` = the reported rank distance (frozen seed-42 min, the value the dossier ranks
on). `coord kind` = the nitrogen achieving it. `aromatic N-Fe` = the closest an **aromatic
ring** N gets, over all poses × 5 seeds — the in-mechanism distance behind the rank. Band =
the validated actives' 2.47–2.88 Å.

| candidate | report N–Fe | coord kind | aromatic N–Fe | verdict |
|---|---:|:--|---:|:--|
| flutrimazole\* | 2.371 | aromatic-ring | 2.371 | in-mechanism |
| **verinurad** | 2.417 | **nitrile** | 2.790 | **dual-mode** |
| vatalanib | 2.479 | aromatic-ring | 2.479 | in-mechanism |
| flucloxacillin | 2.558 | aromatic-ring | 2.517 | in-mechanism |
| **taranabant** | 2.586 | **nitrile** | 6.870 | **out-of-mechanism** |
| **lersivirine** | 2.599 | **nitrile** | 4.052 | **out-of-mechanism** |
| LDN-27219 | 2.603 | amine | 3.466 | out-of-mechanism |
| R-1479 | 2.621 | azide (sp)† | 3.038 | out-of-mechanism |
| ravuconazole\* | 2.653 | aromatic-ring | 2.653 | in-mechanism |
| L-838417 | 2.667 | aromatic-ring | 2.580 | in-mechanism |
| olprinone | 2.680 | nitrile | 3.697 | out-of-mechanism |
| nolatrexed | 2.692 | aromatic-ring | 2.677 | in-mechanism |
| letrozole\* | 2.719 | aromatic-ring | 2.705 | in-mechanism |
| fluconazole\* | 2.740 | aromatic-ring | 2.711 | in-mechanism |
| voriconazole\* | 2.799 | aromatic-ring | 2.657 | in-mechanism |
| ketoconazole\* | 2.921 | aromatic-ring | 2.921 | in-mechanism |

\* known azole control (validates the classifier; exempt from demotion).
† R-1479's short-triple-bond coordinator is its **azide** terminus, not a nitrile — the same
`1 neighbour, bond < 1.25 Å` bucket. Either way it is not an aromatic ring N: out-of-mechanism.
(R-1479 was already Priority-C for PAINS azide/diazo, #75/#119.)

## What this changes — the finding

**The single Priority-A candidate is downgraded, and two of the four "strongest" hits are
out-of-mechanism artifacts.**

- **verinurad — dual-mode; demote from Priority-A.** Its reported 2.417 Å (rank #2, the sole
  automatic-A) is set by its **nitrile** N in all 5 seeds. Its aromatic pyridine N — the
  coordinator the applicability domain and the "coordinator-fragile" flag reasoned about —
  reaches only **2.790 Å**, at the far edge of the validated band. Re-ranked on the
  in-mechanism distance, verinurad falls from #2 to ≈14th, next to the reference azoles. An
  azole-like binding mode is *available* (hence dual-mode, not out-of-mechanism), but the
  geometry that made it the top pick is a nitrile interaction the method never validated.
  **The automatic-A rested on an out-of-mechanism nitrogen.**
- **taranabant — out-of-mechanism.** Reported 2.586 Å is entirely its nitrile; its lone
  aromatic (pyridine) N never comes within **6.9 Å** of the iron. Its rank carries no
  azole-like coordination at all. (It was one of the four "strongest, cleanly human-selective"
  hits, #73.)
- **lersivirine — out-of-mechanism.** Reported 2.599 Å is one of its **two** nitriles; its
  pyrazole ring N stays **4.06 Å** away. Also one of the four "strongest" (#73), and the
  pose-convergence standout (#71, coord rank 1) — but the pose it converged on coordinates
  through a nitrile, not the ring. Both prior wins were reading the same out-of-mechanism pose.
- **LDN-27219 / olprinone** were already demoted (single-dock instability, #68); this axis
  adds the mechanism reason (amine / nitrile coordinator).

**The in-mechanism survivors** — reported rank *is* the aromatic ring N: **vatalanib**
(2.479), **flucloxacillin** (2.517), **L-838417** (2.580), **nolatrexed** (2.677). Of these,
vatalanib is also seed-stable (#68) and cleanly human-selective (#73); flucloxacillin carries
other flags (structural alert, conformer-fragility #84); nolatrexed hits human CYP51 (#73);
L-838417 is stable-provisional. **vatalanib is the one candidate clean on this axis *and* the
prior confidence axes.**

## The deeper cause — a criterion-level caveat this surfaces

`min_nitrogen_iron_distance` measuring the nearest nitrogen of any type is correct for the
azole training/validation set (imidazole/triazole drugs whose only pocket-facing N is the
coordinating ring N), but on a **novel** library it can silently reward nitrile-, azide-, and
amine-bearing molecules for pointing a non-aromatic N at the iron. The applicability-domain
`≥1 aromatic N` gate (#87) admits these molecules on the presence of an aromatic N, then the
criterion ranks them on a *different* nitrogen. This axis is the check that catches the gap.

**The obvious fix — restrict the criterion to the aromatic ring N — was pre-registered and
tested, and it FAILED** ([RESULTS_aromatic_criterion.md](RESULTS_aromatic_criterion.md)):
re-grading the run-2 gate under the aromatic-only criterion *lowers* the validated AUC from
0.794 to 0.729, because the active **luliconazole** itself coordinates the iron via its nitrile
in-pose (aromatic N 5.3 Å). So a non-aromatic coordinator is **unvalidated, not disproven**, and
the frozen criterion is left unchanged. The coordinator identity is therefore used as a
per-candidate **flag** (block automatic-A, surface prominently), never a global re-rank and
never a hard demotion — exactly how the dossier now applies it
([RESULTS_dossier.md](RESULTS_dossier.md)).

## What this establishes / does not

- **Establishes:** for every shortlist molecule, whether the reported iron-coordinating
  nitrogen is the validated aromatic mechanism or an out-of-mechanism nitrile/amine/azide;
  recovers the correct aromatic call on all 6 azole controls; downgrades verinurad from the
  sole Priority-A to dual-mode; and reclassifies taranabant and lersivirine — two of the four
  "strongest" hits — as out-of-mechanism.
- **Does not:** re-run the gate, re-derive AUC, promote any molecule, edit the frozen protocol,
  or relax any standing caveat (single rigid *C. albicans* conformer, coordination
  necessary-not-sufficient, per-dock variance). A `dual-mode` or `in-mechanism` verdict is a
  cleaner *hypothesis*, never a hit.

## Reproduce

    conda run -n openafr python scripts/prep_receptor.py
    conda run -n openafr python scripts/prep_ligands.py \
        work/candidates/shortlist_focus.smi work/candidates/coord_ligands
    conda run -n openafr bash scripts/screen_consensus.sh \
        work/candidates/coord_ligands work/candidates/dock_coordinator "42 1 2 3 4"
    conda run -n openafr python scripts/coordinator_identity.py \
        --docks work/candidates/dock_coordinator \
        --ligands work/candidates/shortlist_focus.smi --json \
        > work/candidates/coordinator_identity.json

Per-seed dock trees are gitignored; the classifier, runner, tests and the JSON summary persist.
