# Docking all 8 known azoles + VT-1161 control into 5TZ1 (C. albicans CYP51)

Date: 2026-07-30
Receptor: 5TZ1 chain A + heme, rigid (work/receptor.pdbqt)
Engine: AutoDock Vina 1.2.7, seed 42 throughout (reproducible)
Constraint: an azole nitrogen must sit within 3.0 A of the heme iron

## Headline

**8 of 9 ligands can achieve iron coordination — but NOT under a single fixed protocol.**
Three different parameter settings were needed. One ligand fails outright.

| ligand | heavy atoms | rot. bonds | protocol that worked | best N-Fe | outcome |
|---|---|---|---|---|---|
| VT1161 (control) | 37 | — | exh 16, box 26 | 2.64 A | iron-bound |
| clotrimazole | 25 | 4 | exh 16, box 26 | 2.47 A | iron-bound |
| voriconazole | 25 | 7 | exh 16, box 26 | 2.67 A | iron-bound |
| miconazole | 25 | 6 | exh 16, box 26 | 2.86 A | iron-bound |
| isavuconazole | 31 | 8 | exh 16, box 26 | 2.72 A | iron-bound |
| ketoconazole | 36 | 8 | exh 16, box 26 | 2.87 A | iron-bound |
| fluconazole | 22 | 6 | **exh 64** needed | 2.67 A | iron-bound (sampling-limited) |
| itraconazole | 49 | 13 | **box 38, recentered** needed | 2.88 A | iron-bound (box-limited) |
| posaconazole | 51 | 15 | none tried worked | 3.96 A | **FAILS** |

## Finding 1 — Vina never ranks the correct pose first

Across every ligand that produced an iron-bound pose, that pose was **never** Vina's
top-ranked result. The constrained score is always worse than the unconstrained top score.
This reproduces the single-ligand VT-1161 finding across 8 independent ligands.

**Implication:** ranking on docking score alone is systematically wrong on this target.
The metal-coordination constraint is not a nice-to-have; it is load-bearing.

## Finding 2 — one fixed protocol produces systematic false negatives

Under the initial fixed protocol (exhaustiveness 16, box 26 A) only 6/9 succeeded.
The three failures split by cause:
- fluconazole: **sampling** — 4x more search found the pose (small ligand, just unlucky).
- itraconazole: **box geometry** — its ~25 A length cannot reach the iron and still fit
  inside a 26 A box (max room from the iron in any direction was under 20 A). A larger,
  recentered box fixed it.
- posaconazole: **unresolved** — still fails at exh 64 and box 38. Most likely induced fit:
  5TZ1 was solved with the short ligand VT-1161, so the access channel is not modeled in a
  conformation that accommodates posaconazole's long tail. A rigid receptor cannot open it.

**Implication for the screen:** running 10,000 library molecules through ONE fixed protocol
will systematically fail large/flexible molecules. Enrichment measured that way is biased
against exactly the chemotypes that include two of the most potent real drugs. The screen
protocol must be validated per size/flexibility class, or the bias must be reported.

## Finding 3 — docking score is mostly a size measurement (r = -0.91)

Correlation between heavy-atom count and raw Vina score across the 9 ligands: **r = -0.91**.
Bigger molecules score better almost automatically.

Ranking by raw score vs by ligand efficiency (score per heavy atom) inverts the list:
posaconazole and itraconazole are 1st/2nd by raw score and last by efficiency.

**Implication for decoys:** this is the circular-enrichment trap made concrete. If decoys
are not property-matched on size, the pipeline will "enrich" by preferring heavy molecules
and the validation gate will pass for the wrong reason. Property-matched decoys are
mandatory, not optional. Consider ranking on ligand efficiency rather than raw score.

## Design decisions this forces

1. Validation gate ranks on **constraint + score**, never score alone. (confirmed x8)
2. Decoys **must** be property-matched on heavy-atom count / MW. (r = -0.91)
3. The screen protocol must be **validated per ligand class**, and any class that cannot
   be docked reliably (posaconazole-like) must be reported as a known blind spot, not
   silently scored.
4. Consider ligand efficiency as the ranking metric instead of raw affinity.
5. Long-tailed azoles may require a flexible receptor or ensemble docking — out of scope
   for the spike, but it is a named limitation.

## Reproduce

    conda run -n openafr --no-capture-output python scripts/prep_ligands.py \
      data/ligands/actives.smi work/ligands
    # dock_batch.sh was the training-run spike script (box 70.61,66.28,4.18 size 26
    # exhaust 16), removed when docking was consolidated into scripts/screen.sh.
    # Those exact constants — not the current frozen protocol — define this run.
    conda run -n openafr --no-capture-output python scripts/analyze_poses.py \
      work/docked work/receptor_A.pdb
