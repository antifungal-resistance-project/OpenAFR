# Pre-registration — human-CYP51 positive control, redone as a native-pose-recovery test (issue #73)

**Written 2026-08-26, BEFORE reading any nitrogen-to-iron distance from this control run.**
Everything below is fixed in advance; deviations are recorded as deviations, not edited in.
Inherits the run-2 protocol and the consensus criterion — read
[PREREGISTRATION_run2.md](PREREGISTRATION_run2.md),
[PREREGISTRATION_consensus.md](PREREGISTRATION_consensus.md) and
[PREREGISTRATION_human_selectivity.md](PREREGISTRATION_human_selectivity.md) first.

## Why this control exists (and why the first one did not fire)

The human-CYP51 selectivity counter-screen ([RESULTS_human_selectivity.md](RESULTS_human_selectivity.md),
#73) pre-committed to a single positive control: re-dock ketoconazole into the prepared human
receptor (`work/receptor_human.pdbqt`, from 3LD6) and require a coordinating pose (N-Fe < 3.0 Å)
near the crystallographic 2.41 Å. That control **did not fire** — over 5 seeds ketoconazole's
best human N-Fe was 3.675 Å, never < 3.0 Å.

The post-hoc diagnosis in that result was that the failure is a **ligand global-search** problem,
not a receptor problem: ketoconazole is 36 heavy atoms with 8 rotatable bonds and a long flexible
imidazole-bearing tail, and it *also* failed to coordinate on the already-validated **fungal**
receptor (1/5 seeds). Two independent facts said the human receptor+box were sound anyway — small
rigid azoles coordinate the human iron reproducibly (fluconazole 5/5, letrozole 5/5), and the iron
sits well inside the box.

That diagnosis was written **after** the number was read, so it cannot promote the #73 numbers on
its own. This control tests it as a fixed, pre-committed hypothesis, and it does so with a test
that isolates the one variable actually in question.

**The distinction this control draws (fixed in advance).** A global redock conflates two very
different failures: (a) the receptor/scoring function does not recognize the true coordinating pose
as favorable — a *receptor* fault that would invalidate every selectivity number; and (b) the
receptor is right but Vina's global Monte-Carlo search cannot *find* the coordinating pose for a
floppy ligand — a *search-power* fault that says nothing about the receptor and does not touch the
selectivity margins (which are read on small, rigid ligands that search finds easily). The original
control could not tell these apart. A native-pose-recovery (local-optimization) control can: it
asks only whether the **crystallographic** ketoconazole pose is a coordinating, favorably-scored
minimum in this receptor — exactly question (a), decoupled from search power.

## Fixed inputs

- **Receptor:** `work/receptor_human.pdbqt` and iron anchor `work/receptor_human_A.pdb`,
  unchanged from #73 (built by `scripts/prep_receptor_human.py` from 3LD6). Iron position read
  from the anchor is (36.841, 3.766, -5.549); crystallographic ketoconazole N2 is 2.405 Å from
  it (verified from the deposited coordinates, matching 3LD6's LINK record of 2.41 Å).
- **Ligand (the crystallographic conformer, NOT a fresh embed):** the bound ketoconazole (residue
  `KKK`, chain A, 36 atoms) extracted verbatim from `data/structures/3LD6.pdb`, then converted to
  PDBQT with **the identical ligand-prep obabel step used by every docked ligand in this project**
  (`obabel … -p 7.4`, AutoDock atom types). The ONLY difference from the #73 ketoconazole input is
  the source of the heavy-atom coordinates: the crystal pose here vs. an RDKit ETKDG embed there.
  Protonation (obabel -p 7.4) is held identical, so the isolated variable is conformer/pose, not
  charge state.
- **Box:** center from `work/receptor_human_box.json` (42.288, 4.970, 1.219), size / num_modes
  from the frozen protocol (26 Å, 20) — unchanged from #73.

## The two runs (both fixed here)

1. **Native-pose recovery — the decisive control (`vina --local_only`).** Local-optimize the
   crystallographic ketoconazole pose in the human receptor. Search power is removed by
   construction; this measures only whether the true pose is a coordinating, negative-affinity
   minimum. This is the pre-committed test of receptor/scoring correctness.
2. **Global redock from the crystal conformer (context, not the decider).** Re-run the frozen
   global protocol (exhaustiveness 32, num_modes 20, seeds {42, 1, 2, 3, 4}) but starting from the
   crystal conformer's fragment geometry. Reported to show whether a better starting conformer,
   alone, lets global search recover the pose. Because Vina randomizes torsions at the start of
   every Monte-Carlo run, this is still fundamentally a search-power test and **may still fail;
   that outcome would not change the verdict below** — it would only re-confirm that the #73
   failure was search power, exactly the hypothesis.

## Pass criterion (frozen, decides the promotion)

The receptor is validated as geometrically correct — and the #73 selectivity numbers are promoted
from *provisional* to *validated* — **iff run 1 (`--local_only`) yields a coordinating pose with
N-Fe < 3.0 Å (the established `FE_CUTOFF`) and a negative Vina affinity.** As a fidelity read
(reported, not gating), the recovered coordinating nitrogen's N-Fe is expected near the crystal
2.405 Å.

If run 1 does **not** coordinate, the receptor itself is suspect and every #73 selectivity number
stays provisional (and the human counter-screen would need re-building, not re-interpreting).

Run 2's outcome is reported for completeness and does not, by itself, pass or fail this control.

## What a PASS does and does not license (fixed in advance)

- **Does:** confirm the human receptor+box+scoring recognize the true ketoconazole–iron
  coordination, removing the one unmet precondition behind #73, so the #73 margins (vatalanib,
  flucloxacillin, taranabant, lersivirine cleanly fungal-selective; nolatrexed a human liability;
  etc.) may be reported as validated rather than provisional.
- **Does not:** change any selectivity margin (none is re-docked here), nor claim ketoconazole is
  findable by global search, nor extend beyond CYP51 to other human CYPs (issue #74). The
  single-conformer, structural-necessary-not-sufficient limitations of #73 all still bind.

## Method provenance

Run by `scripts/human_positive_control.py` (deterministic given the frozen inputs above); it
extracts KKK, prepares the ligand, invokes Vina for both runs, and measures N-Fe with the repo's
canonical `openafr/pdbqt` parsers. Result: `work/RESULTS_human_control.md`.
