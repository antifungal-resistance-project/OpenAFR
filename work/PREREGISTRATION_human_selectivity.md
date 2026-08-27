# Pre-registration — human-CYP51 selectivity counter-screen (issue #73)

**Written 2026-08-26, BEFORE reading any nitrogen-to-iron distance from the human dockings.**
Everything below is fixed in advance; deviations are recorded as deviations, not edited in.
Inherits the run-2 protocol and the consensus criterion —
read [PREREGISTRATION_run2.md](PREREGISTRATION_run2.md) and
[PREREGISTRATION_consensus.md](PREREGISTRATION_consensus.md) first.

## Background — why this run exists

The method scores how closely a nitrogen reaches the **fungal** heme iron. But humans have
their own CYP51 (lanosterol 14α-demethylase), and a molecule tuned to coordinate the fungal
heme iron can just as easily coordinate the human one — a core antifungal toxicity mechanism
(hepatotoxicity, drug–drug interactions). A candidate is only interesting if it reaches the
**fungal** iron more readily than the **human** iron. Nothing in the repo checked this before
issue #73. This run docks the shortlist against human CYP51 under the identical criterion and
protocol, so the only difference between a fungal and a human score is the protein.

## Fixed inputs

- **Ligand set:** identical to the seed-stability run — `work/candidates/shortlist_focus.smi`
  (10 novel survivors + reference azoles + ketoconazole). Same prepared 3D ligands.
- **Human receptor:** `work/receptor_human.pdbqt`, built by `scripts/prep_receptor_human.py`
  from RCSB **3LD6** (human CYP51 in complex with ketoconazole), chain A protein + chain-A
  heme, prepared byte-for-byte the same way as the fungal receptor (obabel -xr -p 7.4). Anchor
  `work/receptor_human_A.pdb`.
- **Docking box:** centered on the crystallographic ketoconazole (KKK) centroid — the human
  active site — derived from the structure into `work/receptor_human_box.json`
  (center 42.288, 4.970, 1.219). Box **size, exhaustiveness, num_modes and seed list are
  unchanged** from the frozen protocol; only the box CENTER moves to the human pocket, because
  a different receptor must be scored on its own active site. This is the exact analogue of the
  fungal box being centered on 5TZ1's bound VT1.
- **Seeds:** {42, 1, 2, 3, 4}, so the human score is a 5-seed consensus min, apples-to-apples
  with the fungal consensus.

## Positive control (validates the human receptor + box before any candidate is trusted)

3LD6's own LINK record ties the ketoconazole N2 to the heme iron at **2.41 Å**. Re-docking
ketoconazole into this prepared receptor must reproduce a coordinating pose (N-Fe < 3.0 Å) near
that distance. If it does not, the human receptor/box is wrong and **no selectivity number from
this run is valid** — this is a hard precondition, fixed here in advance.

## The measure (fixed)

Per molecule, consensus min N-Fe on the human receptor `human_consensus` (same criterion,
FE_CUTOFF = 3.0 Å). Then:

    selectivity_margin = human_consensus − fungal_consensus   (Angstrom)

**Sign convention (frozen):** a **positive** margin means the nitrogen reaches the FUNGAL iron
closer than the human iron = fungal-selective = desirable. `reaches_human_iron` =
(human_consensus < 3.0 Å): a candidate that coordinates the fungal iron but **not** the human
iron is cleanly selective; one that coordinates both is a selectivity concern to be read
against the reference azoles.

Computed by `scripts/shortlist_confidence.py`.

## Pre-committed interpretation

- The reference azoles (fluconazole, voriconazole, ketoconazole, …) set the baseline: real
  clinical azoles are expected to reach the human iron too (that is their known toxicity), so
  their margins define "no better than an existing drug." A candidate is only a selectivity
  *win* if it is cleaner than that baseline (does not reach the human iron, or has a clearly
  larger positive margin than the reference azoles).
- A candidate that reaches the human iron **as well as or better than** the fungal iron
  (margin ≤ 0) is **flagged as a selectivity liability**, regardless of its fungal rank.
- This is a **counter-screen, not a gate on fungal activity**: it demotes/annotates, it never
  promotes. A clean human margin does not make a molecule active; it only removes one specific
  toxicity red flag.

## Limitations (fixed in advance)

1. Rigid single-conformer human receptor (one crystal form, 3LD6); the human pocket is also
   flexible. Same single-conformer caveat as the fungal side.
2. Docking-geometry selectivity is a **structural** signal only — it says nothing about the
   other human CYPs (3A4/2C9/…, issue #74), permeability, or metabolism.
3. Coordination is necessary-not-sufficient on the human side too: a short human N-Fe indicates
   a *possible* off-target interaction, not a measured one.
4. Chain A of 3LD6 is used, matching the fungal chain-A choice; the B copy is not scored.
