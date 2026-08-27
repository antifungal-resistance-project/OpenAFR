# Pre-registration — candidate resistance-pocket retention (issues #69, #77)

**Written 2026-08-27, BEFORE reading any nitrogen-to-iron distance from the *C. auris*
Y132F candidate dockings.** Everything below is fixed in advance; deviations are recorded
as deviations, not edited in. Inherits the run-2 protocol and the consensus criterion —
read [PREREGISTRATION_run2.md](PREREGISTRATION_run2.md),
[PREREGISTRATION_consensus.md](PREREGISTRATION_consensus.md) and
[PREREGISTRATION_auris_Y132F.md](PREREGISTRATION_auris_Y132F.md) first.

## Background — why this run exists

The consolidated scorecard ([RESULTS_scorecard.md](RESULTS_scorecard.md)) ranks and vets
its 10 novel candidates entirely against the **wild-type *C. albicans* 5TZ1 stand-in**. But
the mission is the **resistant *C. auris*** pocket. The scorecard itself prints
"auris ERG11 + Y132F/K143R mutant retention (#69, #77)" as `not-yet-checked` — it is the
first mission-shaped hole on every card. The held-out-azole transfer test already ran on the
Y132F pocket ([RESULTS_auris_Y132F.md](RESULTS_auris_Y132F.md)); this run applies the *same
receptor and protocol* to the **candidate path** — the molecules we would actually forward.

This is a **per-candidate confidence axis, not a gate.** It does not justify or block a new
screen. It adds one honest column per candidate to the scorecard: does the iron-coordinating
geometry that placed a molecule on the shortlist survive the dominant resistance substitution?

## Fixed inputs

- **Ligand set:** identical to the seed-stability and human-selectivity runs —
  `work/candidates/shortlist_focus.smi` (10 novel candidates + 6 known-azole controls). Same
  prepared 3D ligands (`work/candidates/cand_ligands/`).
- **Receptor:** `work/receptor_auris_Y132F.pdbqt`, prepared by
  `scripts/prep_receptor.py --mutant work/receptor_auris_Y132F.pdb` (sha256 368585eae998…)
  from the Y132F anchor `work/receptor_auris_Y132F.pdb` (sha256 924cf7b8…, the same file the
  held-out transfer test used, built by `scripts/mutate_receptor.py Y132F` — a pure para-
  hydroxyl deletion, heme iron unchanged). Same coordinate frame as `work/receptor_A.pdb`.
- **Docking box + search:** the frozen `protocol.yaml` **unchanged** — center 70.61/66.28/4.18,
  size 26, exhaustiveness 32, num_modes 20. The Y132F receptor shares 5TZ1's frame, so the box
  does **not** move (unlike the human counter-screen, where the pocket is elsewhere). The
  **only** variable versus the WT candidate dock is the residue-132 side chain.
- **Seeds:** {42, 1, 2, 3, 4} — a 5-seed consensus min, apples-to-apples with the WT fungal
  consensus already in `work/candidates/shortlist_confidence.json`.

## Wild-type baseline (no re-dock)

The WT candidate baseline is the **already-frozen** `per_seed_fungal` consensus in
`work/candidates/shortlist_confidence.json` — the same 16 ligands, the same 5 seeds, the same
frozen protocol, receptor `work/receptor_A.pdb`. It is **not** re-docked here; re-docking the
WT would change the baseline and break the controlled comparison. Retention is measured as the
shift of the Y132F consensus relative to that fixed WT consensus.

## The measure (fixed)

Per molecule, `y132f_consensus` = MIN nitrogen-to-iron distance over the seeds that produced a
coordinating pose (< FE_CUTOFF = 3.0 Å), with `n_coord/n_seeds` and `spread` — the **exact**
`consensus()` logic in `scripts/shortlist_confidence.py`, unchanged. Then:

    retention_shift = y132f_consensus − fungal_consensus   (Angstrom; + = moved AWAY from iron)

## Pre-committed per-candidate classification

Frozen here, before any Y132F distance is read. `RETENTION_TOL = 0.5 Å` — chosen as
"materially unchanged": half the `SPREAD_MAX = 1.0 Å` stable-consensus tolerance, and larger
than the <0.1–0.3 Å top-rank shifts the full-set Y132F run reported, so a candidate flagged
SHIFTED has moved by more than plausible run-to-run noise.

- **RETAINED** — coordinates the Y132F iron in **≥ 3/5 seeds** AND `|retention_shift| ≤ 0.5 Å`.
- **SHIFTED** — coordinates in **≥ 3/5 seeds** but `|retention_shift| > 0.5 Å` (direction
  reported; a molecule that moves *closer* is reported as SHIFTED-closer, not silently praised).
- **LOST** — coordinates in **< 3/5 seeds** on the Y132F pocket: it reliably fails to reach
  the resistant iron. This is the informative demotion this run exists to surface.

## Pre-committed interpretation (read the addendum lesson — do NOT overclaim)

- **Retention ≠ resistance-breaking.** The held-out Y132F run showed the geometric criterion
  is *nearly insensitive* to the Y132F change (AUC 0.794 → 0.805), so broad retention is close
  to a **tautology**, not evidence a molecule *defeats* resistance. A **RETAINED** verdict
  means only: the coordination geometry that put this molecule on the shortlist is **not
  destroyed** by the dominant resistance substitution — a robustness check.
- **The LOST verdict is the more informative outcome.** It flags a shortlist molecule whose
  geometric basis collapses in the pocket the mission actually targets — a real demotion.
- This axis **demotes/annotates, never promotes.** It is added to the scorecard as one column
  with a `not-yet-checked` → measured state; it does not collapse into any composite score and
  does not make any molecule a hit. Genuine resistance-breaker discovery needs actives
  *measured against the resistant strain* (the T1 active-set expansion), which this is not.

## Limitations (fixed in advance)

1. **A rigid single-atom deletion under-models real resistance** (same caveat as the held-out
   Y132F run): clinical Y132F also involves conformational change, H-bond rearrangement, and
   co-occurring efflux/expression changes. This captures pocket geometry only.
2. **Y132F only.** K143R (clade I) and F126L (clade III) need a side-chain repacker not in the
   pinned environment; they remain `not-yet-checked`.
3. **Single conformer, rigid receptor**, coordination necessary-not-sufficient — every caveat
   of the underlying screen still binds.
4. **Controls are controls.** The 6 known azoles are Y132F-defeated clinically; ranking them
   well on the mutant pocket is a robustness demonstration, never a resistance claim.
