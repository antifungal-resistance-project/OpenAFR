# Pre-registration — the within-azole question under an ENSEMBLE receptor

**Status: FROZEN 2026-08-25, before any active/inactive label was scored.** This document
defines the one permitted next attack on the within-azole ceiling: a *real experimental
conformer ensemble* instead of the single rigid 5TZ1 pose. The three receptors were prepared
first (deterministic superposition, no labels involved) and their anchor hashes are pinned
below; this document is then frozen by its own SHA-256. The grader (`validate_gate_ensemble.py`)
re-checks this document's hash, the `protocol.yaml` hash, and the three receptor-anchor hashes
before it will certify anything, and refuses to certify a pass if any moved. No AUC on the
azole pool had been read at freeze time.

Pinned receptor anchors (produced by `scripts/prep_receptor_ensemble.py`, label-free):

    work/receptor_A.pdb       (5TZ1) sha256 = existing project anchor, unchanged
    work/receptor_5FSA_A.pdb  (5FSA) sha256 9531a1531234072dae1f35df89a8411479393a91e170ff8d8e2b56af6bb777e3
    work/receptor_5V5Z_A.pdb  (5V5Z) sha256 93a261010162782c625b36480568eb1b939b28414d9290817c2b4009e25103cd

Alignment quality (reported, not tuned): 5FSA superposed onto the 5TZ1 frame over 484 shared
CA at backbone RMSD 1.819 A; 5V5Z over 468 shared CA at 1.055 A. The three heme irons land
within 0.3 A of each other (5TZ1 64.163/71.316/2.711; 5FSA 64.230/71.622/3.008; 5V5Z
63.919/71.353/2.741), so the frozen 26 A box at (70.61, 66.28, 4.18) covers each pocket
identically.

## Why this look is permitted (and why it is the *only* permitted one)

Two pre-registered within-azole attacks have now failed on the **single rigid 5TZ1 conformer**:

    look #6  mode C (N–Fe distance)        AUC 0.650   (vs azole-bearing inactives)
    look #7  mode F (channel engagement)   AUC 0.590   (orthogonal to C, r = -0.301)

Both FAILs carried the *same* pre-committed escape clause, written before either result was
known (`work/PREREGISTRATION_channel_engagement.md`, Stopping rule; `work/TODOS.md`):

> "Any further within-class attempt requires a genuinely new dataset (an independent active
> set, or **ensemble/flexible-receptor poses**), separately pre-registered."

This document is exactly that — and it is disciplined in the way the stopping rule demands:
**new docking, a new dataset, one aggregation rule, and the *already-validated* criterion
re-used unchanged.** No new feature is invented on the old rigid poses. The single new
ingredient is the receptor: three conformations instead of one. This is the cleanest possible
test of the hypothesis that both prior FAILs pre-committed as their leading explanation — that
the ceiling is a property of the *rigid single conformer*, not of the geometric criterion.

## The ensemble is EXHAUSTIVE, not selected (this is the key anti-fishing property)

An ensemble introduces an obvious new researcher degree of freedom — *which* conformers, and
*how many*. This design removes it entirely: it uses **every** experimental *C. albicans* CYP51
(Erg11) structure in the PDB. There is no subset to choose, so there is nothing to tune.

The complete population, confirmed by RCSB query on 2026-08-25 (organism = *Candida albicans*,
protein = sterol 14α-demethylase / CYP51, structures carrying protein + heme):

| conformer | resolution | co-crystal azole | ligand code | role |
|---|---:|---|---|---|
| **5TZ1** | 2.0 Å | VT-1161 (oteseconazole precursor) | VT1 | the project's *current* receptor + reference ligand |
| **5FSA** | 2.86 Å | posaconazole | X2N | new conformer |
| **5V5Z** | 2.9 Å | itraconazole | 1YN | new conformer |

The RCSB query returned three further *C. albicans* hits that are **not CYP51** and are excluded
by construction (not by inspection of any result): 4ESW (Thi5 thiamine enzyme), 5UCC (ENTH
domain of ENT2), 9IUL (Cdr1 efflux pump). No apo *C. albicans* CYP51 structure exists; all three
conformers are azole-bound (holo) — see Limitation 2.

**These three, and only these three, are the ensemble. Adding a fourth structure — an *S.
cerevisiae* surrogate (5HS1, 5EQB, 4WMZ), an *Aspergillus* CYP51 (5FRB), or an MD snapshot —
is a different, separately-registered experiment, forbidden under the stopping rule of this one.**

## The self-docking hazard, and why the GATE is clean of it by construction

Cross-docking into holo structures risks the crudest possible leakage: docking a drug back into
the induced-fit pocket that was *solved with that very drug bound*. Two of the three co-crystal
azoles — **posaconazole (5FSA)** and **itraconazole (5V5Z)** — are in the project's actives set.
But they are in the **training** actives (`data/ligands/actives.smi`), and the primary gate uses
the **7 held-out** azoles (`data/ligands/actives_holdout_final.smi`: efinaconazole, luliconazole,
sertaconazole, oxiconazole, bifonazole, fenticonazole, butoconazole). **None of the 7 gate
actives was co-crystallised in any conformer.** So the primary gate is self-docking-leakage-free
by construction, not by a correction. (VT-1161 in 5TZ1 is the reference ligand, not a test
molecule, so it poses no leakage either.)

The n=15 powered secondary *does* re-introduce posaconazole and itraconazole. Pre-committed
guard: in the ensemble-minimum for those two molecules, their self-pair conformer is **excluded**
(posaconazole's ensemble-min is taken over {5TZ1, 5V5Z}; itraconazole's over {5TZ1, 5FSA}). The
n=15 result is reported *both* with and without this exclusion, and the exclusion is the headline.

## The feature — mode C, UNCHANGED, aggregated by a single pre-committed rule

No new feature. The criterion is the project's validated primary metric: **mode C = the minimum
ligand-nitrogen-to-heme-iron distance** — the same metric that passed the overall gate
(AUC 0.716, look #6) and defines the applicability domain. It is computed *per conformer* exactly
as `openafr.scoring` computes it today (ligand N to that conformer's own heme Fe, in that
conformer's coordinates), then aggregated across the 3 conformers by:

> **Ensemble-minimum (fixed, primary).** A molecule's ensemble mode C = the *smallest* N–Fe
> distance it achieves in *any* of the 3 conformers.

Mechanistic justification, fixed in advance: induced fit means the *productive* channel
conformation differs per ligand, so the right question is "**does there exist** an experimentally
observed CYP51 conformation in which this molecule coordinates the iron like a real azole?" That
is ensemble-minimum. It is the standard ensemble-docking aggregation, has **zero free
parameters** beyond the frozen `protocol.yaml`, and is chosen for that reason — not because it
maximises separation (it has never been evaluated). Ensemble-mean is reported as a diagnostic
(Secondary 2) but does **not** gate and will **not** be promoted if min fails (Stopping rule).

## The population — IDENTICAL to look #6 / look #7, so the comparison is a clean head-to-head

| | |
|---|---|
| actives (gate) | the **same 7 held-out azoles** graded in looks #6 and #7, all azole-bearing |
| inactives (gate) | the **azole-bearing** verified inactives only (~162), the *same* `azole_bearing()` selection used in `validate_gate_verified.py` — count is whatever that function reproduces |
| full inactive set | all **279** verified inactives are docked into every conformer too, so the novel-chemotype secondary (below) can be read under the ensemble |
| new docking | the 7 actives, 8 training actives, and 279 inactives are docked into **5FSA and 5V5Z** under the frozen protocol (seed 42). The 5TZ1 poses are reused from look #6 unchanged. |
| protocol | `protocol.yaml` UNCHANGED — box 26 Å @ (70.61, 66.28, 4.18), exhaustiveness 32, 20 modes, seed 42 |
| receptor prep | 5FSA/5V5Z chain A protein + chain-A HEM extracted (as `prep_receptor.py` does for 5TZ1), **superposed onto 5TZ1 chain A by Cα least-squares (Kabsch) fit** so the frozen box lands in each pocket, co-crystal azole stripped, `obabel -xr -p 7.4` → pdbqt. Each conformer's prepared receptor anchor is hashed and pinned in the header above, before any docking or scoring. |

The 5TZ1-alone column of this run must reproduce look #6's 0.650 on the identical pool — a
built-in pipeline sanity check (Secondary 1). If it does not, the run is void and the ensemble
result is not read: the discrepancy is a prep/scoring bug, not science.

## Primary hypothesis (H1) — the only thing that gates

Ranking the pooled molecules (7 held-out azole actives + the ~162 azole-bearing measured
inactives) by **ensemble-minimum mode C over {5TZ1, 5FSA, 5V5Z}** separates working azoles from
failed azole analogues at **AUC ≥ 0.70** (`protocol.yaml` `min_auc`, unchanged, not lowered).
AUC is the gate because it is threshold-free.

**EF does not gate** — the same arithmetic as looks #6/#7: 7 actives in a ~169-molecule pool
makes the top-1% slice ~2 molecules (effectively Bernoulli). EF@1/5/10% and BEDROC(α=20) are
reported for continuity only.

## Secondary — all pre-specified, none gating

1. **Per-conformer AUC** on the gate pool: 5TZ1 alone (must ≈ 0.650, the sanity check), 5FSA
   alone, 5V5Z alone. Shows whether any single new conformer already carries the signal, or
   whether the lift (if any) is genuinely an ensemble effect.
2. **Ensemble-mean aggregation** AUC on the gate pool — the diagnostic alternative to min.
   Reported, never promoted to the gate.
3. **Marginal-conformer increment**: AUC of {5TZ1}, then {5TZ1,5FSA}, then all three, in that
   fixed order. Answers "does each added conformer help, and is 3 enough or plateauing?"
4. **n=15 powered form** with the self-pair exclusion (see hazard section), reported with AUC and
   bootstrap CI. Does not gate (no matched single-conformer comparator at n=15, and it carries
   the residual self-docking caveat).
5. **Novel-chemotype guard (critical honesty check).** The ensemble is only worth adopting if it
   does not *break* the thing that already works. Report ensemble-min mode C AUC on the pool
   where the tool passes today — 7 (or 15) actives vs the **non-azole** verified inactives
   (comparator: AUC 0.810, look #6). An ensemble that lifts within-azole but collapses the 0.810
   is not an improvement; it has just made every molecule coordinate more easily. Both numbers
   are reported together.
6. **Statistical support**, exactly as look #6: permutation p (20,000 shuffles), bootstrap 95%
   AUC CI, active rank positions, H2 tip percentile. Descriptive; the gate is the point estimate.
7. **mode C single-conformer (0.650) and mode F (0.590)** re-printed alongside, so the ensemble
   result is visible against both prior FAILs in one view.

## Analysis plan

    Primary (gates):  ensemble-min mode C, {5TZ1,5FSA,5V5Z}, pooled 7 actives + ~162
                      azole-bearing inactives, AUC >= 0.70
    Reported:         per-conformer AUCs (5TZ1 must ~= 0.650); ensemble-mean; marginal-conformer
                      increment; n=15 powered (self-pair excluded); novel-chemotype guard AUC
                      (comparator 0.810); EF@1/5/10%, BEDROC(20); permutation p; bootstrap CI;
                      mode C=0.650 and mode F=0.590 reprinted
    Unrankable actives (no docked pose with a nitrogen in any conformer) ranked LAST, never dropped.
    Metrics from openafr.scoring — the same code that graded run 2, look #6, look #7.

Graded by `scripts/validate_gate_ensemble.py` (to be written), which re-checks this document's
SHA-256, the `protocol.yaml` hash, and the three pinned receptor hashes before reporting, and
refuses to certify a pass if any moved.

## Pre-committed interpretation

- **PASS — AUC ≥ 0.70, AND the novel-chemotype guard (Secondary 5) still ≥ ~0.78.** A real
  experimental conformer ensemble lifts within-azole ranking above the single-rigid-conformer
  ceiling *without* breaking novel-chemotype triage. This is the concrete evidence that the mode
  C and mode F FAILs were **induced-fit / single-conformer artefacts**, not failures of the
  geometric criterion — the explanation both prior preregs pre-committed as leading. The
  within-class triage product re-opens, still at n=7 with a wide CI; a separately pre-registered
  confirmation on an *independent* active set is required before any within-class wet-lab claim.
- **PASS on within-azole but the novel-chemotype guard collapses.** The ensemble did not add
  discrimination; it just made coordination easier for everything, inflating actives and
  decoys alike. Reported honestly as **not an improvement**, with both numbers. The single rigid
  conformer remains the product receptor.
- **PARTIAL — AUC between 0.650 and 0.70.** The ensemble helps but does not clear the bar: three
  crystallographic conformers are not enough flexibility. Pre-committed next step, if pursued, is
  an **MD-generated ensemble** (separately registered) — *not* re-tuning this one. Reported as a
  PARTIAL, never spun as a pass.
- **FAIL — AUC ≤ 0.650 (no better than the single conformer).** The within-azole ceiling
  survives a *complete, exhaustive experimental conformer ensemble*. This is strong evidence the
  limitation is **not** captured by crystallographic receptor flexibility, and materially
  strengthens the preprint's §4.3 induced-fit limitation: even the real observed spread of the
  CYP51 channel does not let rigid-family docking rank within the azole class. It also argues the
  within-class question may be beyond geometry-over-a-fixed-target entirely. This **closes the
  crystallographic-ensemble line**; MD / flexible-sidechain is a further separately-registered
  step only if independently justified, and the defensible product remains novel-chemotype
  triage (AUC 0.810).

## Stopping rule (binding)

**One ensemble, one aggregation rule, one criterion, one look.** The ensemble is the exhaustive
3-structure experimental set — no conformer added or removed after the fact. The aggregation is
ensemble-minimum — not swapped to ensemble-mean if min disappoints. The criterion is mode C,
unchanged — no new feature is defined on these new poses in this look (a new feature on the
ensemble poses would be a new pre-registration). No threshold is tuned (the gate is AUC). Any
further attack — MD ensemble, flexible sidechains, an independent active set — is a genuinely new
dataset and a separate pre-registration, with this result left standing as graded.

## Limitations (fixed in advance)

1. **Only 3 experimental conformers exist.** A thin ensemble. A FAIL refutes only that *these
   three crystallographic conformers* help — not induced-fit flexibility in general (which an MD
   ensemble could still probe). Stated as such, never overclaimed.
2. **All 3 are holo (azole-bound) structures.** The channel is already in azole-accommodating
   conformations, which if anything biases the test *toward* helping actives — so a FAIL is the
   more informative outcome, and a PASS must be read against this thumb on the scale. No apo
   *C. albicans* CYP51 structure exists to counterbalance it.
3. **Cross-docking into lower-resolution structures.** 5FSA (2.86 Å) and 5V5Z (2.9 Å) are softer
   than 5TZ1 (2.0 Å); their coordinates carry more noise, which can only *hurt* discrimination —
   so it is a conservative addition, not a flattering one.
4. **Ensemble-minimum inflates apparent coordination for decoys too** (more conformers = more
   chances at a short N–Fe). The fair comparator is therefore the single-conformer 0.650 on the
   *identical* pool, reprinted alongside (Secondary 7), not an absolute AUC in isolation.
5. **Superposition affects only the shared box, not the measured distance.** N–Fe is computed
   within each conformer's own coordinate frame; the Cα fit only places the frozen box in each
   pocket. Box occupancy per conformer is reported so a mis-seated box is visible.
6. **Dataset caveats carry over unchanged** from look #6: whole-cell MIC ≠ CYP51 enzyme assay,
   ChEMBL inactivity is heterogeneous, publication bias, and the `azole_bearing()` heuristic.
7. **One docking program, rigid-within-each-conformer, single seed.** "Ensemble" here means
   discrete rigid conformers, not continuous flexibility; sidechains do not relax during docking.

## Implementation manifest (BUILT 2026-08-25, before scoring)

1. `data/structures/5FSA.pdb`, `data/structures/5V5Z.pdb` — raw RCSB deposits (sha256 recorded
   by the prep script's stdout; provenance is this document + the RCSB accessions).
2. `scripts/prep_receptor_ensemble.py` — extracts chain A protein + chain-A HEM per conformer,
   **Kabsch-superposes CA onto the 5TZ1 anchor**, strips the co-crystal azole, writes the tracked
   `work/receptor_{5FSA,5V5Z}_A.pdb` anchor + `.pdbqt`, and prints each hash (pinned above).
3. `openafr/ensemble.py` — ensemble aggregation (ensemble-min / ensemble-mean over per-conformer
   mode C), pure functions locked by `tests/test_ensemble.py`; no change to the pose scorer.
4. `scripts/validate_gate_ensemble.py` — grader that hash-gates this doc + `protocol.yaml` + the
   3 receptor anchors, reads the 3 conformer screens, and emits the analysis-plan table.
5. Docking: the 279 inactives + 15 actives are docked into 5TZ1, 5FSA and 5V5Z (seed 42, frozen
   protocol). Est. ~880 dockings (5TZ1 poses were gitignored/gone, so all three are re-docked);
   by look #7's rate (~294 dockings / 2.5 h on 12 cores) ≈ **7-8 h wall on Apple Silicon —
   local, no VM required.**

## Reproduce (once built + frozen)

```
conda activate openafr
# receptors: 5TZ1 as today, plus the two new conformers superposed onto it
python scripts/prep_receptor.py                                   # 5TZ1 (unchanged)
python scripts/prep_receptor_ensemble.py 5FSA 5V5Z               # new: align + prep
# ligands: reuse look #6 prep; dock the two new conformers (5TZ1 poses reused)
python scripts/prep_ligands.py data/ligands/actives_holdout_final.smi work/verified_ligands
python scripts/prep_ligands.py data/ligands/verified_inactives.smi    work/verified_ligands
python scripts/prep_ligands.py data/ligands/actives.smi               work/verified_ligands
scripts/screen.sh work/verified_ligands work/screen_5FSA work/receptor_5FSA.pdbqt
scripts/screen.sh work/verified_ligands work/screen_5V5Z work/receptor_5V5Z.pdbqt
# one look:
python scripts/validate_gate_ensemble.py \
    work/screen_verified work/screen_5FSA work/screen_5V5Z       # exits 0 pass / 1 fail
```
</content>
</invoke>
