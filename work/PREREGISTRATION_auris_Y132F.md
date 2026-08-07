# Pre-registration — *C. auris* ERG11 Y132F resistance-pocket transfer test

**Written 2026-08-05, BEFORE docking anything into the mutant receptor.** Everything below
is fixed in advance. Deviations must be recorded as deviations, not edited into this
document. This mirrors the discipline of [PREREGISTRATION_run2.md](PREREGISTRATION_run2.md);
read that first for the base protocol this inherits.

## Background and what this run does *not* claim

The run-2 gate ([RESULTS_run2_holdout.md](RESULTS_run2_holdout.md)) validated, on a blind
held-out set, that ranking azoles by minimum nitrogen-to-heme-iron distance separates them
from property-matched decoys (mode C: AUC 0.794, EF@1% 12.68x) in the **wild-type
*C. albicans* 5TZ1 pocket**. That result is not re-litigated here and is not re-used as
evidence for this run.

This run asks a **different, narrower question**: does that geometric criterion still
work when the pocket carries the dominant *C. auris* azole-resistance mutation **Y132F**
(ERG11; clades I and IV)? Per [VISION.md](../VISION.md), the resistant pocket — not the
wild type — is the target the project actually cares about. This is the first test of
whether the WT-validated criterion **transfers** to a resistant target.

It is a **controlled single-variable experiment**: the same held-out actives, the same
decoys, the same frozen docking protocol as run 2 — the *only* change is the receptor side
chain at position 132 (Tyr → Phe).

## Mutant receptor (fixed, deterministic)

    source        work/receptor_A.pdb  (run-2 WT anchor, sha256 a1a8410868d1)
    mutation      Y132F  (Tyr132 -> Phe132), by deletion of the para-hydroxyl OH only
    tool          scripts/mutate_receptor.py  (deterministic; no rotamer choice)
    anchor        work/receptor_auris_Y132F.pdb
                  sha256 924cf7b81b206a386f8007b4b54b2add14add89665393d45b4f648f6c1b2355f
    docking input work/receptor_auris_Y132F.pdbqt, built by:
                  python scripts/prep_receptor.py --mutant work/receptor_auris_Y132F.pdb \
                      work/receptor_auris_Y132F.pdbqt

Every atom except the deleted Tyr132 OH keeps its exact wild-type coordinate; the heme
iron is unchanged (verified at 64.163, 71.316, 2.711). The docking box, seed, search
effort and applicability domain are inherited **unchanged** from the frozen
`protocol.yaml` / run-2 pre-registration — this is not a re-tuned protocol.

## Primary hypothesis (H1-auris)

Ranking molecules by **minimum nitrogen-to-heme-iron distance** (smaller = better) in the
**Y132F mutant pocket** will separate the run-2 held-out azole antifungals from their
property-matched decoys with:

    AUC     >= 0.70
    EF@1%   >= 5.0x

Both must hold. Same bar as run 2 — it is not being lowered or raised for this target.

## Held-out actives and decoys (unchanged from run 2)

The identical 7 held-out actives (`data/ligands/actives_holdout_final.smi`, sha256
0d964c88…) and their 50-per-active property-matched decoys are re-docked into the mutant
pocket. Nothing about the ligand set is regenerated; that is what makes the receptor the
only variable. Decoys are presumed inactive (a known limitation, unchanged).

## Analysis plan (unchanged from run 2)

    Primary (gates):  mode C — rank by minimum N-Fe distance across all poses
    Secondary (report only):  mode A (Vina score), mode D (combined rank)
    Metrics: AUC, EF@1/5/10%, BEDROC(alpha=20) from rdkit.ML.Scoring.
    Any active producing no pose is ranked LAST, never dropped.

Graded by `scripts/validate_gate_auris.py`, which re-checks this document's SHA-256 and the
frozen protocol hash before reporting, exactly as the run-2 gate does.

## Pre-committed interpretation — BOTH outcomes fixed in advance

The result of mode C on this mutant pocket will be reported as the outcome, whatever it is.
Its meaning is committed **now**, before the numbers exist, so neither result can be
spun after the fact:

- **PASS (AUC >= 0.70 and EF@1% >= 5.0x).** The geometric criterion still identifies
  azole iron-coordination in the Y132F pocket. The tool's validated applicability
  **extends to this resistant target**, and screening a novel library against Y132F is
  justified under this protocol. This is a positive transfer result — *not* a claim that
  the tool "defeats resistance"; see limitations.

- **FAIL (AUC < 0.70 or EF@1% < 5.0x).** The criterion does **not** transfer: in the
  Y132F pocket the known azoles no longer produce iron-coordinating poses distinguishable
  from decoys. This is a **publishable, informative negative** — consistent with Y132F
  reshaping the pocket so the coordinating geometry is lost (a structural correlate of the
  clinical resistance), and it means the WT-validated tool **must not** be silently applied
  to the resistant target. Per the run-2 rule, no ranking method invented after seeing
  these numbers may be reported as a pass.

## Declared limitations (fixed in advance, true of either outcome)

1. **A rigid single-point deletion under-models real resistance.** Clinical Y132F
   resistance also involves conformational change, H-bond-network rearrangement, and often
   co-occurring efflux/expression changes. A rigid Tyr→Phe deletion captures the
   pocket-geometry change only. So a PASS does not prove clinical efficacy against
   resistant strains, and a FAIL is not proof the criterion is worthless — each is
   evidence about *this specific structural perturbation*, no more.
2. **Same receptor as run 2 otherwise.** One rigid *C. albicans* 5TZ1 scaffold carrying
   one *C. auris* mutation — not a full *C. auris* ERG11 structure (none exists
   experimentally; see `.context/auris_port_scope.md`, Path B).
3. **Y132F only.** K143R (clade I) and F126L (clade III) are not tested here; they require
   a side-chain repacker not yet in the pinned environment.
