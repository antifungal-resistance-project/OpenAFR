# Pre-registration — run 3: VT-1598 resistance-active probe in the Y132F pocket

> **STATUS: DRAFT — NOT FROZEN, and DE-PRIORITIZED (2026-08-22).** This document is not
> sha-pinned and no docking has been run against it. It is kept frozen-ready but the T1
> widening pass concluded it should **not** be run now: the verifiable in-domain
> resistance-active pool is n=1 (VT-1598), VT-1598 sits at the ≤45-heavy-atom domain edge, and
> resistance-retention is confounded with over-domain molecular size
> (`data/ligands/PROVENANCE.md`). A single-probe run is therefore low-yield — a null is
> confounded by the domain limit, a hit is an n=1 anecdote. **Do not dock against this or treat
> its hash as a commitment** unless a flexible-receptor / induced-fit protocol is later added
> (which would lift the domain ceiling that makes the probe uninformative today), at which
> point replace this banner with a "frozen <date>" line. The project's full pre-registration
> discipline (fixed hypotheses, both outcomes pre-committed, no post-hoc metric swaps) applies
> the moment it freezes. The T1 deliverable is the confound finding, which feeds the T3
> wet-lab package.

Inherits the base protocol from [PREREGISTRATION_run2.md](PREREGISTRATION_run2.md) and the
mutant receptor from [PREREGISTRATION_auris_Y132F.md](PREREGISTRATION_auris_Y132F.md); read
both first. Only the deltas are stated here.

## Why this run is scoped as a single-probe case study, not a differential validation

T1 (`docs/designs/refocus-resistant-cyp51-triage.md`) set out to expand the active set with
CYP51 inhibitors that **retain measured activity against azole-resistant *C. auris***, so a
resistance-aware criterion could be tested by ranking *still-active* compounds above the
*defeated* azole panel in the mutant pocket. The 2026-08-22 literature pass found the
verifiable pool for that is **n = 1**:

- **VT-1598** (PubChem CID **126715974**; PMID 30530603, modal MIC 0.25 µg/mL across 100
  clinical *C. auris* isolates incl. azole-resistant) — the one CYP51-mechanism compound with
  a primary-literature *C. auris* MIC **and** a verifiable PubChem structure.
- Oteseconazole (VT-1161, CID 77050711) and quilseconazole (VT-1129, CID 91886002) — dropped,
  no *C. auris*-specific MIC. The 2024–25 resistance-active hits are unnamed research compounds
  with no PubChem registration → fail discipline constraint #3 (verified SMILES).

One resistance-active compound **cannot** support a differential-activity permutation test (you
cannot test whether a *class* ranks above another with a single member). So this run does the
fallback the refocus plan explicitly named (unresolved decision #2): **it validates broad
enrichment on the mutant pocket and carries VT-1598 as a single, pre-registered probe.** The
deliverable is an honest case study — "does the geometry criterion place the one known
resistance-active compound among the top coordinators in the Y132F pocket?" — **never** a
statistical resistance-breaker claim.

## New active (fixed, verified)

    name          VT-1598
    PubChem CID   126715974
    formula       C31H20F4N6O2   MW 584.5   heavy atoms 43
    SMILES        C1=CC(=CC=C1COC2=CC=C(C=C2)C#CC3=CN=C(C=C3)C([C@](CN4C=NN=N4)(C5=C(C=C(C=C5)F)F)O)(F)F)C#N
    novelty       Morgan r=2, 2048-bit max Tanimoto 0.312 (nearest fluconazole) vs all 15
                  existing actives → PASS (< 0.70); rdkit 2025.09.5 under conda env openafr

VT-1598 is docked into the **Y132F mutant receptor** (`work/receptor_auris_Y132F.pdbqt`,
sha256 924cf7b8…, unchanged from the Y132F run) under the **frozen run-2 protocol** (26 Å box
@ 70.61/66.28/4.18, exhaustiveness 32, seed 42) — no re-tuning.

## Applicability-domain risk — DECLARED IN ADVANCE (this is the load-bearing caveat)

VT-1598 has **43 heavy atoms / MW 584** — the largest molecule in the entire active set (every
existing active is ≤ 36) and at the edge of run-2's declared domain (valid ≤ 45 heavy atoms;
posaconazole/itraconazole class ≥ 52 failed reproducibly in run 1 from induced fit a rigid
receptor cannot model). VT-1598's long ethynyl-phenoxy-benzonitrile arm is exactly the kind of
extended flexibility a rigid pocket handles poorly.

**Pre-committed rule for a no-pose / poor-pose VT-1598 result:** if VT-1598 produces no
iron-coordinating pose (or ranks in the bottom tier), that is reported as **"result lies
outside / at the edge of the validated applicability domain"** — it is **NOT** counted as
evidence the criterion fails on resistance-active chemistry, and it is **NOT** spun as a pass
either. A rigid-receptor false negative on an over-domain ligand is uninformative about the
criterion; saying so up front is what keeps a null honest instead of convenient.

## Primary readout — VT-1598 probe (fixed before docking)

For a single compound the enrichment metrics (AUC/EF) do not apply; the pre-committed readouts
are its **placement** among a property-matched decoy pool in the mutant pocket:

    (a) minimum N–Fe distance of VT-1598's best pose (Å)
    (b) VT-1598's percentile rank by min N–Fe among its property-matched decoys

Decoy pool: 50 property-matched presumed-inactive decoys generated for VT-1598 by
`scripts/make_decoys.py` under the unchanged criteria (MW ±25, cLogP ±1.5, rot ±2, HBD ±1,
HBA ±2, ≥1 aromatic N, Morgan Tanimoto < 0.35 to VT-1598). Decoy generation is rdkit-only
(local); docking is cloud-gated.

Reference band for "coordinating," from the graded Y132F run: real azoles reach ~2.7–3.2 Å
(efinaconazole 2.72 → oxiconazole 3.27); the top decoys there reached 2.62–2.72 Å.

## Pre-committed interpretation — BOTH outcomes fixed now

- **Probe COORDINATES (VT-1598 min N–Fe ≤ ~3.2 Å AND top-decile among its decoys).** The
  geometry criterion places the one known resistance-active compound among the top coordinators
  in the resistant pocket. Reported as a **supportive single-compound case study** — evidence
  the criterion is *consistent with* (not proof of) resistance-aware ranking; n=1, anecdote,
  stated as such.
- **Probe DOES NOT COORDINATE (min N–Fe > ~3.2 Å OR below-median among its decoys), in-domain.**
  The criterion does not rank the known resistance-active compound highly. Reported as a
  **negative case study** — informative, and a reason not to claim resistance-awareness.
- **Probe returns no pose / over-domain (the applicability-domain branch above).** Uninformative;
  reported as a domain-limit result, neither pass nor fail.

No metric or ranking mode invented after seeing these numbers may be reported as a pass
(the run-2 / Y132F rule carries over).

## Secondary (report-only, does not gate): broad-enrichment confirmation

Optionally re-report the 7-azole broad-enrichment on the mutant pocket on the **durable
metrics only** — AUC and EF@5% (per [RESULTS_reliability.md](RESULTS_reliability.md), EF@1% is
a single-shot docking artifact on this decoy set and is **not** used). This is confirmation of
the already-graded Y132F result (AUC 0.805 / EF@5% 5.63x), not a new claim, and it does not
gate run 3.

## Declared limitations (fixed in advance)

1. **n = 1 resistance-active.** This is a single-compound probe, not a validation of a
   resistance-active class. It cannot support a permutation test or any statistical
   resistance-aware claim, and must never be pitched as one.
2. **Applicability-domain edge** (see above) — the dominant threat to a clean read.
3. **Rigid Y132F receptor** — one *C. albicans* 5TZ1 scaffold + one deletion mutation, not a
   real *C. auris* ERG11 structure; inherits every Y132F-run limitation.
4. **Decoys presumed inactive** (unchanged).

## Reproduce (once frozen; docking is cloud-gated)

```
conda activate openafr
# 1. write the new active + generate its decoys (local, rdkit only)
printf '%s\tVT-1598\n' 'C1=CC(=CC=C1COC2=CC=C(C=C2)C#CC3=CN=C(C=C3)C([C@](CN4C=NN=N4)(C5=C(C=C(C=C5)F)F)O)(F)F)C#N' \
    > data/ligands/actives_run3.smi
python scripts/make_decoys.py data/ligands/actives_run3.smi data/ligands/decoys_run3.smi
# 2. dock into the Y132F mutant pocket (cloud Linux/x86 VM — vina/obabel absent locally)
# 3. grade the probe placement (metric = VT-1598 percentile + min N-Fe), verifying this
#    doc's sha256 + protocol hash first, as every gate does.
```
