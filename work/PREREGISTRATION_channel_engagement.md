# Pre-registration — the within-azole question: channel engagement, an orthogonal pose feature

**Written 2026-08-24, BEFORE the feature was computed on a single pose.** The channel-shell
residue set (below) is fixed from the co-crystal structure with **no active/inactive label
involved**, and the ranking feature has never been evaluated on any docked azole. Everything
below — the feature, the population, the gate, the interpretation of both outcomes — is fixed in
advance and frozen by this document's own SHA-256, which the grading script re-checks before it
will certify anything.

## Background — the gap this closes, and why one more look at these poses is permitted

Look #6 (`work/RESULTS_verified_inactives.md`) tested the geometric criterion against **279
compounds measured not to inhibit *C. albicans***. It passed overall (AUC 0.716) but located the
ceiling precisely, in a *pre-specified* subgroup split:

    vs non-azole chemotypes             AUC 0.810   (the tool works)
    vs azole analogues measured inactive AUC 0.650   (below its own 0.70 bar)

Most of what mode C (minimum ligand-nitrogen-to-heme-iron distance) detects is **chemotype, not
potency**: of the 39 molecules inside the validated iron-bound band, 33 are failed azoles and 2
are true actives. Every downstream ambition — a candidate worth an email, the parked wet-lab
package ([[outreach-hold-build-first]]) — needs ranking *within* the azole class, and that is the
0.650. This is now, for the first time, a **well-posed supervised problem**: azole actives vs a
large set of azole-bearing compounds *measured* not to work, all docked under the frozen protocol.

**Why a second feature-look at the verified-inactive poses is not fishing.** The one-look
discipline forbids *fishing* — computing many features and reporting the one that separates. It
does not forbid a single, mechanistically-fixed, orthogonal feature named before it is computed
and reported whatever it yields. That is exactly the standard the axial pre-registration
(`work/PREREGISTRATION_axial.md`, mode E) met when it reused already-docked poses for one new
feature. This document meets it identically: **one** feature, **no** new docking of the frozen
inactive set, **no** threshold to tune (the gate is AUC, which is threshold-free), and a hard
stopping rule (below). The feature is not a function of iron distance or iron *direction* — the
two axes already foreclosed — so it is a genuinely new signal, not a re-aggregation of an old one.

I confirm mode F has not been computed on any verified-inactive pose. The grader's SHA-gate is
what makes that claim checkable: a prereg edited after the result is worth nothing, and this
document's hash is verified at run time.

## The feature (mode F — channel engagement; fixed, and demonstrably not iron distance)

Azole inhibition needs iron coordination, but *every* azole coordinates the iron — that is why
mode C cannot separate working from failed azoles. What azole SAR actually turns on is the **N1
substituent / tail**: the working drugs thread a lipophilic tail through the CYP51 substrate-
access channel and make contacts a stubby, failed analogue does not. That occupancy is the
signal mode C is blind to. Concretely, per docked pose:

1. **Channel-shell residues (fixed once, label-blind).** From `work/receptor_A.pdb`, take every
   protein residue with ≥ 1 heavy atom within **4.5 Å** (`openafr.mapping.CONTACT_CUTOFF`,
   unchanged) of the co-crystallised azole **VT1** (`work/ref_VT1.pdb`, the reference ligand
   already used throughout the mapping stack). This is the set of residues the *real, bound
   drug* touches — defined from the crystal structure alone, using **no** knowledge of which
   test molecules are active. The heme is excluded automatically (it is a HETATM group;
   `read_residues` reads protein ATOM records only), so this is the *protein channel*, not the
   coordination shell.
2. **Which pose (fixed, not a free choice).** Engagement is measured in the molecule's
   **coordinating pose** — the single docked pose with the minimum ligand-nitrogen-to-iron
   distance, i.e. the *exact same pose* mode C selects. Tail occupancy is only meaningful in a
   pose that actually coordinates the iron, and tying mode F to the mode-C pose removes any
   min-vs-max-over-poses aggregation degree of freedom. A molecule with no nitrogen-bearing pose
   (no coordinating pose to select) is unrankable and placed last, exactly as under mode C.
3. **Engagement count.** In that pose, count the **distinct channel-shell residues** the ligand
   contacts — a residue counts if any ligand heavy atom lies within 4.5 Å of any of that
   residue's heavy atoms. Counting *residues* (each at most once), not ligand atoms, is
   deliberately less sensitive to sheer molecular size than an atom count would be.
4. **Rank.** Descending — more channel residues engaged = more drug-like tail occupancy = ranked
   higher. No threshold, no fitted weight, no per-metric search. The feature has exactly zero
   free parameters beyond the repo-wide 4.5 Å contact cutoff, which is fixed and pre-existing.

**Orthogonality is asserted in advance and will be reported, not assumed.** The grader reports
Pearson r between mode F and mode C (the N–Fe distance) on the pool. The mechanistic claim is
that a molecule can coordinate the iron tightly (short N–Fe) yet engage few channel residues
(no tail), or vice-versa, so the two are not redundant. If |r| turns out high, mode F adds
nothing and that is itself a reportable, committed finding — not a reason to swap features.

## The population — matched to look #6 so the comparison is exact

| | |
|---|---|
| actives | the **same 7 held-out azoles** graded in look #6 (`data/ligands/actives_holdout_final.smi`), all 100% azole-bearing |
| inactives | the **azole-bearing** verified inactives only (~162), selected by the *same* `azole_bearing()` function already in `validate_gate_verified.py` — the count is whatever that function reproduces, not a number fixed by hand here |
| poses | reproduced deterministically (seed 42) from `work/screen_verified/`; **no new docking of the inactive set** |
| protocol | `protocol.yaml` unchanged — box 26 Å @ (70.61, 66.28, 4.18), exhaustiveness 32, 20 modes, seed 42 |
| receptor | `work/receptor_A.pdb` → `work/receptor.pdbqt`, unchanged |

This is the *identical pool* on which mode C scored **AUC 0.650** in look #6. So the result is a
clean head-to-head: on the exact molecules where iron-distance fails to rank within the class,
does channel engagement do better? That comparator (0.650) is quoted, not recomputed.

## Primary hypothesis (H1) — the only thing that gates

Ranking the pooled molecules (7 held-out azole actives + the ~162 azole-bearing measured
inactives) by **mode F** separates working azoles from failed azole analogues at **AUC ≥ 0.70**
(`protocol.yaml` `min_auc`, unchanged, not lowered). AUC is the gate because it is threshold-free
— there is no cutoff to tune, which removes the last researcher degree of freedom.

**EF does not gate**, for the arithmetic reason look #6 already established: with 7 actives in a
~169-molecule pool the top-1% slice is ~2 molecules (effectively Bernoulli) and EF@1% is not
interpretable in either direction. EF@1/5/10% and BEDROC(α=20) are reported for continuity only.

## Secondary, all pre-specified, none gating

1. **Powered form (n=15).** The 8 *training* azoles (`data/ligands/actives.smi`) are held-out
   with respect to *this* feature — mode F was built from the co-crystal, never from them — so
   adding them (15 actives vs the same ~162 inactives) is a legitimate, better-powered read of
   the same question. Reported with AUC and bootstrap CI. It does **not** gate: there is no
   matched mode-C comparator at n=15, so the 7-active matched pool is the honest gate.
2. **Orthogonality diagnostics.** Pearson r(mode F, mode C) and r(mode F, ligand heavy-atom
   count), both on the pool, reported before any pass/fail is declared.
3. **The size-adjusted read.** Docking score correlates with heavy-atom count at r = −0.91 on
   this system, so "channel engagement" could be re-measuring molecular weight. Reported: the
   AUC of the *residuals* of mode F after regressing it on ligand heavy-atom count. This is what
   distinguishes "a size-independent fit signal" from "bigger molecules touch more of the
   channel," and it governs which PASS interpretation applies (below).
4. **Statistical support**, exactly as look #6: permutation p (20,000 shuffles), bootstrap 95%
   AUC CI, active rank positions, and the H2 tip percentile. Descriptive; the gate is the AUC
   point estimate against the frozen bar.
5. **mode C on this exact pool** re-printed alongside, so mode F vs 0.650 is visible in one view.

## Analysis plan

    Primary (gates):  mode F, pooled 7 actives + ~162 azole-bearing inactives, AUC >= 0.70
    Reported:         n=15 powered form; r(F,C) and r(F, heavy-atoms); size-residualized AUC;
                      EF@1/5/10%, BEDROC(20); permutation p; bootstrap 95% CI; H2 tip percentile;
                      mode C on the same pool (comparator 0.650)
    Unrankable actives (no docked pose / no channel contact) are ranked LAST, never dropped.
    Metrics from openafr.scoring — the same code that graded run 2 and look #6.

Graded by `scripts/validate_gate_channel.py`, which re-checks this document's SHA-256 and the
frozen `protocol.yaml` hash before reporting and refuses to certify a pass if either moved.

## Pre-committed interpretation

- **PASS (AUC ≥ 0.70) with the size-residualized AUC also ≥ 0.70.** Pose geometry *can* rank
  within the azole class; the 0.650 ceiling was a property of the *iron-distance* criterion, not
  of the rigid pose itself. The within-class triage product is re-opened — carefully, and still
  at n = 7 with a wide CI. A separately pre-registered confirmation on an independent active set
  would be the required next step before any wet-lab shortlist is claimed within-class.
- **PASS on raw mode F but the size-residualized AUC collapses below the bar.** The honest claim
  narrows to *"channel occupancy, largely tracking molecular size, enriches actives"* — real, but
  not a size-independent fit signal, and of limited use for triaging same-size analogues. Stated
  with the residualized number, never hidden behind the raw AUC.
- **FAIL (AUC < 0.70).** The within-azole ceiling is **not specific to iron distance**: an
  orthogonal, mechanistically-motivated tail-fit feature also cannot separate working from failed
  azoles on a single rigid conformer. This strengthens the preprint's §4.3 induced-fit /
  single-conformer limitation and is the decisive evidence for **Option B** — re-scope the
  product to novel-chemotype triage, where the measured AUC is 0.810. The next attack on the
  within-class question would then have to be a *flexible / ensemble receptor* (new poses = a new
  dataset), which is the parked direction in [[outreach-hold-build-first]] — not another feature
  on these rigid poses.

## Stopping rule (binding)

**One look, one feature, no re-filtering.** If mode F disappoints, it will not be re-defined with
a different contact cutoff, a different channel-shell definition, atom-counts instead of
residue-counts, or a fitted weight, and re-run. Exactly one orthogonal feature is spent on the
verified-inactive poses here; the axial precedent spent one on the presumed-decoy poses. Any
further within-class attempt requires a genuinely new dataset (an independent active set, or
ensemble/flexible-receptor poses), separately pre-registered, with this result left standing as
graded.

## Limitations (fixed in advance)

1. **n = 7 actives on the gate.** Every small-count caveat from look #6 carries over, including a
   bootstrap AUC CI that may straddle the bar. The n=15 secondary mitigates power but is not the
   gate (no matched comparator).
2. **A single rigid conformer may not even contain the tail interactions the feature is meant to
   read** (§4.3). A tail that would thread the channel under induced fit can be modelled folded
   back by rigid docking, which would blunt mode F through no fault of the hypothesis. A FAIL is
   therefore evidence about *rigid-pose geometry*, not a refutation of the SAR mechanism itself.
3. **The channel-shell is defined from one co-crystal ligand (VT1).** A different reference azole
   might touch a slightly different residue set. VT1 is chosen because it is the structure already
   pinned throughout the mapping stack, not because it maximises separation (it has never been
   scored against the test molecules).
4. **Whole-cell MIC ≠ CYP51 assay; ChEMBL inactivity is heterogeneous; publication bias** — the
   three dataset caveats from look #6's pre-registration carry over unchanged.
5. **One docking program, one conformer, one rigid receptor, single seed.** Nothing here revisits
   induced fit or resistance-breaking.

## Reproduce

```
conda activate openafr
# The verified-inactive poses reproduce exactly as in look #6 (no new inactive docking):
python scripts/make_verified_inactives.py fetch
python scripts/make_verified_inactives.py build
python scripts/prep_receptor.py
python scripts/prep_ligands.py data/ligands/actives_holdout_final.smi work/verified_ligands
python scripts/prep_ligands.py data/ligands/verified_inactives.smi   work/verified_ligands
scripts/screen.sh work/verified_ligands work/screen_verified
# The n=15 powered secondary additionally needs the 8 training actives docked under the
# identical protocol (deterministic, seed 42):
python scripts/prep_ligands.py data/ligands/actives.smi work/verified_ligands
scripts/screen.sh work/verified_ligands work/screen_verified
# One look:
python scripts/validate_gate_channel.py work/screen_verified   # exits 0 pass / 1 fail
```
