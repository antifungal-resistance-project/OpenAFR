# Results — per-candidate stereochemistry axis (#85)

Date: 2026-08-27. Pre-registered in
[PREREGISTRATION_candidate_stereo.md](PREREGISTRATION_candidate_stereo.md)
(sha256 cdd6882a…, frozen BEFORE any per-isomer distance was read;
`work/PREREG_candidate_stereo.sha256`). Protocol hash unchanged from run 2.

> **This is a per-candidate confidence axis, not a gate.** It annotates the scorecard; it
> never promotes. It answers issue #85: is each shortlist rank tied to a single **defined**
> stereoisomer, or could an unspecified stereocenter be silently making or breaking the
> iron-coordinating geometry the rank rests on?

## What was tested

Every one of the 16 scorecard molecules was audited straight from its supplied SMILES
(`work/candidates/shortlist_focus.smi`) with RDKit — tetrahedral stereocenters (assigned +
unassigned) and unspecified E/Z bonds — and its distinct isomers from the **unassigned**
centers were enumerated (`scripts/enumerate_stereo.py`, deterministic, docking-free).
`stereo_explicit` ≡ the supplied SMILES pins a single isomer.

For the molecules that were **not** stereo-explicit, every enumerated isomer was embedded
(ETKDGv3, fixed seed) and docked into the **same wild-type receptor the shortlist rank was
measured on** (`work/receptor.pdbqt`), under the frozen protocol, 5 seeds {42,1,2,3,4} — the
**only** variable versus the shortlist dock is the ligand stereochemistry. Per-isomer
consensus N–Fe uses the frozen `consensus()` logic unchanged.

Verdicts (frozen): **EXPLICIT** = single defined isomer, no dock needed; **ROBUST** = every
isomer coordinates (≥3/5 seeds) and the per-isomer consensus spread ≤ 1.0 Å;
**ISOMER-DEPENDENT** = the isomers disagree (some coordinate, some don't) or the spread > 1.0 Å.

## The result

**14 of 16 molecules are stereo-explicit — including all 10 novel candidates.**

| candidate | isomers | unassigned centers | verdict |
|---|:--:|:--:|---|
| verinurad, vatalanib, flucloxacillin, taranabant, lersivirine, LDN-27219, R-1479, olprinone, nolatrexed, L-838417 | 1 | 0 | **EXPLICIT** |
| ravuconazole*, voriconazole*, letrozole*, fluconazole* | 1 | 0 | **EXPLICIT** |
| flutrimazole* | 2 | 1 | **ROBUST** |
| ketoconazole* | 4 | 2 | **ISOMER-DEPENDENT** |

*known-azole controls.

The two molecules with unspecified stereochemistry are **both controls**, docked per-isomer to
exercise the axis:

| control | isomer | consensus N–Fe (Å) | seeds coord |
|---|---|---:|:--:|
| flutrimazole | iso01 | 2.541 | 5/5 |
| flutrimazole | iso02 | 2.563 | 3/5 |
| ketoconazole | iso01 | 2.548 | 3/5 |
| ketoconazole | iso02 | 2.720 | 2/5 |
| ketoconazole | **iso03** | **3.342** | **0/5** |
| ketoconazole | iso04 | 2.841 | 4/5 |

## What this shows

- **Every novel candidate's rank is stereochemistry-explicit.** The supplied SMILES pins a
  single defined isomer for all 10 — seven are achiral (verinurad, vatalanib, lersivirine,
  LDN-27219, L-838417, olprinone, nolatrexed) and three carry fully-assigned centers
  (flucloxacillin 3, taranabant 2, R-1479 4). So the specific failure #85 names — a rank
  silently averaged over, or pinned to the wrong member of, an *unspecified* stereo set —
  **does not occur for any molecule on the shortlist.** No isomer docking was needed for the
  candidate path; the finding *is* the result.
- **The axis is not vacuous — it distinguishes the two ambiguous controls.**
  **flutrimazole is ROBUST**: its two enantiomers coordinate to within 0.022 Å (2.541 vs
  2.563 Å), so its rank is isomer-insensitive. **ketoconazole is ISOMER-DEPENDENT**: one of
  its four enumerated isomers (iso03) *never reaches the iron* (0/5 seeds, 3.342 Å) while the
  others coordinate at 2.5–2.8 Å. Had ketoconazole been a novel candidate supplied without
  configuration, *which* isomer happened to be docked would have changed whether it looked
  like a coordinator at all — exactly the risk this axis exists to catch. (ketoconazole is
  marketed as the *cis* diastereomer pair; the enumeration is deliberately blind to that
  external knowledge, which is why it also generates the non-coordinating *trans* forms.)
- **The scorecard gains a `stereo` column.** Novel candidates read `explicit`; flutrimazole
  `robust/2iso`; ketoconazole `iso-dep/4iso`. The ISOMER-DEPENDENT flag is a demotion for a
  *novel* candidate (concern: "rank depends on unspecified stereochemistry") but
  control-exempt — the ambiguous molecules here are the racemic azole controls, present to
  validate the axis, not discoveries to strike.

## What this establishes / does not

- **Does:** confirm every shortlist rank refers to one defined isomer, remove "unspecified
  stereochemistry" as a hidden confound for the novel candidates, and demonstrate on the
  controls that the isomer-docking machinery genuinely separates an isomer-robust molecule
  (flutrimazole) from an isomer-dependent one (ketoconazole).
- **Does not:** second-guess a *supplied* configuration. A molecule whose SMILES asserts
  `@`/`@@` is taken at its word — this axis closes the "unspecified stereo" hole, not the
  "supplied the wrong configuration" one. It also does not claim the docked isomer is the
  bioactive one, nor relax any caveat of the underlying screen (single conformer, rigid
  receptor, coordination necessary-not-sufficient). Protonation/tautomer state is the separate
  axis #84 and remains `not-yet-checked`.

## Reproduce

    conda run -n openafr python scripts/prep_receptor.py            # -> work/receptor.pdbqt
    conda run -n openafr python scripts/enumerate_stereo.py         # audit + isomer SDFs
    conda run -n openafr bash -c 'RECEPTOR=work/receptor.pdbqt \
        scripts/screen_consensus.sh work/candidates/stereo_ligands \
        work/candidates/dock_stereo "42 1 2 3 4"'                   # dock ambiguous isomers
    conda run -n openafr python scripts/candidate_stereo.py --json \
        > work/candidates/shortlist_stereo.json                     # grade
    conda run -n openafr python work/candidates/build_scorecard.py  # fold into the scorecard

The enumeration + grader are committed; the per-seed dock tree
(`work/candidates/dock_stereo/`) is regenerated by the screen and gitignored like the other
dock trees.
