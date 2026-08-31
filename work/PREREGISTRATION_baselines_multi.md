# Pre-registration — multi-baseline comparison (issue #83)

**Written 2026-08-30, BEFORE running any baseline scorer of this set.** Everything below is
fixed in advance. It reuses the run-2 protocol, receptor and held-out poses unchanged and adds
a *table* of independent scorers applied to the SAME poses. Read
[PREREGISTRATION_run2.md](PREREGISTRATION_run2.md), [RESULTS_run2_holdout.md](RESULTS_run2_holdout.md)
and [PREREGISTRATION_baselines.md](PREREGISTRATION_baselines.md) (the single-opponent GNINA
pre-registration, #63) first — this widens that one, it does not re-open it.

## The objection this closes

Run 2 showed the iron-approach geometric criterion (AUC 0.794) beats **AutoDock Vina's own
score** (AUC 0.471) on identical poses, and #63 set up **GNINA** as the one metal-aware learned
opponent. The remaining rebuttal:

> Beating Vina — and even one CNN — is not enough. Quantify the criterion against *several* of
> the tools a practitioner would actually reach for on a metalloenzyme, head-to-head, or the
> "geometry beats the scoring function" claim is undersupported.

Issue #83 answers exactly that: ≥ 2 independent metal-aware baselines, run on the identical
held-out poses, reported as one comparison table with CIs.

## The opponents (fixed; chosen in advance)

Both are drop-in **rescorers of the frozen run-2 poses** — no re-docking, no pose regeneration,
no per-molecule tuning — so only the scorer changes. Both are metal-aware and free/open, the
strongest readily-available such opponents; each is Linux/x86 native (see "run" below).

1. **GNINA — CNNaffinity** (McNutt et al. 2021), the learned metal-aware rescorer. `--score_only`
   with the default shipped CNN model, exactly as fixed in
   [PREREGISTRATION_baselines.md](PREREGISTRATION_baselines.md); its saved outputs
   (`work/gnina_scores/`) feed straight in.
2. **PLANTS — ChemPLP** (Korb et al. 2009), an ant-colony docking engine whose ChemPLP scoring
   function carries **explicit metal-acceptor terms**. `--mode rescore` over the same poses;
   the per-pose `TOTAL_SCORE` (lower = better) is the ranking value
   (`openafr.baselines.parse_plants_features`, saved to `work/plants_scores/`).

Considered and **excluded** (documented, not silently dropped): smina/Vinardo and AutoDock4Zn —
smina/Vinardo is metal-blind (already the class Vina represents) and AutoDock4Zn is
zinc-parameterised, not heme-Fe; both are weaker, less informative opponents for this metal
site, so neither is part of this pre-registration.

## The rankings compared (all on the same poses)

    C  (ours)     iron filter (N-Fe < 3.0 A), then best Vina affinity among survivors
                  — the run-2 PRIMARY criterion, reproduced on these poses
    G  (GNINA)    best CNNaffinity over ALL poses  — GNINA with no mechanistic prior
    P  (PLANTS)   best (lowest) ChemPLP TOTAL_SCORE over ALL poses

Each baseline row answers the practitioner's question *"should I just use this tool instead of
your whole pipeline?"* Metrics for every method: AUC, EF@5%, and a percentile bootstrap 95% CI
over the actives (10,000 resamples, seed 42), all from `openafr.scoring`. Any active with no
rankable pose/score is ranked **LAST**, never dropped — identical to every prior gate.

## Held-out set, protocol and endpoints (fixed)

- Held-out set: the run-2 set unchanged — 7 held-out azoles (`actives_holdout_final.smi`) +
  348 property-matched decoys, rigid 5TZ1 chain A + heme (`work/receptor_A.pdb`).
- Poses: the exact run-2 Vina poses (`work/screen2`), regenerated bit-for-bit from the frozen
  `scripts/screen.sh` under `protocol.yaml`
  (sha256 `6eca16b1be2e5e33ae42d4b1b7a63a0e9b891a5dd6cbfb19642dd3e9a149aba4`). Same box,
  exhaustiveness, seed, receptor.
- **Endpoints: AUC and EF@5%.** EF@1% is excluded — [RESULTS_reliability.md](RESULTS_reliability.md)
  established it is not robust under seed sampling.

## Pre-committed interpretation (fixed before seeing numbers)

Read per baseline on **AUC and EF@5%**, with bootstrap CIs qualifying every call. "Beats" =
strictly higher point estimate on **both** endpoints AND non-overlapping bootstrap CIs on AUC;
otherwise it is a **TIE** (overlap = tie, either way — the same rule #63 fixed).

- **C beats every baseline.** The mechanistic criterion beats multiple competent, metal-aware
  scorers, not just broken vanilla Vina. The run-2 framing stands and strengthens: "geometry
  beats the scoring function" survives contact with the strongest readily-available opponents.
- **C ties a baseline.** The free, interpretable geometric criterion *matches* an opaque
  learned/physics scorer at zero training cost and full transparency. Reported as a **feature,
  not a loss** — the geometry recovers what the tool finds.
- **A baseline beats C.** On both endpoints with separated CIs. Per pre-registration, retire
  "geometry beats the scoring function" as the headline and lead with the self-validating gate
  instead; the criterion remains a cheap, interpretable first-pass filter. Do NOT re-tune the
  criterion to recover the comparison.

## Run (cloud-gated; grading is not)

The scorer binaries are Linux/x86 native — the same wall the ERG11 re-caller and the GNINA
baseline (#63) hit — so the numeric run happens on a cloud Linux host and the graded,
hash-checked scoring stays reproducible on any host (mirrors `data/ligands/PROVENANCE.md`
discipline constraint #4):

    scripts/run_gnina_rescore.sh   work/screen2 work/gnina_scores  work/receptor.pdbqt
    scripts/run_plants_rescore.sh  work/screen2 work/plants_scores work/receptor_plants.mol2
    python scripts/rescore_baselines.py work/screen2 work/receptor_A.pdb   # runs anywhere

Until both output dirs exist the grader reports each missing baseline as
`not-yet-run (cloud-gated)` and headlines nothing — the honest current state. The exact tool
builds/versions are recorded by the run scripts (prereg reproducibility).

## Honest scope (fixed)

1. **Same held-out set, same rigid *C. albicans* receptor, presumed-inactive decoys** — the
   run-2 limitations carry over unchanged; this compares *scorers*, not datasets.
2. **Rescoring, not re-docking.** Every method ranks the identical Vina poses, so a baseline
   that would have found a better pose by re-docking is not credited that — a deliberate,
   single-variable comparison (only the scoring number changes), stated so no baseline is
   handicapped silently.
3. **This compares ranking criteria, not drugs.**
