# Runbook — the multi-baseline comparison run (#83)

The operational companion to `work/PREREGISTRATION_baselines_multi.md`. The harness,
grader and pre-registration are merged (#109); what remains is executing the two
metal-aware baseline rescores on a capable host and recording the filled table. This is
the order to run it in, and how to read each step.

**Run on a Linux/x86 host.** GNINA is a compiled CNN binary and PLANTS/SPORES are
closed-source Linux/x86 executables — the same osx-arm64 wall the ERG11/FKS1 re-callers
hit. The grader (`rescore_baselines.py`) and the preflight run anywhere; only the two
scorers are host-bound.

**One human-in-the-loop dependency you cannot script around:** PLANTS and SPORES are a
free *academic* download behind registration at the University of Tübingen
(<http://www.tcd.uni-konstanz.de/plants_download/>). They are not on conda/pip and cannot
be fetched by any container build. Download them yourself once and put them on `PATH`.

## 0. Prerequisites (once, on the run host)

```bash
git checkout main && git pull                 # must include #109 (harness) + this runbook
conda env create -f environment.yml           # or: mamba env create -f environment.yml
conda activate openafr
python -m pytest -q                            # sanity: the tested core is green here
```

Tools + pinned versions (the run scripts capture each tool's exact `--version` into a
`*_VERSION.txt` beside its outputs — that captured string is the frozen provenance the
pre-registration refers to; pin at least these):

| tool | version | source | used for |
|------|---------|--------|----------|
| AutoDock Vina | 1.2.7 | conda-forge `vina` | regenerate the frozen run-2 poses (step 1) |
| Open Babel | 3.1.1 | conda-forge `openbabel` | SMILES→pdbqt, pose→mol2 for PLANTS |
| GNINA | 1.3 | github.com/gnina/gnina (release binary or image) | baseline G (CNNaffinity) |
| PLANTS | 1.2 (`PLANTS1.2_64bit`) | Tübingen academic download | baseline P (ChemPLP) |
| SPORES | 1.3 | Tübingen academic download | prepare the PLANTS receptor mol2 |

- **Disk:** the ~355-ligand pose tree is ~40 MB; GNINA/PLANTS outputs are small text.
  No multi-GB downloads here (unlike the re-caller).
- **GPU is optional.** GNINA `--score_only` runs on CPU; a GPU only speeds it up.

## 1. Regenerate the frozen run-2 poses (the identical poses every scorer rescores)

The ~40 MB Vina pose tree is gitignored and regenerates bit-for-bit from the frozen
protocol. This is the standard run-2 reproduction (see `work/RESULTS_run2_holdout.md`):

```bash
python scripts/prep_receptor.py                 # -> work/receptor.pdbqt (also the GNINA receptor)
python scripts/prep_ligands.py <run-2 ligands>  # -> work/run2_ligands/*.pdbqt  (per RESULTS_run2)
scripts/run_screen2.sh                          # docks -> work/screen2/*.pdbqt under protocol.yaml
```

The protocol is hash-frozen (`protocol.yaml` sha `6eca16b1…`); do not tune it. The
preflight (step 3) re-checks that sha, so a drifted protocol is caught before scoring.

## 2. Prepare the two baseline receptors

Both rescore the same poses but each scorer wants the receptor in its own format:

```bash
# GNINA: a pdbqt receptor — the same one docking used.
#   work/receptor.pdbqt is produced by prep_receptor.py in step 1.

# PLANTS: a SPORES-prepared mol2, protonated + typed for ChemPLP.
SPORES --mode complete work/receptor_A.pdb work/receptor_plants.mol2
```

`work/receptor_A.pdb` is the committed rigid 5TZ1 chain A + heme (the grader also reads
the heme Fe from it for mode C). Keep the PLANTS binding site on the frozen box centre
(70.61, 66.28, 4.18); `run_plants_rescore.sh` already sets that.

## 3. Preflight — go/no-go before spending time (runs anywhere, no scorers)

```bash
python scripts/preflight_baselines.py
```

- Verifies the two frozen hashes (pre-registration + `protocol.yaml`). A **DRIFT** verdict
  (exit 1) means the rules or protocol moved after freezing — stop and resolve; any number
  graded against a drifted frozen surface is not credible.
- Reports poses, the grading receptor, and per-baseline: run-state, receptor presence, and
  which tools are on `PATH`. **RUN-STEPS-REMAIN** (exit 2) is the normal pre-run state;
  **READY-TO-GRADE** (exit 0) means both baselines are already scored.

## 4. Run — one command

```bash
scripts/run_baselines_all.sh
```

Runs the preflight, then the GNINA rescore, then the PLANTS rescore, then the grader —
printing the comparison table. Each rescore self-gates on its binary, so on a host missing
one scorer the wrapper skips it (not fatal) and the grader marks it `not-yet-run`. Both
rescores are resumable (a ligand whose output already exists is skipped), so a re-run after
an interruption is cheap.

Equivalent explicit steps, if you prefer to run them one at a time:

```bash
scripts/run_gnina_rescore.sh  work/screen2 work/gnina_scores  work/receptor.pdbqt
scripts/run_plants_rescore.sh work/screen2 work/plants_scores work/receptor_plants.mol2
python scripts/rescore_baselines.py work/screen2 work/receptor_A.pdb
```

## 5. Read + record the result

The grader prints one row per method — mode C plus each baseline — with AUC, EF@5%, and a
bootstrap 95% CI, on the identical poses. Apply the **pre-registered interpretation**
(`work/PREREGISTRATION_baselines_multi.md`), which was fixed before any number was seen:

- **"Beats"** = strictly higher on BOTH AUC and EF@5% AND non-overlapping AUC bootstrap
  CIs. Otherwise it is a **TIE** (CI overlap = tie, either direction).
- C beats every baseline → "geometry beats the scoring function" survives the strongest
  readily-available metal-aware opponents. C ties → reported as a feature (the free,
  interpretable criterion matches an opaque scorer). A baseline beats C → retire that
  headline and lead with the self-validating gate; **do not re-tune C to recover it.**

Then record the table in **`work/RESULTS_baselines_multi.md`** (create it), commit the
`*_VERSION.txt` provenance files alongside, and close #83 with the table. The pose trees
and per-ligand score files stay gitignored; the results markdown + version files are the
durable record, exactly as every other `RESULTS_*.md`.

## Honest scope (carried from the pre-registration)

Same held-out set, same rigid *C. albicans* receptor, presumed-inactive decoys — this
compares **scorers on identical poses**, not datasets and not re-docking. A baseline that
would have found a better pose by re-docking is deliberately not credited that; only the
scoring number changes. See `work/PREREGISTRATION_baselines_multi.md` §"Honest scope".
