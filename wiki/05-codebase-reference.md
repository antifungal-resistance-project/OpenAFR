# Codebase Reference

*A file-by-file map so you can find the right code fast. For the narrative of how the pieces
fit, read the pipeline deep dives ([track 1](03-drug-discovery-pipeline.md),
[track 2](04-early-warning-pipeline.md)) first; this page is the index you come back to.*

*The repo has grown a lot: track 1 now has a whole **candidate-confidence layer** (a stack of
per-candidate axes fused into an A/B/C dossier) and a **battery of pre-registered validation
gates** beyond the original run-2 holdout; track 2 now has the ERG11 re-caller **built and run**
plus a parallel **FKS1/echinocandin detection** track. This page is the current inventory.*

---

## The mental model

```
openafr/     ← the shared LIBRARY: deterministic, stdlib-or-rdkit-only, unit-tested. The core.
scripts/     ← the runnable ENTRY POINTS (CLIs + shell). Thin glue over openafr/, plus the
               un-reproducible tool/network orchestration kept OUT of the tested core.
tests/       ← known-answer tests locking the library's math, parsers, and orchestration logic.
data/        ← inputs: structures, ligands, ortholog + reference sequences, the snapshot store.
work/        ← run outputs + every write-up (RESULTS_*.md) and pre-registration (PREREGISTRATION_*.md).
```

Rule of thumb: **logic that must be correct lives in `openafr/` and has a test; a `scripts/`
file mostly parses arguments and calls into `openafr/`.** The exceptions are the parts that
shell out to external binaries or the network (docking, the reads-to-consensus half of both
re-callers, the network precedent/synth probes) — those live only in `scripts/` (or behind an
injectable `opener` in the library) and are honestly labelled as untested-in-CI.

---

## `openafr/` — the library (28 modules)

### Shared core

| File | What it is |
|---|---|
| `__init__.py` | Package docstring; explains why the shared parsers/scoring were consolidated here (previously copy-pasted across scripts). |
| `pdbqt.py` | The canonical parsers: `read_poses`, `read_scores`, `iron_position`, and `min_nitrogen_iron_distance` — the **N–Fe measurement the whole track-1 result rests on**. Heavy atoms only. |
| `scoring.py` | The ranking math: `metrics()` (AUC / EF / BEDROC via `rdkit.ML.Scoring`) and `report()` (sorts, places unrankable molecules **last**, prints). |
| `protocol.py` | Loads + hash-verifies `protocol.yaml`. Tiny YAML-subset parser (no PyYAML dep). `verify()`, `--freeze`, `--shell`. |

### Track 1 — geometry criterion, its stress-tests, and the confidence axes

| File | What it is |
|---|---|
| `coordinator.py` | **Which** nitrogen reaches the iron. Classifies a pose's coordinating N by bonded-neighbour topology (ring vs. nitrile vs. amine/amide) — the "is this rank azole-mechanism or an out-of-mechanism artifact?" check. Also holds `min_aromatic_nitrogen_iron_distance` (the tested-but-not-adopted in-mechanism twin of the criterion). |
| `ensemble.py` | Aggregates the mode-C N–Fe distance over multiple receptor conformers (ensemble-min primary, ensemble-mean diagnostic). Pure functions over already-computed distances. |
| `pose_convergence.py` | Did one docking search settle on one geometry or scatter? Unaligned index-matched pose RMSD clustering; thresholds frozen in a pre-registration. |
| `conservation.py` | Cross-species CYP51 **pocket conservation** map (C. albicans ↔ C. auris ↔ A. fumigatus) from pinned orthologs via a deterministic Needleman–Wunsch alignment. Same pocket definition as `mapping.py`. |
| `inactives.py` | **Verified-inactive** decoy construction: compounds *measured* not to inhibit C. albicans (single-variable swap vs. the presumed-inactive decoys). Encodes the CLSI-anchored evidence rule. |
| `potency.py` | Continuous-potency correlation core: does tighter N–Fe track *higher measured potency*? `median_pchembl` + `spearman`/`permutation_p`/`bootstrap_ci` (numpy only). |
| `baselines.py` | Multi-baseline comparison scaffold: turns each rival scorer's saved score files into per-pose ranking keys so the criterion is quantified against real alternatives, not just vanilla Vina. |
| `gnina.py` | Parser for GNINA `--score_only` output — the metal-aware learned baseline scorer. Pure/stdlib; the binary runs off-box. |
| `admet.py` | ADMET / physicochemical + structural-alert profiling (Ro5, Veber, PAINS/BRENK/NIH, ESOL, a coarse hERG risk-factor hint). All informational, never a gate. Local RDKit, no network. |
| `synth.py` | Synthesizability / availability: SAscore (Ertl, RDKit-Contrib, local) + an optional best-effort PubChem presence probe behind an injectable `opener`. |
| `cyp_offtarget.py` | Human drug-metabolizing-CYP off-target liability (3A4/2C9/2D6/2C19/1A2): the type-II heme-ligation warhead the azole mechanism *rewards on target* is the liability a chemist *dreads* on the human panel. Coarse mechanistic + pharmacophore flags, not IC50s. |
| `precedent.py` | Assay precedent — has this candidate already been measured against CYP51 (and failed)? Merges ChEMBL + PubChem + BindingDB; `no-precedent` ≠ untested. Network behind an injectable `opener`. |
| `filedrawer.py` | Executable file-drawer model: bounds the *leak* L (fraction of existing negative measurements that evade the queried stores) as a band, and shows adding sources shrinks it. Makes the `no-precedent` caveat quantitative. |
| `dossier.py` | **Fuses** a molecule's already-computed evidence (geometry + precedent + ADMET + CYP + synth + coordinator identity) into one machine-readable dossier and assigns an **A/B/C** priority from a documented table. Computes nothing new; never re-ranks; never drops a row. |

### Track 2 — genomic early-warning surveillance

| File | Track | What it is |
|---|---|---|
| `earlywarning.py` | 2 | Snapshot store: NCBI pull → normalized byte-stable per-release snapshot; `baseline_as_of`, `diff_snapshots`, `update_index`. Now carries `fks1_call` alongside `erg11_call`. Stdlib only. |
| `emergence.py` | 2 | Emergence detection: `detect_emergence(call_field=…)` flags first-appearing / rising substitutions over *called* isolates only, with the backlog guard — now gene-parameterized (serves ERG11 **and** FKS1). |
| `mapping.py` | 2 | Mutation → pocket mapping: places a residue in the modeled pocket, measures drug/iron distance, decides buildability. |
| `structural.py` | 2 | Structural so-what: `estimate_fit_consequence` → measured / estimate / no-call fit verdict. Holds `EMPIRICAL_DOCKS` (Y132F). *ERG11/azole only — deliberately not extended to FKS1.* |
| `alert.py` | 2 | Alert composition: joins emergence + structural, assigns ACT-NOW / WATCH / CONTEXT, renders Markdown + JSON. Adds `compose_fks1_alerts` / `render_fks1_markdown` (**WATCH/CONTEXT only, never ACT-NOW** — detection-only). |
| `delivery.py` | 2 | Delivery decision: deliver only on genuine new news; per-gene `STATE.json` memory. Routes off the result's `gene` marker so ERG11 and FKS1 deliver independently. |
| `backtest.py` | 2 | The credibility test: H1-arrival / H2-null / H3-mechanism, plus `panel_prevalence` (azole/echinocandin event frequency) with Wilson CIs. |
| `recaller.py` | 2 | **ERG11 re-caller core:** consensus CDS → substitution tokens. Deterministic, hash-pinned B8441 reference, four honesty constraints. (The reads→consensus half lives in `scripts/recall_erg11.py`.) |
| `fks1_caller.py` | 2 | **FKS1 re-caller core** (v2): *windowed* consensus → echinocandin-resistance tokens (HS1 ≈ S639, HS2 ≈ R1354). Detection-only; same honesty bar as `recaller.py`, applied per hot-spot window. |
| `runlog.py` | 2 | Append-only provenance log for the re-callers' un-reproducible tool/network runs. Best-effort; never breaks the run. |

## `scripts/` — the entry points

### Track 1 — the discovery path (front door + filter → prep → screen → rank → report)

| Script | What it does |
|---|---|
| `applicability_triage.py` | **The front door, step 1.** Per-molecule verdict — does the validated method even *cover* this molecule? (`in-envelope-novel` / `near-known` / `out-of-domain` / `out-of-scope` / `unparseable`). Runs anywhere, no docking. |
| `score_molecule.py` | **The front door, full.** SMILES → triage → prep → dock → N–Fe geometry, one command, ~9 s/molecule on a laptop. `--triage-only`, `--json`. Never fakes a number. |
| `prep_receptor.py` | `data/structures/5TZ1.pdb` → `work/receptor.pdbqt` (chain A + heme); asserts atoms match the tracked `work/receptor_A.pdb`. |
| `prep_ligands.py` | SMILES → 3D structures for docking. Deterministic; failures logged, not dropped. |
| `make_decoys.py` | Property-matched **presumed-inactive** decoys for the known azoles (defends the circular-decoy and size-shortcut traps). |
| `filter_library.py` | Enforce the pre-registered applicability domain on a novel library (> 45 heavy atoms / no aromatic N excluded, with reasons). |
| `screen.sh` | **The docking engine.** One frozen protocol for every molecule, all cores, resumable. Reads the box from `protocol.yaml`. `RECEPTOR=` to dock a mutant/ensemble/human receptor. |
| `run_screen2.sh` | Thin wrapper: run-2 held-out **validation** screen through `screen.sh`. |
| `screen_consensus.sh` | Dock every ligand under R independent Vina searches (the multi-seed reliability experiment). |
| `rank_candidates.py` | Screened library → ranked candidate list by the validated N–Fe criterion. Unrankable go last. |
| `report_candidates.py` | Ranked list → human-readable Markdown report with per-candidate evidence + provenance. |
| `check_precedent.py` | Flag candidates ChEMBL / PubChem / BindingDB has already measured against the target (the file-drawer guard). |
| `build_dossier.py` | Screened+ranked library → one machine-readable **dossier** per molecule (A/B/C), fusing every confidence axis. |

### Track 1 — validation gates (each implements one pre-registration; exits nonzero on fail)

| Script | Grades |
|---|---|
| `validate_gate2.py` | **THE gate.** The run-2 held-out N–Fe criterion vs the frozen pass bar (AUC ≥ 0.70 AND EF@1% ≥ 5.0×). |
| `validate_gate.py` | The original run-1 validation gate. |
| `validate_gate_reliability.py` | Multi-seed reliability (the check that found EF@1% was inflated). |
| `validate_gate_verified.py` | The criterion vs **measured** inactives (look #6 — closes the presumed-inactive caveat; AUC 0.716). |
| `validate_gate_external.py` | A fresh, never-scored external active/decoy holdout (#82). |
| `validate_gate_temporal.py` | A publication-date split — does the frozen criterion generalize to the newest-literature azoles (#82; AUC 0.754)? |
| `validate_gate_potency.py` | Continuous potency: does tighter N–Fe track higher pChEMBL (#81, a pre-registered null)? |
| `validate_gate_channel.py` | Mode F (channel engagement), the one permitted orthogonal within-azole feature (look #7 — FAIL). |
| `validate_gate_ensemble.py` / `validate_gate_ensemble_confirm.py` | Within-azole ranking on the crystallographic ensemble (look #8 PASS) and its independent-active confirmation (look #9 FAIL — closes the line). |
| `validate_gate_aromatic.py` | Does restricting the coordinator to an **aromatic** ring N preserve the gate? (Tested → **lowers** AUC to 0.729; criterion left unchanged.) |
| `validate_gate_auris.py` | The *C. auris* Y132F resistant-pocket port gate. |
| `validate_gate_axial.py` / `validate_gate_consensus.py` | Axial-geometry and consensus-docking pre-registered gates. |

### Track 1 — datasets, baselines, receptors, and the heme audit

| Script | What it does |
|---|---|
| `make_verified_inactives.py` | Measured-inactive set from ChEMBL bioactivity (consistent inactivity, incl. the CYP51 enzyme target). |
| `make_verified_actives.py` | Independent measured-active azole set — the mirror of verified inactives (for look #9). |
| `make_external_holdout.py` / `make_temporal_holdout.py` | The fresh never-scored and the publication-date-split holdouts (#82). |
| `make_potency_set.py` | CYP51 ligands with measured pChEMBL potency (#81). |
| `domain_lifted_pool_audit.py` / `domain_sensitivity.py` | Stage-0 pool size if the domain ceiling is lifted; sensitivity of the domain thresholds (#87). |
| `preflight_baselines.py`, `rescore_baselines.py`, `rescore_gnina.py` | Offline preflight + gates for the multi-baseline / GNINA metal-aware comparison (#63, #83). |
| `run_baselines_all.sh`, `run_gnina_rescore.sh`, `run_plants_rescore.sh` | The cloud-gated runners for the rival scorer binaries (Linux/x86). |
| `prep_receptor_ensemble.py` | Extra C. albicans CYP51 conformers → aligned docking receptors. |
| `prep_receptor_human.py` / `screen_human.sh` | Build + screen the **human** CYP51 counter-screen receptor (3LD6). |
| `human_positive_control.py` | Human-CYP51 native-pose-recovery control (#73). |
| `mutate_receptor.py` | `receptor_A.pdb` + a **deletion** mutation (Y132F, V125A) → a mutant heavy-atom receptor. |
| `repack_mutant.py` | `receptor_A.pdb` + an **add-atom** mutation (K143R, F126L) → a mutant receptor, side chain rebuilt by a pinned PyMOL repacker. |
| `pocket_conservation.py` | Cross-species CYP51 pocket-conservation map (#76). |
| `audit_heme_model.py` | Audit the heme/iron model the whole criterion rests on (#72). |

### Track 1 — per-candidate confidence axes (fed into the scorecard/dossier)

| Script | Axis |
|---|---|
| `scorecard.py` | Pure decision logic for the consolidated per-candidate confidence scorecard (#88). |
| `shortlist_confidence.py` | Seed-stability (#68) + human-CYP51 selectivity (#73). |
| `coordinator_identity.py` | Is a candidate's iron-reach azole-like or an artifact (nitrile/amine)? (#120). |
| `admet_filter.py`, `synth_filter.py`, `cyp_offtarget_screen.py` | Annotate a library with ADMET (#75), synthesizability (#86), human-CYP liability (#74). |
| `candidate_stereo.py` / `enumerate_stereo.py` | Stereochemistry axis + its enumeration/prep (#85). |
| `candidate_protomer.py` / `enumerate_protomers.py` | Protonation/tautomer (pH-7.4 microstate) axis + prep (#84). |
| `candidate_ensemble.py` | Per-candidate ensemble-consistency across the CYP51 conformers (#70). |
| `pose_convergence.py` | Per-candidate within-search pose convergence (#71). |
| `candidate_resistance.py` / `candidate_resistance_addatom.py` | C. auris resistance-pocket retention for deletion and add-atom mutants (#69, #77). |

### Track 1 — pose analysis, visualization, and the drift harness

| Script | What it does |
|---|---|
| `analyze_poses.py`, `pose_analysis.py`, `make_pose_view.py`, `render_views.sh`, `metal_filter_demo.py` | Metal-coordination analysis and pose visualization helpers. |
| `repro.py` | **One-command reproducibility/drift harness** (#89): verifies `protocol.yaml`, every `work/PREREG_*.sha256`, and that the consolidated scorecard re-derives byte-for-byte. Runs on every PR. |

### Track 2 — early-warning surveillance

| Script | What it does |
|---|---|
| `probe_ncbi_auris.py` | The original spike: *is NCBI Pathogen Detection a usable feed?* (Answer: metadata yes, resistance calls no.) |
| `snapshot_ncbi_auris.py` | Snapshot store CLI: `pull` a release into a versioned snapshot; baseline/diff. |
| `detect_emergence.py` | Emergence-detection CLI (`demo` mode); gene-selectable (ERG11/FKS1). |
| `map_mutation.py` | Mutation → pocket mapping CLI. |
| `structural_sowhat.py` | Structural so-what CLI (`estimate <tokens…>`). ERG11/azole only. |
| `compose_alert.py` | Alert-composition CLI. `--gene fks1` for the detection-only FKS1 surface. |
| `deliver_alert.py` | Delivery CLI (`demo`, `--markdown`, `--gene {erg11,fks1}` with its own digest/state). |
| `backtest_earlywarning.py` | Backtest CLI: `mechanism`, `replay`, `prevalence`, plus the H1 arrival budget. |
| `recall_erg11.py` | ERG11 re-caller CLI. `call` = deterministic. `plan` = offline preflight. `recall`/`fill` = reads→consensus orchestration (external binaries). **Run: n=199 representative fill, 80.4% azole event frequency.** |
| `recall_fks1.py` | FKS1 re-caller CLI (v2): windowed `call`/`recall`/`fill`/`prevalence`, tool-gated like ERG11. Prevalence still "NOT MEASURED YET" (awaiting a real fill). |
| `recaller_sanity.py` / `recaller_sanity_fks1.py` | Sanity-run each re-caller against a fixed published known-genotype truth set (positive controls). |
| `cloud_recaller_setup.sh` | Provision a fresh Linux/x86_64 host with the bioconda toolchain for the re-callers' real runs. |
| `README_recaller_sanity.md` | Notes for the sanity run. |

## `tests/`

Every test is a **known-answer** test — it locks a specific numeric or structural output so a
change that alters behavior fails loudly. 60 test files. Highlights, by area:

| Area | Tests |
|---|---|
| Criterion + parsers + protocol | `test_gate_metrics.py`, `test_pdbqt.py`, `test_protocol.py`, `test_prep_receptor.py`, `test_mutate_receptor.py`, `test_repack_mutant.py` |
| Discovery path | `test_report_candidates.py`, `test_filter_library.py`, `test_applicability_triage.py`, `test_score_molecule.py` |
| Validation experiments | `test_reliability_gate.py`, `test_axial_gate.py`, `test_consensus_gate.py`, `test_verified_gate.py`, `test_external_gate.py`, `test_temporal_gate.py`, `test_potency.py`, `test_channel_gate.py`, `test_ensemble.py` |
| Datasets + baselines + audit | `test_inactives.py`, `test_baselines.py`, `test_gnina.py`, `test_preflight_baselines.py`, `test_heme_audit.py`, `test_domain_sensitivity.py`, `test_pocket_conservation.py` |
| Confidence axes | `test_admet.py`, `test_synth.py`, `test_cyp_offtarget.py`, `test_precedent.py`, `test_check_precedent.py`, `test_filedrawer.py`, `test_candidate_stereo.py`, `test_candidate_protomer.py`, `test_candidate_ensemble.py`, `test_pose_convergence.py`, `test_candidate_resistance.py`, `test_scorecard.py`, `test_shortlist_confidence.py`, `test_coordinator.py`, `test_dossier.py`, `test_repro.py` |
| Track 2 | `test_earlywarning_snapshot.py`, `test_emergence.py`, `test_mapping.py`, `test_structural.py`, `test_alert.py`, `test_delivery.py`, `test_backtest.py`, `test_recaller.py`, `test_recall_orchestration.py`, `test_recall_preflight.py`, `test_recall_fill_logging.py`, `test_runlog.py` |
| FKS1 (v2) | `test_fks1_caller.py`, `test_fks1_alert.py`, `test_fks1_delivery.py`, `test_fks1_prevalence.py`, `test_recall_fks1_orchestration.py` |

CI runs three workflows (`.github/workflows/`): `tests.yml` (with a coverage gate — the openafr
core sits at ~97%), `repro.yml` (the drift harness), and `earlywarning.yml` (the scheduled
surveillance run).

## `data/`

| Path | Contents |
|---|---|
| `data/structures/5TZ1.pdb` | The real *C. albicans* CYP51 crystal structure (+ heme + bound azole) from RCSB — the docking receptor's source. |
| `data/ligands/actives.smi` | 8 known azoles used to form the hypothesis (training). |
| `data/ligands/actives_holdout_final.smi` | 7 held-out azoles — the blind test set. |
| `data/ligands/decoys*.smi` | Property-matched (presumed-inactive) decoys (~50 per active). |
| `data/ligands/verified_inactives.smi` | Compounds **measured** not to inhibit C. albicans (n=279). |
| `data/ligands/temporal_actives.smi`, `temporal_decoys.smi` | The publication-date holdout set (#82); `PROVENANCE.md` records the split. |
| `data/cyp51_orthologs/` | Pinned CYP51 ortholog sequences (C. albicans / C. auris / A. fumigatus A+B) for the pocket-conservation map. |
| `data/earlywarning/erg11_reference/` | The pinned B8441 ERG11 reference CDS (XM_085597798.1) + `PROVENANCE.md`. |
| `data/earlywarning/fks1_reference/` | The pinned B8441 GSC1/FKS1 reference CDS (XM_085597048.1) + `PROVENANCE.md` (numbering verified by translation). |
| `data/earlywarning/snapshots/` | The versioned NCBI snapshot store (big TSVs gitignored; `INDEX.tsv` committed). Has its own `README.md`. |
| `data/earlywarning/digest/` | The delivery digests + `STATE.json` / `STATE_FKS1.json`. Has its own `README.md`. |
| `data/earlywarning/runlog/` | The append-only run-provenance `.jsonl` files (incl. `fill.jsonl`). |
| `data/earlywarning/recaller_sanity/` | The fixed known-genotype truth sets for the sanity runs (ERG11 + FKS1). |

## `work/` — outputs, results, and pre-registrations

`work/` holds run outputs (much of it gitignored) plus the committed documents central to the
project's credibility:

- **`PREREGISTRATION_*.md`** (+ `PREREG_*.sha256`) — the rules and pass bars, written and
  SHA-256-frozen *before* data was generated. ~30 of them now, from `run2` through the
  confidence axes (`coordinator`, `pose_convergence`, `candidate_*`), the holdouts (`external`,
  `temporal`, `potency`), the baselines, and the aromatic-criterion test.
- **`RESULTS_*.md`** — the write-up of every run, *including the ones that failed or that
  qualified an earlier claim*. Notable reads:
  - `RESULTS_run2_holdout.md` — the headline result; `RESULTS_reliability.md` — the self-inflicted EF@1% caveat.
  - `RESULTS_verified_inactives.md` — measured inactives (largest caveat closed) + the located ceiling.
  - `RESULTS_channel_engagement.md`, `RESULTS_ensemble.md`, `RESULTS_ensemble_confirm.md` — the within-azole line, closed.
  - `RESULTS_temporal.md`, `RESULTS_external.md` — generalization to never-scored / newest-literature data.
  - `RESULTS_scorecard.md`, `RESULTS_dossier.md`, `RESULTS_coordinator.md`, `RESULTS_aromatic_criterion.md` — the confidence layer and the coordinator-identity finding.
  - `RESULTS_candidate_shortlist.md` — the worked repurposing shortlist; `RESULTS_repro.md` — the drift harness.
  - `RESULTS_auris_Y132F.md` — the only `measured` structural verdict (track 2).
  - `RESULTS_prevalence.md` — **the measured azole event frequency (80.4%)**, the number track 2 was gated on; `RESULTS_backtest.md` — the pre-registered credibility test.
  - `RESULTS_earlywarning_spike.md` — the "NCBI has no resistance calls" discovery.
- **Runbooks / methods / spec:** `RUNBOOK_recaller_run.md`, `RUNBOOK_fks1_run.md`,
  `RUNBOOK_baselines_multi.md`, `METHODS_file_drawer.md`, `SPEC_candidate_dossier.md`,
  `PREPRINT_geometry_ceiling.md`. The per-candidate dossiers live under `work/candidates/`.

## Root files

| File | What it is |
|---|---|
| `README.md` | What the pipeline does + how to run it. |
| `VISION.md` | Strategy, novelty, where it goes next. |
| `TODOS.md` | Deferred work — now headed by the coordinator-identity finding and the FKS1 real-fill run. |
| `protocol.yaml` | The frozen docking box + pass bar. Hash-pinned in `openafr/protocol.py`. |
| `environment.yml` | The version-pinned conda toolchain (rdkit, openbabel, vina, pytest, pytest-cov, numpy). |
| `pytest.ini` | Test configuration. |
| `LICENSE` | PolyForm Noncommercial 1.0.0. |

---

Recurring patterns that show up across many of these files are collected in the
[Project Overview](02-project-overview.md#the-recurring-design-patterns-learn-these-once).
Stuck on a term? → **[Glossary](06-glossary.md)**. Ready to build? → **[Contributing](07-contributing.md)**.
