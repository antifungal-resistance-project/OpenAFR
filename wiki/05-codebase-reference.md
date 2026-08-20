# Codebase Reference

*A file-by-file map so you can find the right code fast. For the narrative of how the pieces
fit, read the pipeline deep dives ([track 1](03-drug-discovery-pipeline.md),
[track 2](04-early-warning-pipeline.md)) first; this page is the index you come back to.*

---

## The mental model

```
openafr/     ← the shared LIBRARY: deterministic, stdlib-or-rdkit-only, unit-tested. The core.
scripts/     ← the runnable ENTRY POINTS (CLIs + shell). Thin glue over openafr/, plus the
               un-reproducible tool/network orchestration kept OUT of the tested core.
tests/       ← known-answer tests locking the library's math, parsers, and orchestration logic.
data/        ← inputs: structures, ligands, the ERG11 reference, the snapshot store.
work/        ← run outputs + every write-up (RESULTS_*.md) and pre-registration (PREREGISTRATION_*.md).
```

Rule of thumb: **logic that must be correct lives in `openafr/` and has a test; a `scripts/`
file mostly parses arguments and calls into `openafr/`.** The exceptions are the parts that
shell out to external binaries or the network (docking, the reads-to-consensus half of the
re-caller) — those live only in `scripts/` and are honestly labelled as untested-in-CI.

---

## `openafr/` — the library

| File | Track | What it is |
|---|---|---|
| `__init__.py` | shared | Package docstring; explains why the shared parsers/scoring were consolidated here (they were previously copy-pasted across scripts). |
| `pdbqt.py` | 1 | The canonical parsers: `read_poses`, `read_scores`, `iron_position`, and `min_nitrogen_iron_distance` — the N–Fe measurement the whole track-1 result rests on. Heavy atoms only. |
| `scoring.py` | 1 | The ranking math: `metrics()` (AUC / EF / BEDROC via `rdkit.ML.Scoring`) and `report()` (sorts, places unrankable molecules **last**, prints). |
| `protocol.py` | 1 | Loads + hash-verifies `protocol.yaml`. Tiny YAML-subset parser (no PyYAML dep). `verify()`, `--freeze`, `--shell`. |
| `earlywarning.py` | 2 | Snapshot store: NCBI pull → normalized byte-stable per-release snapshot; `baseline_as_of`, `diff_snapshots`, `update_index`. Stdlib only. |
| `emergence.py` | 2 | Emergence detection: `detect_emergence` flags first-appearing / rising ERG11 substitutions over *called* isolates only, with the backlog guard. |
| `mapping.py` | 2 | Mutation → pocket mapping: `map_mutation` places a residue in the modeled pocket, measures drug/iron distance, decides buildability. |
| `structural.py` | 2 | Structural so-what: `estimate_fit_consequence` → measured / estimate / no-call fit verdict. Holds `EMPIRICAL_DOCKS` (Y132F). |
| `alert.py` | 2 | Alert composition: `compose_alerts` joins emergence + structural, assigns ACT-NOW / WATCH / CONTEXT priority, renders Markdown + JSON. |
| `delivery.py` | 2 | Delivery decision: `new_alerts` / `is_noop` / `next_state` — deliver only on genuine new news; `STATE.json` memory. |
| `backtest.py` | 2 | The credibility test: H1-arrival / H2-null / H3-mechanism, plus `panel_prevalence` (azole event frequency) with Wilson CIs. |
| `recaller.py` | 2 | ERG11 re-caller **core** (the missing piece): consensus CDS → substitution tokens. Deterministic, hash-pinned reference, four honesty constraints. |
| `runlog.py` | 2 | Append-only provenance log for the re-caller's un-reproducible tool/network runs. Best-effort; never breaks the run. |

## `scripts/` — the entry points

### Track 1 — drug discovery

| Script | What it does |
|---|---|
| `prep_receptor.py` | `data/structures/5TZ1.pdb` → `work/receptor.pdbqt` (chain A + heme); asserts atoms match the tracked `work/receptor_A.pdb`. |
| `prep_ligands.py` | SMILES → 3D structures for docking. Deterministic; failures logged, not dropped. |
| `make_decoys.py` | Generate property-matched decoys for the known azoles (defends the circular-decoy and size-shortcut traps). |
| `filter_library.py` | Enforce the pre-registered applicability domain on a novel library (> 45 heavy atoms / no aromatic N excluded, with reasons). |
| `screen.sh` | **The docking engine.** One frozen protocol for every molecule, all cores, resumable. Reads the box from `protocol.yaml`. |
| `run_screen2.sh` | Thin wrapper: runs the run-2 held-out **validation** screen through `screen.sh`. |
| `screen_consensus.sh` | Dock every ligand under R independent Vina searches (the multi-seed reliability experiment). |
| `validate_gate2.py` | **THE gate.** Grades run-2 by the N–Fe criterion against the frozen pass bar; exits nonzero on fail. |
| `validate_gate.py` | The original run-1 validation gate. |
| `validate_gate_reliability.py` | Multi-seed reliability gate (the check that found EF@1% was inflated). |
| `validate_gate_auris.py` | The *C. auris* Y132F resistant-pocket port gate. |
| `validate_gate_axial.py` | Axial-geometry pre-registered gate. |
| `validate_gate_consensus.py` | Consensus-docking pre-registered gate. |
| `rank_candidates.py` | Screened novel library → ranked candidate list (by the validated N–Fe criterion). Unrankable go last. |
| `report_candidates.py` | Ranked list → human-readable Markdown report with per-candidate evidence + provenance. |
| `mutate_receptor.py` | `work/receptor_A.pdb` + a deletion mutation (Y132F, V125A) → a mutant heavy-atom receptor PDB. |
| `analyze_poses.py`, `pose_analysis.py`, `make_pose_view.py`, `render_views.sh`, `metal_filter_demo.py` | Metal-coordination analysis and pose visualization helpers. |

### Track 2 — early-warning surveillance

| Script | What it does |
|---|---|
| `probe_ncbi_auris.py` | The original spike: *is NCBI Pathogen Detection a usable feed?* (Answer: metadata yes, resistance calls no.) |
| `snapshot_ncbi_auris.py` | Snapshot store CLI: `pull` a release into a versioned snapshot; baseline/diff. |
| `detect_emergence.py` | Emergence-detection CLI (has a `demo` mode). |
| `map_mutation.py` | Mutation → pocket mapping CLI. |
| `structural_sowhat.py` | Structural so-what CLI (`estimate <tokens…>`). |
| `compose_alert.py` | Alert-composition CLI. |
| `deliver_alert.py` | Delivery CLI (`demo`, `--markdown`). |
| `backtest_earlywarning.py` | Backtest CLI: `mechanism`, `replay`, `prevalence`, plus the H1 arrival budget. |
| `recall_erg11.py` | ERG11 re-caller CLI. `call` = deterministic (no tools). `plan` = offline preflight. `recall`/`fill` = the reads→consensus orchestration (needs external binaries). |
| `recaller_sanity.py` | Sanity-run the re-caller against a fixed published known-genotype truth set (the wild-type + known-mutant positive controls). |
| `cloud_recaller_setup.sh` | Provision a fresh Linux/x86_64 host with the bioconda toolchain for the re-caller's real run. |
| `README_recaller_sanity.md` | Notes for the sanity run. |

## `tests/`

Every test is a **known-answer** test — it locks a specific numeric or structural output so a
change that alters behavior fails loudly. Highlights:

| Test | Locks |
|---|---|
| `test_gate_metrics.py` | The AUC / EF / BEDROC math (`openafr.scoring`). |
| `test_pdbqt.py` | The pose/score/iron parsers. |
| `test_protocol.py` | The protocol loader + hash verification. |
| `test_prep_receptor.py`, `test_mutate_receptor.py` | Receptor prep + deletion mutations. |
| `test_report_candidates.py`, `test_filter_library.py` | The discovery-path report + applicability filter. |
| `test_reliability_gate.py`, `test_axial_gate.py`, `test_consensus_gate.py` | The pre-registered experiment gates. |
| `test_earlywarning_snapshot.py` | Snapshot normalization / byte-stability / diff. |
| `test_emergence.py`, `test_mapping.py`, `test_structural.py`, `test_alert.py`, `test_delivery.py`, `test_backtest.py` | Each track-2 pipeline stage. |
| `test_recaller.py` | The re-caller core (reference assertion, uncalled codons, frame refusal, panel tagging). |
| `test_recall_orchestration.py`, `test_recall_preflight.py`, `test_recall_fill_logging.py` | The reads→consensus orchestration *logic*, offline — so a `fill` bug can't corrupt a snapshot after the expensive downloads. |
| `test_runlog.py` | The provenance logger (best-effort, never breaks the run). |

CI enforces a coverage gate (see `.github/`); the openafr core sits at ~97% coverage.

## `data/`

| Path | Contents |
|---|---|
| `data/structures/5TZ1.pdb` | The real *C. albicans* CYP51 crystal structure (+ heme + bound azole) from RCSB — the docking receptor's source. |
| `data/ligands/actives.smi` | 8 known azoles used to form the hypothesis (training). |
| `data/ligands/actives_holdout_final.smi` | 7 held-out azoles — the blind test set, never used to form hypotheses. |
| `data/ligands/decoys*.smi` | Property-matched decoys (~50 per active). |
| `data/earlywarning/erg11_reference/` | The pinned B8441 ERG11 reference CDS (accession XM_085597798.1) + `PROVENANCE.md`. |
| `data/earlywarning/snapshots/` | The versioned NCBI snapshot store (big TSVs gitignored; `INDEX.tsv` committed). Has its own `README.md`. |
| `data/earlywarning/digest/` | The delivery digest + `STATE.json`. Has its own `README.md`. |
| `data/earlywarning/runlog/` | The append-only run-provenance `.jsonl` files. |
| `data/earlywarning/recaller_sanity/` | The fixed known-genotype truth set for the sanity run. |

## `work/` — outputs, results, and pre-registrations

`work/` holds run outputs (much of it gitignored) plus two families of committed documents that
are central to the project's credibility:

- **`PREREGISTRATION_*.md`** — the rules and pass bars, written and SHA-256-frozen *before* data
  was generated: `run2`, `reliability`, `axial`, `consensus`, `auris_Y132F`, `backtest`.
- **`RESULTS_*.md`** — the write-up of every run, *including the ones that failed or that
  qualified an earlier claim*. Notable ones to read: `RESULTS_run2_holdout.md` (the headline
  result), `RESULTS_reliability.md` (the self-inflicted EF@1% caveat), `RESULTS_auris_Y132F.md`
  (the resistant-pocket port, the only `measured` structural verdict), `RESULTS_backtest.md` (the
  NOT-YET-VALIDATED verdict), `RESULTS_earlywarning_spike.md` (the "NCBI has no resistance calls"
  discovery). `RUNBOOK_recaller_run.md` is the single-host run order for the re-caller's real run.

## Root files

| File | What it is |
|---|---|
| `README.md` | What the pipeline does + how to run it. |
| `VISION.md` | Strategy, novelty, where it goes next. |
| `TODOS.md` | Deferred work — the ERG11 re-caller real run is the top item. |
| `protocol.yaml` | The frozen docking box + pass bar. Hash-pinned in `openafr/protocol.py`. |
| `environment.yml` | The version-pinned conda toolchain (rdkit, openbabel, vina, pytest, pytest-cov). |
| `pytest.ini` | Test configuration. |
| `LICENSE` | PolyForm Noncommercial 1.0.0. |

---

Recurring patterns that show up across many of these files are collected in the
[Project Overview](02-project-overview.md#the-recurring-design-patterns-learn-these-once).
Stuck on a term? → **[Glossary](06-glossary.md)**. Ready to build? → **[Contributing](07-contributing.md)**.
