# The Drug-Discovery Pipeline (Track 1)

*Audience: you understand the [biology](01-antifungal-resistance-primer.md) and have skimmed
the [overview](02-project-overview.md). This page breaks down track 1 down to the level of
individual files and functions.*

---

## The core idea, restated precisely

Azoles kill fungi by reaching a **nitrogen** to the **heme iron** at CYP51's core and jamming
it. So if you want to know whether a candidate molecule *acts like an azole*, the most direct
question is: **when docked into the pocket, does it bring a nitrogen close to the iron?**

Track 1's whole thesis is that this one geometric measurement — the **N–Fe distance** — ranks
real azoles above look-alike decoys *better than the docking program's own affinity score
does*, on the exact same docked poses. It then does the hard part: it **proves** that claim on
molecules it has never seen, under rules fixed in advance, before it will rank anything new.

---

## The end-to-end flow

There are two paths through track 1. They share the same docking engine and the same criterion.

```
  VALIDATION PATH (must pass first — proves the method)
  ────────────────────────────────────────────────────
  known azoles + matched decoys ─▶ prep ─▶ screen (dock) ─▶ validate_gate2.py
                                                                   │
                                                    exit 0 = PASS  │  exit 1 = FAIL (blocks everything)
                                                                   ▼
  DISCOVERY PATH (only justified after the gate passes)
  ─────────────────────────────────────────────────────
  novel library ─▶ filter ─▶ prep ─▶ screen (dock) ─▶ rank ─▶ report ─▶ hand top of list to a wet lab
                                                          │
                                                          ▼
                     CONFIDENCE LAYER (annotates a rank, never re-ranks)
                     ── precedent · ADMET · synth · human-CYP · stereo · protomer ·
                        ensemble · pose-convergence · resistance-retention ·
                        coordinator-identity ─▶ scorecard ─▶ A/B/C dossier
```

The rule from the [overview](02-project-overview.md#the-founding-rule): you do not get to run
the discovery path until the validation path returns a pass. A rank is a hypothesis, never a
hit.

> **A note on where the product landed.** The original headline is enrichment on
> presumed-inactive decoys (AUC 0.794). Two later pre-registered probes sharpened the claim:
> the criterion was re-run against **measured** inactives (closing the paper's largest caveat)
> and the within-azole ranking line was tested to exhaustion and **closed**. The net, honest
> product today is **novel-chemotype triage at AUC ≈ 0.81** — telling a genuinely new
> iron-coordinating chemotype apart from look-alikes — *not* ranking one working azole above a
> failed azole analogue. Parts C–E below tell that story. There is also now a **front door**
> (score your own molecule) and a **candidate-confidence layer** on top of the raw rank.

---

## Part A: The concepts the code rests on

Before the file walkthrough, three ideas you must hold.

### 1. Actives, decoys, and why decoys are the whole ballgame

- **Actives** = molecules known to work (real azole drugs). The pipeline must rank these high.
- **Decoys** = molecules *presumed inactive*, included to see whether the method can tell real
  drugs apart from plausible-looking impostors.

A ranking method looks impressive only if the decoys are genuinely hard to distinguish from
the actives. This project is paranoid about two ways decoys can secretly make the test too
easy — both documented and defended in the code:

- **The "circular decoy" trap.** The N–Fe criterion needs a nitrogen near the iron. If the
  decoys had *no aromatic nitrogen at all*, the method would separate actives from decoys
  perfectly — but really it would only be detecting *"does this molecule contain a nitrogen?"*,
  a triviality. **Defense:** every decoy is required to be a plausible iron-coordinator (it has
  an aromatic nitrogen). See `scripts/make_decoys.py`.
- **The "size shortcut" trap.** On this system the docking score correlates with molecular
  size at r ≈ −0.91 — bigger molecules get better scores almost automatically. If decoys were
  smaller than actives, the pipeline could "succeed" just by preferring heavy molecules.
  **Defense:** decoys are property-matched to actives on molecular weight, cLogP (a measure of
  greasiness), heavy-atom count, and rotatable-bond count.

If you remember nothing else about decoys: **a validation is only as honest as its decoys**,
and this project spends real effort making them hard.

### 2. Held-out actives (the blind test)

Some known azoles are used during development to *form* the hypothesis (the "training"
actives). A **separate** set is locked away and never looked at until the final graded run —
the **held-out** actives (`data/ligands/actives_holdout_final.smi`). Testing on molecules the
method has never seen is what makes the result a *blind* test rather than a memorized one.

### 3. Pre-registration and the frozen protocol

The pass/fail rules are written down and **hashed** before any data is generated:

- `work/PREREGISTRATION_run2.md` — the prose rules (what will be measured, what counts as a
  pass, how to handle edge cases). Its SHA-256 is pinned in `validate_gate2.py`.
- `protocol.yaml` — the machine-readable numbers: the docking box center/size, the search
  effort, and the pass bar. Its SHA-256 is pinned in `openafr/protocol.py`.

At run time, the gate re-checks both hashes and loudly declares any pass "not credible" if
either file changed after freezing. This is the anti-tuning guarantee: you cannot quietly move
the goalposts to make a run pass.

Here is the actual frozen protocol, annotated:

```yaml
docking:
  center_x: 70.61      # the docking box is centered on the crystal ligand's position…
  center_y: 66.28
  center_z: 4.18
  size: 26             # …and is 26 Å on a side (big enough for the long azole drugs)
  exhaustiveness: 32   # how hard Vina searches (doubled vs. an earlier run so search isn't the bottleneck)
  num_modes: 20        # keep the 20 best poses per molecule
  seed: 42             # fixed random seed → reproducible runs
gate:
  min_auc: 0.70        # to PASS: AUC must be ≥ 0.70 …
  min_ef1: 5.0         # … AND EF@1% must be ≥ 5.0×
  ef_fractions: [0.01, 0.05, 0.10]
  bedroc_alpha: 20.0
```

Both the docking script *and* the grading gate read this same file, so **what is docked can
never drift from what is graded.**

---

## Part B: File-by-file walkthrough

The `openafr/` package holds the shared, tested core; `scripts/` holds the runnable steps.
Here is the path in execution order.

### `openafr/pdbqt.py` — the shared parsers (everything depends on this)

Small and load-bearing. Three functions every gate/analysis script needs, in one canonical
place (they used to be copy-pasted across scripts — a bug-prone situation this file fixed):

- `read_poses(path)` — yields each docked pose from a Vina output file as a list of
  `(atom_name, x, y, z)`. **Skips hydrogens** — the geometry is measured on heavy atoms only.
- `read_scores(path)` — returns Vina's own affinity score for every pose (the number the
  project deliberately distrusts, kept for comparison).
- `iron_position(receptor_pdb)` — finds the `(x, y, z)` of the heme iron atom (Fe) in the
  receptor. If there's no iron, it raises an error rather than guessing — the wrong structure
  should fail loudly.
- `min_nitrogen_iron_distance(atoms, fe)` — **the heart of the criterion.** Given one pose's
  atoms and the iron position, returns the closest approach of *any* nitrogen to the iron — or
  `None` if the pose has no nitrogen at all (so it can't be a real azole pose). This single
  number is what the whole ranking is built on.

### `openafr/scoring.py` — the ranking math

Two functions, locked by `tests/test_gate_metrics.py`. All three metrics come from
`rdkit.ML.Scoring`, the canonical implementations — the project doesn't hand-roll its own
statistics.

- `metrics(ordered_flags)` — takes a list of 1s (active) and 0s (decoy) **already sorted
  best-first** and returns `(AUC, [EF@1%, EF@5%, EF@10%], BEDROC)`.
  - **AUC** — overall ranking power. 0.5 = coin flip, 1.0 = perfect.
  - **EF@x%** — enrichment factor: how many times more actives appear in the top x% of the
    ranked list than you'd expect by random chance. EF@1% is the razor-top-1% number that
    later proved fragile.
  - **BEDROC** — an AUC variant that weights the *top* of the list heavily (this is what
    matters when you can only afford to test the top few candidates).
- `report(label, ranked)` — takes `(name, key, is_active)` rows where a lower `key` = better
  (for the N–Fe criterion, `key` is the N–Fe distance). It sorts them, and — critically —
  appends any **unrankable** molecules (`key is None`, i.e. no nitrogen pose) **last, never
  dropping them**, then computes and prints the metrics. This "unrankable go last" rule is the
  pre-registered anti-cheat: dropping molecules with no usable pose would silently inflate
  every metric.

> **The `EF_FRACTIONS`/`BEDROC_ALPHA` constants here are metric *definitions*; the tunable
> pass bar (min AUC, min EF) lives in the hash-frozen `protocol.yaml`, not in the code.** The
> definition of the ruler is separate from where you draw the pass line — and only the pass
> line is frozen.

### `openafr/protocol.py` — load + verify the frozen protocol

Loads `protocol.yaml` and enforces its hash. Notable choices:

- **No PyYAML dependency.** The file is a simple, self-authored key/value map, so this module
  ships a tiny strict parser (`_parse`) rather than dragging in a YAML library — keeping the
  whole pipeline runnable with just rdkit + the docking toolkit.
- `verify()` — recomputes the file's SHA-256 and compares it to the pinned `PROTOCOL_SHA`. A
  *soft* check: it never blocks a run, but it makes a tampered protocol impossible to miss in
  the output. (Changing the protocol legitimately is a deliberate act: edit the file, re-freeze
  with `python -m openafr.protocol --freeze`, update `PROTOCOL_SHA`, and record why.)
- `--shell` mode emits the docking parameters as shell variables so `screen.sh` reads the
  exact same numbers the gate does.

### `scripts/prep_receptor.py` — build the docking target

Turns the raw crystal structure `data/structures/5TZ1.pdb` into the prepared receptor
`work/receptor.pdbqt` (chain A + the heme). It asserts the extracted atoms match a tracked
reference (`work/receptor_A.pdb`) so the receptor can't silently change. `work/receptor.pdbqt`
is gitignored, so a fresh clone rebuilds it from the tracked crystal file. (Verified by
re-docking a known azole, VT-1161, and reproducing its iron-bound pose.)

### `scripts/prep_ligands.py` — SMILES → 3D molecules

Ligands start as **SMILES** strings (a compact text encoding of a molecule's structure, e.g.
in `data/ligands/*.smi`). This converts each to a 3D structure ready for docking,
**deterministically**, and — following the honesty pattern — **logs failures instead of
silently dropping them.** A molecule that vanishes from a run is a molecule that can't quietly
flatter the results.

### `scripts/filter_library.py` — enforce the applicability domain (discovery path only)

Before spending any docking compute on a *novel* library, this drops molecules the gate's
evidence does not cover, and writes every exclusion to a `.tsv` with a reason (nothing is
silently dropped). Two filters, both pre-registered:

- **> 45 heavy atoms** → excluded. This is the declared applicability domain; big
  posaconazole/itraconazole-class ligands are an explicit blind spot the validation didn't
  cover.
- **No aromatic nitrogen** → excluded. The N–Fe criterion is simply undefined for a molecule
  with no nitrogen to measure.

Deliberately, it adds **no** un-registered drug-likeness filters — that would be the same
tuning trap the gate is built to avoid.

### `scripts/screen.sh` — the docking engine (shared by both paths)

The single place any ligand library gets docked. Key properties:

- **One frozen protocol for every molecule.** The box and search parameters come from
  `protocol.yaml` (via `openafr/protocol.py`), never hand-copied — so validation and screening
  can't drift apart.
- **Parallel but deterministic.** Docking runs are independent, so it dispatches one Vina
  process per CPU core. The fixed seed is applied identically and independently to each ligand
  (no shared RNG), so a parallel run is bit-for-bit identical to a serial one, just faster.
- **Resumable.** A ligand whose output already exists is skipped, so an interrupted screen
  picks up where it stopped. Failures are counted and reported, never silently dropped.
- Docks into `work/receptor.pdbqt` by default; set `RECEPTOR=` to dock into a mutant (this is
  how the *C. auris* resistance ports are screened — see track 2's structural stage).

`scripts/run_screen2.sh` is a thin wrapper that runs the specific **validation** screen (the 7
held-out azoles + their decoys) through this same engine.

### `scripts/validate_gate2.py` — THE gate

This is the program the founding rule is about. It implements `PREREGISTRATION_run2.md`
exactly:

1. **Re-checks the hashes** — the pre-registration's and (via `openafr.protocol.verify`) the
   protocol's. Any mismatch prints a loud "not credible" warning.
2. For each docked molecule, computes the **N–Fe distance** (mode C, the primary criterion)
   using `min_nitrogen_iron_distance`. It also computes Vina's own score (mode A) and a
   combined rank (mode D) — but those are **reported only**, they don't gate.
3. Follows the pre-registered handling rule: an active that produced **no usable pose is
   ranked last, never dropped.**
4. Calls `scoring.report`, then checks the result against the frozen pass bar:
   **AUC ≥ 0.70 AND EF@1% ≥ 5.0× on mode C.**
5. **Exits 0 (pass) or nonzero (fail).** A nonzero exit blocks the pipeline — this is the
   "block the screen" mechanism in action.

There are sibling gates for related pre-registered experiments — `validate_gate_reliability.py`
(the multi-seed reliability check), `validate_gate_auris.py` (the Y132F resistant-mutant port),
`validate_gate_axial.py`, `validate_gate_consensus.py` — each implementing its own
pre-registration. See the [Codebase Reference](05-codebase-reference.md).

### `scripts/rank_candidates.py` — turn poses into a ranked shortlist (discovery path)

Only justified *after* the gate has passed. Ranks a screened novel library by the one validated
criterion — mode C, the N–Fe distance (smaller = better). Vina's affinity score rides along as
a secondary column only (on this system it ranked worse than random, so it doesn't decide
order). Same anti-cheat rule: a molecule with no usable pose is ranked **last, never dropped.**
It explicitly does **not** re-run the gate — it trusts that the gate already passed and reminds
you, in its own docstring, that a rank is a hypothesis, never a hit.

### `scripts/report_candidates.py` — the human-readable hand-off (discovery path)

Turns the ranked list into a Markdown report a wet-lab collaborator can actually read:
per-candidate iron-coordination evidence, the frozen-protocol provenance (including a live
check that `protocol.yaml` is unmodified), and the caveats that keep the list honest.

### `scripts/analyze_poses.py`, `pose_analysis.py`, `make_pose_view.py`, `render_views.sh`

Supporting analysis and visualization of the metal-coordination geometry — the tools used to
inspect *why* a pose does or doesn't coordinate the iron.

### The front door — `scripts/applicability_triage.py` + `scripts/score_molecule.py`

The validated method now has a **front door**: give it a SMILES string and it runs the whole
pipeline on your own molecule, honestly.

- `applicability_triage.py` answers the question the ranker can't: *is your molecule even
  inside the region the validation covers?* It returns one verdict — `in-envelope-novel`,
  `near-known`, `out-of-domain`, `out-of-scope`, or `unparseable` — and only the first two are
  inside the validated envelope. This runs anywhere `rdkit` is installed, no docking.
- `score_molecule.py` chains triage → prep → dock → N–Fe read in one command (~9 s/molecule on
  a laptop). Two honest gates are load-bearing: **(1)** only in-envelope molecules are docked;
  **(2)** if the docking tools or receptor are missing, the geometry read is reported `PENDING`
  with the exact command to finish it — a number is never invented. `--json` makes it scriptable
  as a public benchmark. Sanity-check your install with voriconazole (should report N–Fe ≈ 2.69 Å).

---

## Part C: The result, and the caveat the project inflicted on itself

On the pre-registered blind test — 7 held-out azoles hidden among 348 property-matched decoys,
rules SHA-256-frozen in advance:

| Ranking method | AUC | EF@1% | Verdict |
|---|---|---|---|
| **N–Fe distance** (primary) | **0.794** | 12.68× | **PASS** |
| Vina affinity score | 0.471 | 0.00× | worse than random |

- **Permutation test** (20,000 shuffles): p = 0.0028 — the AUC advantage is very unlikely to be
  chance.
- **Bootstrap 95% CI** on the AUC: 0.590–0.946.

**The self-inflicted caveat.** A *later*, also-pre-registered reliability test
(`work/RESULTS_reliability.md`) found that the **EF@1% = 12.68× figure was inflated by
single-shot docking.** When each molecule is docked over 5 random seeds instead of once, that
top-1% number collapses to 0.00× under *every* aggregation method tried — because
property-matched decoys reliably produce genuine iron-coordinating poses right at the very top.
What **survives** is the broad enrichment: AUC stays ~0.79–0.85, EF@5% stays ~5.6–8.5×.

So the honest headline is: **a trustworthy top-~5% shortlist, not razor-top ranking.** The AUC
is the durable number; EF@1% is not. The original result table is kept in the README *as the
graded record of the run as pre-registered*, with this caveat attached — the project does not
delete an inflated number it once reported, it annotates it. **That is the discipline the whole
pipeline exists to enforce, applied to itself.**

---

## Part D: Honest limitations (worth internalizing)

- **n = 7** held-out actives in the original run — a small test set; that AUC's 95% CI lower
  bound dips below the pass bar. (Later probes below use larger, independent sets.)
- **The presumed-inactive caveat is closed, and it located a ceiling** — see Part E.
- **One rigid receptor** — *C. albicans* 5TZ1, not *C. auris*; no induced-fit flexibility.
  The ensemble look (Part E) tested whether extra crystallographic conformers lift this; they
  do not lift the within-azole ceiling.
- **Applicability domain** — ligands ≤ 45 heavy atoms; the posaconazole/itraconazole class is a
  declared blind spot.
- **The defensible product is novel-chemotype triage, not within-class azole ranking** (Part E).
- **This validates a ranking criterion, not a drug.** A high rank is a hypothesis for a wet
  lab. **No wet-lab collaborator yet is the one critical open dependency.**

---

## Part E: The criterion, stress-tested to its ceiling

Everything after the original run-2 result is a series of pre-registered attempts to *break* the
criterion or find where it stops working. Together they are why the honest product is
"novel-chemotype triage," stated precisely.

### Measured inactives close the largest caveat — and locate the ceiling (look #6)

The 348 run-2 decoys were *presumed* inactive. `openafr/inactives.py` +
`scripts/make_verified_inactives.py` build a **measured** decoy set: 279 compounds with
*consistent measured non-inhibition* of *C. albicans* in ChEMBL (single-variable swap — every
other construction rule held fixed, only "inactive" changes from assumed to measured). Re-run
unchanged (`scripts/validate_gate_verified.py`), the criterion scores **AUC 0.716 — the gate
passes**, so the original result is not a decoy-construction artifact.

But the same run *located the ceiling*, and this is the honest headline:

- vs. **non-azole chemotypes**: AUC **0.810**
- vs. **azole analogues measured not to work**: AUC **0.650** — below its own bar.

So the criterion mostly detects **chemotype**, not **potency**. (`work/RESULTS_verified_inactives.md`.)

### The within-azole ranking line — tested to exhaustion, then closed (looks #7–#9)

Can geometry rank a *working* azole above a *failed* azole analogue? Three pre-registered looks
say no, on both the rigid and the ensemble receptor:

- **Look #7 (mode F, channel engagement, `validate_gate_channel.py`)** — the one permitted
  orthogonal feature. AUC 0.590 < 0.70. FAIL. The failure isn't re-measuring size (size-residualized 0.612).
- **Look #8 (ensemble, `validate_gate_ensemble.py`)** — the min N–Fe over an exhaustive
  3-structure crystallographic ensemble (`openafr/ensemble.py`). PASSED at AUC 0.750 on n=7,
  briefly re-opening the line — but its own PASS branch demanded an independent confirmation.
- **Look #9 (independent confirmation, `validate_gate_ensemble_confirm.py`)** — that ensemble
  PASS did **not** replicate on an independent n=200 measured-active azole set: AUC 0.682 < 0.70
  (real signal, p<0.0001, but sub-threshold).

**Settled conclusion:** pose geometry over a fixed *C. albicans* CYP51 target — rigid **or** full
crystallographic ensemble — cannot rank working azoles above failed analogues to the 0.70 bar.
That is a property of the *question*, not of any one feature or conformer. The line is closed;
re-opening it needs genuinely flexible/induced-fit poses (a new dataset and pre-registration).

### Generalization to data the criterion has never touched (holdouts)

Two pre-registered holdouts (#82) test whether the *frozen* criterion generalizes:

- **External holdout** (`make_external_holdout.py` / `validate_gate_external.py`) — a fresh
  active/decoy set no gate has ever scored.
- **Temporal holdout** (`make_temporal_holdout.py` / `validate_gate_temporal.py`) — a
  *publication-date* split that holds out the **newest-literature azoles (2024–25)**. The
  criterion, defined on older chemistry, still generalizes: **AUC 0.754**.

### Metal-aware baselines — a competent opponent (looks #63, #83)

Vina mishandles metal coordination, so the honest next question is whether the criterion also
beats a *metal-aware* scorer. `openafr/gnina.py` + `openafr/baselines.py` turn GNINA's CNN
rescoring (and other rivals) into ranking keys on the **identical** held-out poses, graded by
`rescore_gnina.py` / `rescore_baselines.py`. The scorer binaries are Linux/x86-native, so they
run on a cloud host (`run_*_rescore.sh`) while parsing + grading stay reproducible anywhere.

### The heme model, audited (look #72)

`scripts/audit_heme_model.py` audits the heme/iron model the whole criterion rests on — the one
structural assumption a skeptic would attack first.

---

## Part F: The candidate-confidence layer (rank → dossier)

A geometry rank is necessary but not sufficient. Before a molecule is proposed for a bench, the
pipeline now profiles it along a stack of **orthogonal per-candidate axes**, fuses them into a
scorecard, and emits a machine-readable **A/B/C dossier**. Every axis *annotates* a rank — none
re-ranks, none deletes a row, and none is a hard gate.

| Axis | Module / script | What it flags |
|---|---|---|
| Assay precedent | `precedent.py` / `check_precedent.py` | already measured against CYP51 (and failed)? ChEMBL + PubChem + BindingDB. `no-precedent` ≠ untested. |
| File-drawer bound | `filedrawer.py` | quantifies how much a `no-precedent` could be hiding (a leak *band*, not false precision). |
| ADMET / alerts | `admet.py` / `admet_filter.py` | Ro5/Veber, PAINS/BRENK/NIH, solubility, a coarse hERG hint. |
| Synthesizability | `synth.py` / `synth_filter.py` | SAscore + optional catalogue presence. |
| Human-CYP off-target | `cyp_offtarget.py` / `cyp_offtarget_screen.py` | the type-II heme-ligation warhead the target *rewards* is the liability the human panel *dreads*. |
| Stereochemistry | `candidate_stereo.py` | every candidate isomer-explicit. |
| Protonation/tautomer | `candidate_protomer.py` | the pH-7.4 microstates. |
| Ensemble consistency | `candidate_ensemble.py` | does the pose survive across CYP51 conformers? |
| Pose convergence | `pose_convergence.py` | did one search settle on one geometry or scatter? |
| Resistance retention | `candidate_resistance.py`, `candidate_resistance_addatom.py` | does iron-reach survive in the C. auris resistant pocket? |
| Seed-stability + human-selectivity | `shortlist_confidence.py` | robust to the docking seed; selective over human CYP51. |
| **Coordinator identity** | `coordinator.py` / `coordinator_identity.py` | **is the rank azole-mechanism, or an artifact?** (see below) |

`scorecard.py` is the pure decision logic that consolidates these (#88); `build_dossier.py`
(over `openafr/dossier.py`) fuses geometry + precedent + ADMET + CYP + synth + coordinator
identity into one dossier per molecule with an **A/B/C** priority.

### The coordinator-identity finding — why the sole Priority-A fell

The validated criterion (`min_nitrogen_iron_distance`) measures the closest approach of *any*
nitrogen — it doesn't ask *which*. But the method was validated on the azole mechanism: an
**aromatic ring** nitrogen donating a lone pair to the iron. On a novel library, a molecule can
post an identical tight distance through a **nitrile**, amine, or amide nitrogen — an
out-of-mechanism artifact. `openafr/coordinator.py` classifies the coordinating N purely from
pose geometry (bonded-neighbour topology; the one threshold, 1.25 Å, sits in the empty gap
between a C≡N triple bond and every single/aromatic C–N).

It validates on all 6 azole controls (in-mechanism), then reshapes the shortlist: the sole
Priority-A candidate (verinurad) turned out to coordinate through its **nitrile** in every seed
(its aromatic pyridine N reaches only the band edge), so it was **demoted A→B**; two other
"strong, selective" hits were nitrile artifacts too. **After the mechanism is enforced, no
molecule clears automatic-A** (dossier A=0 / B=14 / C=2); the strongest survivor is vatalanib
(in-mechanism, seed-stable, human-selective).

> **The honest twist, kept on purpose.** The obvious fix — restrict the criterion to an aromatic
> N — was pre-registered and tested (`validate_gate_aromatic.py`). It **lowers** the gate AUC
> (0.729 < 0.794), because the real active *luliconazole* coordinates via its own nitrile in
> pose. Lesson: a non-aromatic coordinator is *unvalidated, not disproven*. So the frozen
> criterion is left **unchanged**, and coordinator identity is used only as a per-candidate flag.

---

## Part G: The drift harness (`scripts/repro.py`)

One command re-checks the entire frozen surface (#89): `protocol.yaml`, **every**
`work/PREREG_*.sha256` (all the pre-registrations plus the frozen datasets they pin), and that
the consolidated candidate scorecard still re-derives byte-for-byte from its committed per-axis
inputs — exiting nonzero on any drift. It runs on every PR (`.github/workflows/repro.yml`).
Re-docking stays out of scope (pose trees are gitignored and take hours); the harness guarantees
nothing that does *not* require docking can drift unnoticed.

---

## How to run it yourself

See the [main README](../README.md#reproducing) for the exact commands and the
[Contributing](07-contributing.md) page for environment setup. The short version:

```bash
# validation path (must pass before screening anything novel)
python scripts/prep_receptor.py
python scripts/prep_ligands.py data/ligands/actives_holdout_final.smi work/run2_ligands
python scripts/prep_ligands.py data/ligands/decoys_holdout.smi        work/run2_ligands
./scripts/run_screen2.sh
python scripts/validate_gate2.py work/screen2 work/receptor_A.pdb   # exit 0 = pass
pytest
```

Next: **[The Early-Warning Pipeline →](04-early-warning-pipeline.md)** (track 2), or the
**[Codebase Reference →](05-codebase-reference.md)** for the full file map.
