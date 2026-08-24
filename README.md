# OpenAFR — a self-validating antifungal screening pipeline

An open-source computational pipeline for finding candidate antifungal molecules against
drug-resistant *Candida auris*, built around one rule:

> **The pipeline must prove it can re-discover the drugs we already know work, on molecules
> it has never seen, before it is allowed to rank anything new.**

If it cannot, the validation gate exits nonzero and blocks the screen.

The repo now holds **two tracks** against the same pathogen and enzyme:

1. **Drug discovery** (the original, validated track) — the geometry-over-docking-score
   triage tool described immediately below.
2. **Genomic early-warning surveillance** (a newer, honestly-incomplete track) — watching
   NCBI for emerging *C. auris* azole-resistance mutations and giving each one a structural
   verdict. See [Second track: genomic early-warning surveillance](#second-track-genomic-early-warning-surveillance).

Everything from here to that section is about track 1.

> **New here?** The [`wiki/`](wiki/README.md) is written to be picked up cold — a
> [biology primer](wiki/01-antifungal-resistance-primer.md) for anyone new to antifungal
> resistance, a [project overview](wiki/02-project-overview.md), code-level deep dives of
> [track 1](wiki/03-drug-discovery-pipeline.md) and [track 2](wiki/04-early-warning-pipeline.md),
> a [file-by-file codebase reference](wiki/05-codebase-reference.md), a
> [glossary](wiki/06-glossary.md), and a [contributing guide](wiki/07-contributing.md).

## Score your own molecule

The validated method has a front door: give it a SMILES string and it tells you, honestly,
whether the method even applies to your molecule — and if so, how closely it coordinates the
heme iron the way real azole drugs do. **Read [`docs/INTERPRETING_A_SCORE.md`](docs/INTERPRETING_A_SCORE.md)
before acting on any number: a score here is a hypothesis for a wet lab, never a hit.**

**1. Install the environment (once).** Resolves on both Apple Silicon and Linux:

```
conda env create -f environment.yml     # or: mamba env create -f environment.yml
```

**2. Triage — does the validated method even cover your molecule?** Runs anywhere in seconds,
no docking, no cloud:

```
conda run -n openafr python scripts/applicability_triage.py \
    --smiles "CC(C)c1nc(cs1)..." --name my_compound
```

You get one verdict: `in-envelope-novel`, `near-known`, `out-of-domain`, `out-of-scope`, or
`unparseable`. Only the first two are inside the validated envelope; the rest are reported and
stopped, because the method makes no claim about them.

**3. Score — the full SMILES → geometry pipeline in one command.** Triage → prep → dock →
N-to-iron geometry read, on your laptop in ~9 s per molecule:

```
conda run -n openafr python scripts/prep_receptor.py                # once, builds work/receptor.pdbqt
conda run -n openafr python scripts/score_molecule.py \
    --smiles "C[C@@H](C1=NC=NC=C1F)[C@](CN2C=NC=N2)(C3=C(C=C(C=C3)F)F)O" --name voriconazole
```

**Confirm your install** with that voriconazole command (a real azole drug): triage returns
`near-known`, and the geometry read should report **closest N–Fe ≈ 2.69 Å**, inside the
validated actives' iron-bound band (2.47–2.88 Å). If you see that, the whole pipeline works on
your machine. To score without docking (triage only, runs anywhere), add `--triage-only`; the
tool also prints the exact commands to finish a run if `vina`/`obabel` are missing, and never
fakes a number.

You can pass a whole library instead of `--smiles`: a `SMILES<TAB>name` file as the first
argument. Add `--json` to get one machine-readable object on stdout (human text goes to
stderr), so the front door is scriptable as a benchmark. What a readout does and does not
establish is spelled out in [`docs/INTERPRETING_A_SCORE.md`](docs/INTERPRETING_A_SCORE.md).

## Why antifungal resistance

*Candida auris* is a WHO critical-priority fungal pathogen: multidrug-resistant, lethal, and
spreading in hospitals. Antifungal R&D is dramatically underfunded relative to antibacterial
R&D. The azole drugs work by jamming a fungal enzyme, **CYP51**, and they are failing.

## Why "self-validating"

Molecular docking produces *hypotheses*, not hits. Published hit rates are low and scoring
functions are noisy. Most amateur docking repositories publish ranked lists with no evidence
the ranking means anything. This project treats that evidence as the deliverable.

## Result

On a **pre-registered** held-out test — rules and pass thresholds written and SHA-256 frozen
before any data was generated (`work/PREREGISTRATION_run2.md`) — the pipeline was asked to
find 7 azole antifungals it had never seen, hidden among 348 property-matched decoys:

| ranking method | AUC | EF@1% | |
|---|---|---|---|
| **iron-approach distance** (primary) | **0.794** | **12.68x** | **PASS** |
| AutoDock Vina affinity score | 0.471 | 0.00x | worse than random |

Permutation test (20,000 shuffles): **p = 0.0028**. Bootstrap 95% CI: 0.590–0.946.

> **Caveat added by a later pre-registered check (kept here on purpose).** The AUC and its
> permutation p-value are robust, but the **EF@1% = 12.68x figure was inflated by single-shot
> docking.** A pre-registered reliability test (`work/RESULTS_reliability.md`) found that
> top-1% number collapses to 0.00x once each ligand is docked over 5 seeds — under *every*
> aggregation tried — because property-matched decoys reliably reveal genuine iron-coordinating
> poses at the very top. Broad enrichment survives (AUC ~0.79–0.85, EF@5% ~5.6–8.5x). The
> honest headline is a **trustworthy top-~5% shortlist, not razor-top ranking.** The table
> above is the graded record of the run as pre-registered; this note qualifies how to read it.

**The finding:** for azoles against CYP51, *how closely a molecule brings a nitrogen to the
heme iron* discriminates real drugs from look-alike decoys far better than the docking
program's own affinity score — on the exact same poses. Geometry beats the scoring function.

> **Confirmed on measured inactives (2026-08-23, look #6).** Those 348 decoys were *presumed*
> inactive — never tested. The criterion has now been re-run, unchanged, against **279 compounds
> measured not to inhibit *C. albicans*** (`work/RESULTS_verified_inactives.md`): **AUC 0.716,
> gate passes**, so the result above is not an artifact of how the decoys were built. Two things
> got *worse* under measurement and are the honest headline: **7 measured inactives outrank the
> best real drug**, and a pre-specified split puts the criterion at **AUC 0.810 against non-azole
> chemotypes but 0.650 against azole analogues measured not to work** — below its own bar. Most
> of what it detects is *chemotype*, not potency.

Full write-ups, including a run that **failed** its gate first, are in `work/`.

## Honest limitations

- n = 7 held-out actives; the 95% CI lower bound sits below the pass bar.
- Decoys are *presumed* inactive, not experimentally verified.
- One rigid receptor (*C. albicans* 5TZ1, not *C. auris*); no induced fit.
- Applicability domain: ligands ≤ 45 heavy atoms. Posaconazole/itraconazole-class ligands
  are a declared blind spot.
- **This validates a ranking criterion, not a drug.** A high rank is a hypothesis for a wet
  lab, never a hit.

## Layout

```
data/structures/5TZ1.pdb          C. albicans CYP51 + heme + bound azole (RCSB)
data/ligands/actives.smi          8 known azoles (training)
data/ligands/actives_holdout_final.smi   7 held-out azoles (never used to form hypotheses)
data/ligands/decoys*.smi          property-matched decoys, 50 per active
data/ligands/verified_inactives.smi  compounds MEASURED not to inhibit (n=279)

protocol.yaml                     the frozen run protocol (box, search effort, pass bar), hash-pinned

openafr/pdbqt.py                  shared Vina/PDBQT/receptor parsers
openafr/scoring.py                metrics() / report() — the ranking math
openafr/protocol.py               loads + hash-verifies protocol.yaml

scripts/prep_receptor.py          5TZ1.pdb -> work/receptor.pdbqt, chain-A+heme, matches the anchor
scripts/prep_ligands.py           SMILES -> 3D, deterministic, failures logged not dropped
scripts/make_decoys.py            property-matched decoy generation from ChEMBL
scripts/make_verified_inactives.py   measured-inactive set from ChEMBL bioactivity
scripts/validate_gate_verified.py    grades the criterion vs measured inactives (look #6)
scripts/filter_library.py         applicability-domain filter — keep only what the gate validated
scripts/screen.sh                 the docking engine — one fixed protocol, all cores, resumable
scripts/run_screen2.sh            thin wrapper: run-2 held-out validation via screen.sh
scripts/rank_candidates.py        screen poses -> candidate list ranked by the validated criterion
scripts/report_candidates.py      candidate list -> human-readable Markdown report with the evidence
scripts/analyze_poses.py          metal-coordination analysis
scripts/validate_gate*.py         the gate — exits nonzero when it cannot recover knowns

tests/                            known-answer tests for the gate math, parsers, and protocol
work/PREREGISTRATION_run2.md      rules fixed in advance, hash-frozen
work/RESULTS_*.md                 every run, including the failure
```

## Reproducing

```bash
# toolkit (Apple Silicon / Linux), version-pinned in environment.yml
conda env create -f environment.yml     # or: mamba env create -f environment.yml
conda activate openafr

# 1. Prepare the docking receptor: 5TZ1 chain A + heme -> work/receptor.pdbqt.
#    work/receptor.pdbqt is gitignored, so a fresh clone must build it first.
python scripts/prep_receptor.py   # -> work/receptor.pdbqt (asserts atoms match work/receptor_A.pdb)

# 2. Prepare the run-2 held-out ligands — 7 azole actives + 350 property-matched
#    decoys — into the directory run_screen2.sh reads from (work/run2_ligands):
python scripts/prep_ligands.py data/ligands/actives_holdout_final.smi work/run2_ligands
python scripts/prep_ligands.py data/ligands/decoys_holdout.smi        work/run2_ligands

# 3. Dock every ligand under the one frozen protocol, then grade the result.
#    The screen docks against work/receptor.pdbqt (step 1); the gate reads the
#    heme-iron position from the tracked work/receptor_A.pdb — two files, both
#    the same 5TZ1 receptor in different formats.
./scripts/run_screen2.sh                                            # -> work/screen2/
python scripts/validate_gate2.py work/screen2 work/receptor_A.pdb   # exit 0 = pass

pytest   # known-answer tests for the gate math, parsers, and protocol
```

The gate script re-checks two hashes — the pre-registration (`work/PREREGISTRATION_run2.md`)
and the frozen protocol (`protocol.yaml`) — and flags any pass as not credible if either was
edited after being frozen. The docking script reads its box and search parameters from the
same `protocol.yaml`, so what is docked and what is graded can never silently drift apart.
Changing the protocol is a conscious act: edit `protocol.yaml`, re-freeze with
`python -m openafr.protocol --freeze`, update `PROTOCOL_SHA`, and record why.

## Score your own molecule

One command takes a SMILES string through the whole validated pipeline — applicability triage
→ ligand prep → docking → the nitrogen-to-iron geometry read — and prints an honest verdict.
No cloud VM: the docking tools (`vina`, `openbabel`) install locally from `environment.yml` on
Apple Silicon and Linux alike.

```bash
conda activate openafr
python scripts/prep_receptor.py                      # once: build work/receptor.pdbqt (gitignored)

# a single molecule (this is voriconazole):
python scripts/score_molecule.py \
  --smiles "C[C@@H](C1=NC=NC=C1F)[C@](CN2C=NC=N2)(C3=C(C=C(C=C3)F)F)O" --name voriconazole

# or a whole library file of `SMILES<TAB>name` rows:
python scripts/score_molecule.py candidates.smi

# just the triage (no docking, runs anywhere rdkit is installed):
python scripts/score_molecule.py --smiles "<SMILES>" --triage-only
```

Two honest gates are load-bearing. **(1)** Only molecules the triage puts *inside* the validated
envelope (in-envelope-novel / near-known) are docked; anything out-of-scope, out-of-domain, or
unparseable is reported and stopped, never docked — the method makes no claim there. **(2)** If
the docking tools or a prepared receptor are absent, the geometry read is reported `PENDING`
with the exact command to finish it; a number is never invented or silently skipped.

Read every geometry readout as a **hypothesis for a wet lab, never a hit**: reaching the iron
like the validated actives is necessary but *not* sufficient (property-matched decoys reach it
too — the decoy ceiling), and the durable validated metric is AUC ≈ 0.79. The pipeline was run
end-to-end on known azoles as a sanity check — see [`work/RESULTS_frontdoor_dogfood.md`](work/RESULTS_frontdoor_dogfood.md).

## Design notes worth knowing

Two traps this pipeline is explicitly built to avoid, both documented in `work/`:

1. **Circular decoys.** The gate requires a nitrogen near the heme iron. If decoys had no
   aromatic nitrogen, the gate would separate actives from decoys perfectly while really only
   asking "does this molecule contain a nitrogen?" Every decoy is therefore required to be a
   plausible iron-coordinator.
2. **Size shortcuts.** Docking score correlates with heavy-atom count at r = −0.91 on this
   system. Unmatched decoys would let the pipeline "succeed" by preferring heavy molecules.
   Decoys are matched on MW, cLogP, heavy atoms, and rotatable bonds.
3. **Presumed inactivity.** A decoy nobody tested is an assumption, not a measurement, and it was
   the largest caveat on the published result. `data/ligands/verified_inactives.smi` (n=279,
   `scripts/make_verified_inactives.py`) is built only from compounds with **consistent measured
   inactivity** in ChEMBL — at least one inactive-grade record and zero active or ambiguous ones,
   including on the CYP51 enzyme target.

## Status

Prove-the-core spike: complete and passing. The `screen/` stage now exists as a single,
parallel docking engine (`scripts/screen.sh`) that docks any prepared ligand library under the
one frozen protocol; `scripts/rank_candidates.py` turns its poses into a candidate list ranked
by the validated iron-approach criterion. Run-2 validation now goes through the same engine, so
validation and screening can never drift apart. The `filter/` stage
(`scripts/filter_library.py`) now enforces the pre-registered applicability domain on
a novel library before any docking compute is spent, so the ranked list only ever
contains molecules the gate's evidence covers. The `report/` stage
(`scripts/report_candidates.py`) then turns that ranked list into a Markdown report a
wet-lab collaborator can read: per-candidate iron-coordination evidence, the frozen-
protocol provenance (including a live check that `protocol.yaml` is unmodified), and
the caveats that keep the list honest. With filter and report in place, the novel-
screening path — filter → prep → screen → rank → report — is complete end to end.

No wet-lab collaborator yet — that is the critical open dependency, since nothing here is
validated until someone puts candidates on real fungus. And screening a novel library is only
justified *after* `validate_gate2.py` passes: a high rank is a hypothesis for a wet lab, never
a hit.

To screen a novel library once the gate has passed:

```bash
python scripts/filter_library.py my_library.smi -o work/in_domain.smi -x work/excluded.tsv
                                                                    # keep only what the gate validated
python scripts/prep_ligands.py work/in_domain.smi work/screen_ligands   # SMILES -> 3D
./scripts/screen.sh                                                  # dock all cores -> work/screen/
python scripts/rank_candidates.py work/screen -o work/candidates.tsv # ranked by validated criterion
python scripts/report_candidates.py work/screen -o work/report.md    # human-readable report + evidence
```

The `filter/` step drops molecules the gate's evidence does not cover *before* any
docking compute — > 45 heavy atoms (the pre-registered applicability domain, the
posaconazole/itraconazole blind spot) and anything with no aromatic nitrogen (the
N–Fe criterion is undefined for it). Nothing is silently dropped: every exclusion
is written to `work/excluded.tsv` with a reason. It deliberately adds *no*
un-registered drug-likeness filters — that would be the same tuning trap the gate
is built to avoid.

Receptor preparation is now scripted end to end (`scripts/prep_receptor.py`), so the
reproduce recipe runs from a fresh clone: `data/structures/5TZ1.pdb` -> `work/receptor.pdbqt`,
with the extracted chain-A + heme atoms asserted identical to the tracked `work/receptor_A.pdb`
the gate grades against. Verified by re-docking VT-1161 (top −12.1 kcal/mol; iron-bound pose
at N–Fe 2.63 Å, matching `work/RESULTS_redock_VT1.md`).

## Second track: genomic early-warning surveillance

A separate pipeline that watches **NCBI Pathogen Detection** for *C. auris* isolates and
tries to flag *emerging* azole-resistance mutations early — before they show up in the
published surveillance record — then attaches a **structural so-what** (does fluconazole
still fit the pocket if this spreads?) to each flag. Track 1 finds new drugs; track 2
watches the enemy evolve against the drugs we have.

### Honest status: scaffold complete, blocked on one missing piece

The whole chain is **built, tested, and pulls real NCBI data** — but it cannot produce a
validated warning yet, and the docs say so plainly. The load-bearing finding, from the
initial spike ([work/RESULTS_earlywarning_spike.md](work/RESULTS_earlywarning_spike.md)):

> **NCBI runs no AMR pipeline on *C. auris*.** It is a metadata + genome-pointer feed, not
> a resistance feed. `AMR_genotypes` is empty for all ~29k isolates — there is no
> resistance call to read off the feed.

So the resistance signal has to be *manufactured by us*, via an **ERG11 re-caller** (SRA
reads → azole-resistance mutation call). **That re-caller has not been built** — it is the
single blocker for the whole track, and it is the top open item in
[TODOS.md](TODOS.md). Because `erg11_call` is empty on every real isolate, the
pre-registered backtest ([work/RESULTS_backtest.md](work/RESULTS_backtest.md)) came back
**NOT-YET-VALIDATED**: a walk-forward over 8 years of real snapshots flagged the target
mutation 0 times (it correctly refuses to invent a warning from absent calls), while the
same detector fed synthetic calls fires **397 days before** the reference date — proving
the machinery works and the gap is the input, not the method. What *is* real: the arrival
budget clears the bar (median 97-day lag from sample collection to NCBI visibility), so the
re-caller is worth building.

### The pieces (all shipped and tested; each has a `work/RESULTS_*.md` write-up)

| Stage | Library | CLI | What it does |
|---|---|---|---|
| snapshot store | `openafr/earlywarning.py` | `scripts/snapshot_ncbi_auris.py` | pull → normalized, hash-versioned per-release snapshot + baseline/diff |
| emergence detection | `openafr/emergence.py` | `scripts/detect_emergence.py` | flags first-appearing / rising ERG11 substitutions (the signal NCBI can't give) |
| mutation → pocket mapping | `openafr/mapping.py` | `scripts/map_mutation.py` | a flagged token → where the residue sits in the modeled azole pocket; can build the mutant |
| structural so-what | `openafr/structural.py` | `scripts/structural_sowhat.py` | does fluconazole still fit if this spreads, and how much worse — or an honest decline |
| alert composition | `openafr/alert.py` | `scripts/compose_alert.py` | the one-screen brief a clinician acts on |
| delivery + schedule | `openafr/delivery.py` | `scripts/deliver_alert.py` | delivers only on genuine new news; `.github/workflows/earlywarning.yml` runs it |
| backtest / validation | `openafr/backtest.py` | `scripts/backtest_earlywarning.py` | the pre-registered credibility test for the whole track |

The initial feasibility probe is `scripts/probe_ncbi_auris.py`. Persistence and delivery
state have their own READMEs: [data/earlywarning/README.md](data/earlywarning/README.md)
(snapshot store) and [data/earlywarning/digest/README.md](data/earlywarning/digest/README.md)
(delivery log). Every stage has offline tests under `tests/` (`test_earlywarning_snapshot`,
`test_emergence`, `test_mapping`, `test_structural`, `test_alert`, `test_delivery`,
`test_backtest`).

### Try it without NCBI

Each stage has a synthetic `demo` mode so a fresh clone can see the chain fire end to end
without pulling reads or waiting on the missing re-caller:

```bash
python scripts/detect_emergence.py demo
python scripts/structural_sowhat.py estimate Y132F K143R F126L V125A TR34
python scripts/deliver_alert.py demo --markdown   # deliver -> re-run no-op -> escalation
python scripts/backtest_earlywarning.py mechanism # detector produces lead time once fed calls
```

### Coverage caveat (say it out loud)

The feed is **US-centric** (~88% USA). "Global early warning" is honestly "US early warning
+ thin international." The alert copy and scope claims are written to reflect that.

## License

**[PolyForm Noncommercial License 1.0.0](LICENSE)** — © The Antifungal Resistance Project.

This is a *source-available* license, not an OSI open-source license:

- **Free** for noncommercial use — academic and public research, nonprofits,
  educational, health, and government institutions, and personal/hobby use.
- **Commercial use requires a separate license** from the foundation. If you want
  to use OpenAFR for commercial purposes, contact us — see the site.

Copyright is held by the foundation so it can offer commercial licenses; outside
contributions are accepted under the same terms (a contributor license agreement
may be requested) to keep that possible.
