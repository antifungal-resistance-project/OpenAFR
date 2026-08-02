# OpenAFR — a self-validating antifungal screening pipeline

An open-source computational pipeline for finding candidate antifungal molecules against
drug-resistant *Candida auris*, built around one rule:

> **The pipeline must prove it can re-discover the drugs we already know work, on molecules
> it has never seen, before it is allowed to rank anything new.**

If it cannot, the validation gate exits nonzero and blocks the screen.

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

**The finding:** for azoles against CYP51, *how closely a molecule brings a nitrogen to the
heme iron* discriminates real drugs from look-alike decoys far better than the docking
program's own affinity score — on the exact same poses. Geometry beats the scoring function.

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

protocol.yaml                     the frozen run protocol (box, search effort, pass bar), hash-pinned

openafr/pdbqt.py                  shared Vina/PDBQT/receptor parsers
openafr/scoring.py                metrics() / report() — the ranking math
openafr/protocol.py               loads + hash-verifies protocol.yaml

scripts/prep_receptor.py          5TZ1.pdb -> work/receptor.pdbqt, chain-A+heme, matches the anchor
scripts/prep_ligands.py           SMILES -> 3D, deterministic, failures logged not dropped
scripts/make_decoys.py            property-matched decoy generation from ChEMBL
scripts/screen.sh                 the docking engine — one fixed protocol, all cores, resumable
scripts/run_screen2.sh            thin wrapper: run-2 held-out validation via screen.sh
scripts/rank_candidates.py        screen poses -> candidate list ranked by the validated criterion
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

## Design notes worth knowing

Two traps this pipeline is explicitly built to avoid, both documented in `work/`:

1. **Circular decoys.** The gate requires a nitrogen near the heme iron. If decoys had no
   aromatic nitrogen, the gate would separate actives from decoys perfectly while really only
   asking "does this molecule contain a nitrogen?" Every decoy is therefore required to be a
   plausible iron-coordinator.
2. **Size shortcuts.** Docking score correlates with heavy-atom count at r = −0.91 on this
   system. Unmatched decoys would let the pipeline "succeed" by preferring heavy molecules.
   Decoys are matched on MW, cLogP, heavy atoms, and rotatable bonds.

## Status

Prove-the-core spike: complete and passing. The `screen/` stage now exists as a single,
parallel docking engine (`scripts/screen.sh`) that docks any prepared ligand library under the
one frozen protocol; `scripts/rank_candidates.py` turns its poses into a candidate list ranked
by the validated iron-approach criterion. Run-2 validation now goes through the same engine, so
validation and screening can never drift apart. `filter/` and `report/` stages are next.

No wet-lab collaborator yet — that is the critical open dependency, since nothing here is
validated until someone puts candidates on real fungus. And screening a novel library is only
justified *after* `validate_gate2.py` passes: a high rank is a hypothesis for a wet lab, never
a hit.

To screen a novel library once the gate has passed:

```bash
python scripts/prep_ligands.py my_library.smi work/screen_ligands   # SMILES -> 3D
./scripts/screen.sh                                                  # dock all cores -> work/screen/
python scripts/rank_candidates.py work/screen -o work/candidates.tsv # ranked by validated criterion
```

Receptor preparation is now scripted end to end (`scripts/prep_receptor.py`), so the
reproduce recipe runs from a fresh clone: `data/structures/5TZ1.pdb` -> `work/receptor.pdbqt`,
with the extracted chain-A + heme atoms asserted identical to the tracked `work/receptor_A.pdb`
the gate grades against. Verified by re-docking VT-1161 (top −12.1 kcal/mol; iron-bound pose
at N–Fe 2.63 Å, matching `work/RESULTS_redock_VT1.md`).

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
