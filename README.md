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

scripts/prep_ligands.py           SMILES -> 3D, deterministic, failures logged not dropped
scripts/make_decoys.py            property-matched decoy generation from ChEMBL
scripts/run_screen*.sh            docking, one fixed protocol for every molecule
scripts/analyze_poses.py          metal-coordination analysis
scripts/validate_gate*.py         the gate — exits nonzero when it cannot recover knowns

work/PREREGISTRATION_run2.md      rules fixed in advance, hash-frozen
work/RESULTS_*.md                 every run, including the failure
```

## Reproducing

```bash
# toolkit (Apple Silicon / Linux)
mamba create -n openafr -c conda-forge python=3.11 rdkit openbabel vina numpy pandas
conda activate openafr

python scripts/prep_ligands.py data/ligands/actives.smi work/ligands
./scripts/run_screen2.sh
python scripts/validate_gate2.py work/screen2 work/receptor_A.pdb   # exit 0 = pass
```

The gate script re-checks the pre-registration hash and refuses to certify a pass if the
rules file was edited after being frozen.

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

Prove-the-core spike: complete and passing. `screen/`, `filter/`, and `report/` stages are
next. No wet-lab collaborator yet — that is the critical open dependency, since nothing here
is validated until someone puts candidates on real fungus.

## License

TBD.
