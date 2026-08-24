# Front-door dogfood — the end-to-end pipeline runs, on a laptop

Date: 2026-08-23
Command chain: `scripts/prep_receptor.py` → `scripts/score_molecule.py`
Host: Apple Silicon (Darwin arm64), `openafr` conda env
Engine: AutoDock Vina 1.2.7, Open Babel 3.1.0, frozen protocol (box 26 Å @ 70.61/66.28/4.18,
exhaustiveness 32, seed 42, 20 modes)

## Why this run exists

`score_molecule.py` (the SMILES → geometry front door, PR #58) chains four already-validated
stages, but the **subprocess dock** in the middle had never been executed end-to-end — only
its decision/interpretation logic was unit-tested (13 known-answer tests), the same discipline
the ERG11 re-caller used before its real run. This closes that gap: the whole chain now runs,
and it gives the right known-answers.

## Key finding: docking needs no cloud VM

The re-caller's reads→consensus step is Linux/x86-only (bioconda `sra-tools`/`minimap2`/
`samtools`), which is why the runbook provisions a VM. **Docking is different:** its tools
(`vina`, `openbabel`) are conda-forge and install on Apple Silicon, so the full SMILES →
geometry pipeline runs locally in ~9 s wall (≈63 s CPU across 12 jobs) per molecule. The
"needs a VM" language in `score_molecule.py` referred to the recaller toolchain; the docstring
and PENDING-path text were corrected in this change to say "the pinned openafr env," not "VM."

## Receptor prep (reproducible from script)

    conda run -n openafr python scripts/prep_receptor.py
    # extracted chain A: 3971 atoms (3928 protein, 43 heme)
    # anchor check: 3971 atoms match work/receptor_A.pdb exactly
    # heme Fe at (64.163, 71.316, 2.711)  [matches anchor]
    # output: work/receptor.pdbqt  (sha256 d1f52a252f8f...)

## Known-answer results

Two real azoles (both in the validated panel) and one non-azole control:

| molecule     | triage        | docked | closest N–Fe | iron-bound poses | Vina (secondary) | verdict |
|--------------|---------------|--------|--------------|------------------|------------------|---------|
| voriconazole | near-known    | yes    | **2.689 Å**  | 1/20             | −9.3 kcal/mol    | within validated band (2.47–2.88 Å) |
| fluconazole  | near-known    | yes    | **2.715 Å**  | 4/20             | −8.2 kcal/mol    | within validated band (2.47–2.88 Å) |
| aspirin      | out-of-scope  | **no** | —            | —                | —                | no aromatic N → criterion undefined; **not docked** |

Both real drugs coordinate the heme iron inside the band the validated actives occupy —
exactly the expected answer. Aspirin is stopped at triage and never docked, proving the honest
gate that refuses to compute a geometry the method makes no claim about.

## What this does and does NOT establish

- **Does:** the four-stage chain executes end-to-end from a SMILES string, the docking
  subprocess link works, and the output is sane on known-answer inputs. The tool runs locally.
- **Does NOT:** re-adjudicate the science. Every standing caveat still holds — a tight N–Fe is
  a hypothesis for a wet lab, not a hit; coordinating the iron is necessary but NOT sufficient
  (property-matched decoys reach it too, the decoy ceiling); the durable validated metric is
  AUC ≈ 0.79 (the top-1% enrichment was a single-run artifact). No novel candidate was scored.

## Reproduce

    conda run -n openafr python scripts/prep_receptor.py
    conda run -n openafr python scripts/score_molecule.py \
        --smiles "C[C@@H](C1=NC=NC=C1F)[C@](CN2C=NC=N2)(C3=C(C=C(C=C3)F)F)O" --name voriconazole

Outputs (`work/receptor.pdbqt`, `work/score/`) are gitignored — regenerate them from the
commands above.
