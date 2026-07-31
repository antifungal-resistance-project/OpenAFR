# First docking result — VT-1161 redock into 5TZ1 (C. albicans CYP51)

Date: 2026-07-30
Receptor: 5TZ1 chain A + heme (work/receptor.pdbqt)
Ligand:   VT-1161 (oteseconazole), extracted from the crystal
Engine:   AutoDock Vina 1.2.7, exhaustiveness 16, seed 42, box 24 A centered (70.61, 66.28, 4.18)

## Result
| pose | score (kcal/mol) | RMSD vs crystal | N-Fe distance |
|------|------------------|-----------------|---------------|
| 1 (Vina's top) | -12.20 | 2.20 A | 6.21 A (no iron bond) |
| **2** | **-12.12** | **1.53 A** | **2.63 A (iron-bound)** |
| 4 | -11.24 | 1.92 A | 3.26 A |
| others | -10.2 to -11.3 | 2.8 - 12.9 A | > 5 A |

Crystal reference: closest drug nitrogen sits 2.25 A from the heme iron (a real coordination bond).

## Findings
1. Vina CAN reproduce the crystal pose (1.53 A, under the 2.0 A success bar).
2. Vina RANKS it second, behind a pose that misses the iron bond, by 0.08 kcal/mol —
   far below the scoring function's noise floor. Ranking alone is not trustworthy here.
3. Vina's scoring function ignores partial charges: patching the heme iron to +3.000
   produced byte-identical results. Charge-based fixes do not work with Vina.
4. A metal-coordination constraint (azole N within 3 A of heme Fe) filters 9 poses -> 1,
   and that survivor is the correct pose. Domain chemistry beats raw score ranking.

## Implication for the pipeline
The validation gate should not rank on docking score alone. It needs a metal-coordination
filter as a hard constraint. This is deterministic, testable, and provable by construction —
exactly the project's thesis.

## Methodology gotcha (cost an hour, do not repeat)
OpenBabel REORDERS atoms during PDB -> PDBQT conversion. Naive atom-by-atom RMSD gave a
false 6.81 A "failure". Always match atoms by name (or use a graph-matching RMSD).
