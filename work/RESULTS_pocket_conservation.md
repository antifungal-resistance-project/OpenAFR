# Results — cross-species CYP51 pocket conservation (*C. albicans* / *C. auris* / *Aspergillus*)

Date: 2026-08-29. Pre-registered in
[PREREGISTRATION_pocket_conservation.md](PREREGISTRATION_pocket_conservation.md)
(sha256 `bd98a1cd…`, `work/PREREG_pocket_conservation.sha256`, which also pins the four
ortholog sequences). Regenerate:
`conda run -n openafr python scripts/pocket_conservation.py`.

> **This is a shared pocket-conservation MAP, not a gate and not a per-candidate axis.**
> It answers the transferability question a single-species docking score cannot: for each
> residue that lines the modeled azole pocket, is it conserved in the pathogens the mission
> targets, or is it specific to the *C. albicans* stand-in we dock against? It never re-docks
> and never promotes.

## What was tested

The pocket is defined exactly as in the issue-#23 mutation→pocket mapper: a residue lines the
pocket iff a heavy atom falls within **4.5 Å** of the co-crystallized azole (VT1) on the
modeled 5TZ1 chain-A receptor (`work/receptor_A.pdb`). That yields **22 pocket residues**.
Each was aligned (deterministic Needleman–Wunsch, stock BLOSUM62, BLAST-default gaps) across
four CYP51 orthologs, all pinned with SHA-256:

| role | ortholog | source |
|---|---|---|
| reference (structure frame) | *C. albicans* ERG11 | UniProt P10613 |
| verdict | *C. auris* ERG11 | NCBI `XM_085597798.1`, B8441 (translated CDS) |
| verdict | *A. fumigatus* **CYP51A** | UniProt Q4WNT5 |
| supplementary (not scored) | *A. fumigatus* CYP51B | UniProt Q96W81 |

The numbering was **asserted, not assumed**: all 22 structure residues match P10613 at their
positions (the canonical resistance residues V125/F126/Y132/K143 among them), so the
cross-species columns are read in the same frame as every gate.

## The result

Pocket residues sorted by distance to the heme iron (`*` = coordination shell, ≤ 8 Å to iron):

| residue | shell | dFe (Å) | dLig (Å) | auris | CYP51A | CYP51B | verdict |
|---|:--:|---:|---:|:--:|:--:|:--:|---|
| T311 | * | 5.09 | 3.39 | T | S | S | CONSERVATIVE |
| **G307** | * | 6.03 | 3.19 | G | **A** | **A** | **DIVERGENT** |
| G308 | * | 6.21 | 4.16 | G | G | G | IDENTICAL |
| L376 | * | 6.93 | 3.32 | L | I | I | CONSERVATIVE |
| I304 | * | 7.64 | 4.19 | V | L | L | CONSERVATIVE |
| **G303** |  | 8.27 | 3.47 | G | **T** | A | **DIVERGENT** |
| Y132 |  | 8.45 | 3.58 | Y | Y | Y | IDENTICAL |
| I131 |  | 10.09 | 3.73 | I | V | V | CONSERVATIVE |
| F228 |  | 10.09 | 3.38 | F | F | F | IDENTICAL |
| F126 |  | 10.11 | 3.68 | F | F | F | IDENTICAL |
| Y118 |  | 10.52 | 3.74 | Y | Y | Y | IDENTICAL |
| M508 |  | 10.70 | 3.15 | M | L | L | CONSERVATIVE |
| T122 |  | 11.23 | 3.85 | T | T | T | IDENTICAL |
| S378 |  | 11.40 | 3.15 | S | S | S | IDENTICAL |
| H377 |  | 11.70 | 2.79 | H | H | H | IDENTICAL |
| L121 |  | 12.42 | 4.37 | L | L | L | IDENTICAL |
| **F380** |  | 12.99 | 3.40 | F | **M** | I | **DIVERGENT** |
| S507 |  | 14.35 | 4.01 | S | S | S | IDENTICAL |
| F233 |  | 14.47 | 4.28 | F | F | F | IDENTICAL |
| P230 |  | 15.73 | 4.09 | P | P | P | IDENTICAL |
| Y505 |  | 16.94 | 4.27 | Y | Y | Y | IDENTICAL |
| Y64 |  | 17.16 | 3.33 | Y | Y | Y | IDENTICAL |

**Tally: 14 IDENTICAL, 5 CONSERVATIVE, 3 DIVERGENT.** Coordination shell = 5 residues, of
which **one (G307) is divergent**. Divergent pocket positions: **303, 307, 380**.

## What this shows

- **The *C. albicans* → *C. auris* pocket is near-invariant — 21 / 22 identical, the lone
  difference conservative.** Every pocket-lining residue except one is the *same* amino acid in
  the mission's actual target (*C. auris* ERG11) as in the 5TZ1 stand-in we dock against; the
  single difference is **I304V** — a conservative branched-aliphatic swap (BLOSUM62 I/V = +3),
  in the coordination shell at 7.6 Å from the iron. This is the decisive result: the shortlist
  is ranked on a pocket that is a **faithful residue-level proxy** for *C. auris*, so its
  geometry-based transfer to the target species is structurally underwritten at essentially
  every contact. Note the contrast with **global** identity — only 70.6% across the whole
  524-aa protein: the two Candida-clade enzymes diverge substantially *outside* the drug site
  and barely at all *inside* it. (Sequence identity is necessary, not sufficient — the backbone
  is not re-docked — but a conserved contact residue is the precondition for transfer, and here
  it holds at 21/22 with the 22nd conservative.)
- **All *divergent* transferability risk is in the *Aspergillus* (broad-spectrum) direction.**
  The three divergent residues — G303, G307, F380 — differ only in *A. fumigatus*, never between
  the two Candida enzymes (whose sole pocket difference, I304V, is conservative). So the map
  draws a clean line: an albicans/auris-modeled geometry is a strong *Candida* hypothesis and a
  **weaker, flagged** hypothesis for the mold.
- **One divergence sits in the coordination shell: G307 (Gly→Ala in both *Aspergillus*
  paralogs), 3.19 Å from the azole and 6.0 Å from the iron.** Because the shell governs the
  N→Fe anchor that *every* candidate shares, this is a candidate-independent flag on the
  broad-spectrum claim: added β-carbon bulk immediately beside the coordinating channel in
  *A. fumigatus* is exactly the kind of change that can perturb a geometry optimised against a
  glycine there. G303 and F380 are peripheral (> 8 Å) and can only matter for candidates whose
  pose reaches that far.

## Scope — what needs the off-box poses

This map flags **residues** and the **shared coordination shell**; it does not claim which
individual candidate's pose contacts G303/G307/F380, because the per-candidate pose trees are
gitignored and regenerate off-box (same constraint as the resistance and pose-convergence
axes). The candidate-independent conclusions above (auris pocket invariant; G307 shell
divergence toward *Aspergillus*) hold with no pose file. Attributing a specific divergent
contact to a specific candidate — e.g. "candidate X leans on F380, so its mold hypothesis is
weaker" — is the documented refinement, to be run when the candidate poses are on the same host.

## Bottom line

The single-species docking target is **not** a liability for the mission's core claim: at the
residue level the modeled pocket is essentially the *C. auris* pocket (21/22 identical, the 22nd
a conservative I304V), so *C. auris* transfer is underwritten and needs no separate re-dock to
defend. Broad-spectrum (*A. fumigatus*) transfer is the honest caveat — three divergent pocket
residues, one of them (G307) beside the coordinating channel — and is now mapped and flagged
rather than assumed.
