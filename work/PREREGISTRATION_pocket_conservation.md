# Pre-registration — cross-species CYP51 pocket conservation (issue #76)

**Written 2026-08-29.** Everything below is fixed in advance; deviations are recorded as
deviations, not edited in. This axis inherits the modeled receptor and the pocket
definition of the mutation→pocket mapper — read [../openafr/mapping.py](../openafr/mapping.py)
(issue #23) first. It shares 5TZ1 chain A (`work/receptor_A.pdb`), the same heme, iron and
coordinate frame as the frozen protocol and every gate.

## Background — why this run exists

The whole shortlist is ranked, vetted and demoted against **one** species' CYP51: the
*C. albicans* 5TZ1 chain-A stand-in. The mission targets ***C. auris*** and, for any
broad-spectrum claim, **_Aspergillus fumigatus_**. If a residue that *lines* the modeled
azole pocket is different in those species, an albicans-optimized geometry may not transfer,
and a candidate whose binding leans on that contact is a weaker cross-species hypothesis.
No docking score against a single species can see this; it is a fact about the sequences and
the structure.

This is a **shared pocket-conservation MAP, not a gate and not a per-candidate axis.** It
never re-docks, never promotes, invents no composite score. It answers one question per
pocket residue — conserved or species-specific — and reports the heme-coordination shell that
governs every candidate's N→Fe anchor.

## No free parameter (why this is pre-registerable at all)

Every number below is **inherited or stock**, chosen before any conservation verdict was read
and with nothing to tune toward a desired outcome:

- **Pocket definition:** `openafr.mapping.CONTACT_CUTOFF` = **4.5 Å** — a residue lines the
  pocket iff a heavy atom is within 4.5 Å of any atom of the co-crystallized azole
  (`work/ref_VT1.pdb`, the VT1 ligand from 5TZ1). This is the *identical* definition the
  issue-#23 mutation→pocket mapper already uses; it is not re-picked here.
- **Substitution scoring:** stock **BLOSUM62** with BLAST/EMBOSS default affine gaps
  (open −11, extend −1). No custom matrix, no tuned gap penalty.
- **Alignment:** deterministic Needleman–Wunsch (Gotoh affine gaps), ties broken
  diagonal→up→left, so each ortholog aligns to the reference identically every run. No
  external MSA binary, no network at run time.

## Fixed inputs (pinned, hashable)

- **Receptor / pocket source:** `work/receptor_A.pdb` (5TZ1 chain A) + `work/ref_VT1.pdb`
  (co-crystal azole), the same pair `openafr.mapping` reads. Heme iron from
  `openafr.pdbqt.iron_position`.
- **Reference sequence (structure frame):** *C. albicans* ERG11, UniProt **P10613**
  (`data/cyp51_orthologs/P10613.fasta`, sha256 `82621a5b…`). The 5TZ1 frame **is** the P10613
  frame — asserted, not assumed (see honesty constraint 1).
- **Verdict orthologs** (the drug targets that DEFINE the conserved/divergent verdict):
  - *C. auris* (Candidozyma auris) ERG11 — translated from the pinned B8441 reference CDS
    `data/earlywarning/erg11_reference/erg11_cds.fasta` (NCBI `XM_085597798.1`,
    sha256 `9ffff0a4…`), the same reference the ERG11 re-caller uses.
  - *A. fumigatus* **CYP51A**, UniProt **Q4WNT5** (`data/cyp51_orthologs/Q4WNT5.fasta`,
    sha256 `2d8a8e15…`) — the clinically dominant azole target/resistance locus in the mold.
- **Supplementary ortholog (reported, NOT scored):** *A. fumigatus* **CYP51B**, UniProt
  **Q96W81** (`data/cyp51_orthologs/Q96W81.fasta`, sha256 `87e76a70…`).

## Frozen procedure

1. Extract the pocket residues on `work/receptor_A.pdb` (heavy atom ≤ 4.5 Å of the VT1
   ligand); record each residue's closest approach to the heme iron.
2. **Assert the frame:** every pocket residue's structure identity must equal P10613 at that
   position. Any mismatch voids the run (raise), because the cross-species columns would be
   read at the wrong positions.
3. Globally align each ortholog to P10613; read off, per pocket position, the aligned residue
   in each ortholog ('-' = gap/deletion).
4. Classify each pocket residue over the **verdict set** {C. albicans, C. auris, A. fumigatus
   CYP51A}:
   - **IDENTICAL** — every verdict ortholog carries exactly the reference residue.
   - **CONSERVATIVE** — not identical, no gaps, and *every* substitution scores BLOSUM62 ≥ 1.
   - **DIVERGENT** — any verdict ortholog has a gap, or a substitution scoring BLOSUM62 ≤ 0.

## Pre-committed reporting + interpretation rules

- A **DIVERGENT** pocket residue is flagged a **transferability risk**.
- The **coordination shell** = pocket residues within **8.0 Å** of the heme iron (a reporting
  split, not a verdict input). A divergence *in the shell* is the strongest flag: the shell
  governs the N→Fe anchor **every** candidate shares, so a shell divergence is a
  candidate-independent transferability concern; a peripheral divergence can only affect
  candidates whose pose reaches that far.
- **Per-candidate resolution is explicitly out of scope here.** Which residues a given
  candidate's pose actually contacts needs the per-candidate pose trees (gitignored,
  regenerated off-box). This map therefore flags residues and the shared shell; attributing a
  specific divergent contact to a specific candidate is the documented refinement, to be run
  when the pose files are on the same host.

## What a result cannot claim (limitations, fixed in advance)

- Sequence identity at a pocket position is **necessary, not sufficient**, for transfer:
  identical residues can still shift under a different backbone (no cross-species structures
  are docked here). A conserved verdict says "the contact residue is the same," not "the
  affinity is the same."
- The map is **not** evidence of broad-spectrum activity; it bounds where an
  albicans-modeled geometry is *structurally plausible* to transfer, and flags where it is not.
- CYP51B is shown for completeness only; a divergence unique to CYP51B does not by itself flag
  a residue (A. fumigatus azole action/resistance runs through CYP51A).

## Deliverables

- `openafr/conservation.py` (the frozen logic) + `scripts/pocket_conservation.py` (CLI).
- `work/candidates/pocket_conservation.json` (the map) + `work/RESULTS_pocket_conservation.md`.
- `tests/test_pocket_conservation.py`.
