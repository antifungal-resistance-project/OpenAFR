# Pre-registration — FKS1 structural proximity (moonshot #4, Rung 2)

> **STATUS: FROZEN 2026-10-17, before the contact map was wired into the live weather
> page.** This document fixes the structure source, alignment method, tier cutoffs, and
> the six things this module refuses to claim. It crosses the FKS1 track's declared
> "detection-only, no structural verdict" boundary by adding a 3D proximity measurement,
> so under this project's proof-before-claim discipline it must declare its scope and
> non-claims in advance (`docs/DIAGNOSTIC_GO_NO_GO.md`). The repro gate
> (`scripts/repro.py`) re-checks this document's SHA-256 via
> `work/PREREG_fks1_structural_proximity.sha256` and fails on any drift.

## Why this exists — Rung 1 is 1D; we can do better with a real structure

Rung 1 (`openafr/fks1_position.py`) gave a sequence-based geography tier for novel FKS1
hot-spot variants: which window, offset from the anchor, nearest known panel residue. This is
strictly one-dimensional (primary sequence). It cannot distinguish a residue that physically
contacts the echinocandin drug from one that is merely sequentially adjacent.

Rung 2 adds the third dimension: for each hot-spot position, the minimum heavy-atom distance
from the aligned *S. cerevisiae* residue to any atom of the bound enfumafungin ligand in the
9WZU cryo-EM structure. Enfumafungin binds the same FKS1 site as clinical echinocandins
(caspofungin, micafungin, anidulafungin); a residue at drug-contact distance in this structure
is at drug-contact distance in the echinocandin binding site.

## Structure and drug

| Field | Value |
|---|---|
| PDB ID | **9WZU** |
| Title | Cryo-EM structure of Fks1 in complex with enfumafungin |
| Organism | *Saccharomyces cerevisiae* S288C |
| Drug | Enfumafungin (comp_id A1E06, chain A / HETATM residue 1906) |
| Resolution | 2.97 Å |
| Source SHA-256 | `1de4c46b84f4caa958a82caa52aab1fff034d0d1b339ea10d62df5119dc45492` |
| Pre-derivation script | `scripts/derive_fks1_contact_map.py` |
| Derived contact map | `data/structures/fks1_contact_map.json` |
| Contact map SHA-256 | `6a2033ac5bd28b821e5cc3a47b759932c5b728c7f85b69958a82f1f89902acea` |

## Alignment method

The *S. cerevisiae* Fks1 sequence in 9WZU and the *C. auris* FKS1 reference protein
(XM_085597048.1, B8441 strain, CDS SHA-256 `d50ececa…`) are aligned with a pure Smith–Waterman
local alignment (inline, stdlib only, no BioPython) over a 50-residue window around each
hot-spot (HS1 635–643, HS3 688–698, HS2 1350–1358). The most-common Sc−Ca offset across
matched pairs in the window is used as the hot-spot offset. Verified offsets:

- HS1: Sc position = Ca position + **4** (no indels in this window)
- HS3: Sc position = Ca position + **4** (no indels in this window)
- HS2: Sc position = Ca position + **3** (one-residue insertion in Ca relative to Sc before HS2)

The alignment is reproducible from the pinned CDS and the pinned 9WZU SHA-256; any change to
either fails `repro.py`.

## Distance metric

For each *C. auris* hot-spot residue at position P, the mapped *S. cerevisiae* residue at
position P + offset is looked up in the structure. The **minimum Euclidean distance** from any
heavy atom (non-hydrogen) of that residue to any heavy atom of the enfumafungin ligand is
recorded as `min_dist_to_drug_A`. Using minimum heavy-atom distance (not Cβ-only) captures the
actual extent of each side chain toward the drug.

## Proximity tiers (pre-registered cutoffs, NOT tuned post-hoc)

| Tier | Cutoff | Meaning |
|---|---|---|
| `drug_contact` | ≤ 5.0 Å | Residue heavy atoms within direct-contact distance of the drug |
| `pocket_adjacent` | 5.0–8.0 Å | Immediately adjacent to the drug contact zone |
| `pocket_proximate` | 8.0–12.0 Å | In the wider binding-site cavity |
| `distal` | > 12.0 Å | Within hot-spot window but not at the drug-contact face |
| `no_structure` | — | Position absent from the contact map (residue not resolved in structure) |

The cutoffs are derived from standard vdW/contact conventions (≤5 Å = direct contact;
≤8 Å = pocket inner shell; ≤12 Å = outer binding-site cavity). They are NOT tuned to
optimize any phenotype classification — phenotype prediction is an explicit non-claim.

## The six things this module DOES NOT claim

1. **No phenotype (R/S) prediction.** Proximity says nothing about whether a substitution at
   this position will confer echinocandin resistance. The only structural claim is distance.
2. **No change to the verdict enum.** `UNCHARACTERIZED_VARIANT` stays; a `drug_contact`
   proximity tier does not convert a variant to `RESISTANCE_MARKER_DETECTED`.
3. **No binding-affinity or ΔΔG claim.** A residue at drug-contact distance does not imply
   weaker or stronger binding — that requires a docking protocol (Rung 3).
4. **No resistance probability.** This is consistent with the concordance-only contract
   (`docs/DIAGNOSTIC_VERDICT_CONTRACT.md`) and the blocked calibration track (#136).
5. **`structural` field for FKS1 stays null.** Proximity lives in
   `position_context[*].proximity`; the `structural` field is reserved for a full
   fit-consequence verdict (requires a docking protocol — deferred to Rung 3).
6. **S. cerevisiae distances are a proxy.** Side-chain conformations may differ between Sc and
   Ca. The HS1/HS3 hot-spot motifs are highly conserved (identical or near-identical residues
   in the aligned window), so the structural proxy is well-founded, but exact Angstrom values
   carry ±1–2 Å uncertainty from species divergence.

## Anti-circularity and sanity checks

The contact zone is derived from the 9WZU cryo-EM drug coordinates, NOT from the C. auris
resistance panel. The panel residues are NOT used to define which positions are "near the drug".

**Pre-committed falsifiable predictions** (verified by `tests/test_fks1_structural.py`,
`test_panel_residue_sites_in_expected_tiers`):
- Ca F635 must land in `drug_contact` (2.94 Å measured → PASS)
- Ca L638 must land in `drug_contact` (2.82 Å measured → PASS)
- Ca S639 must land in `pocket_adjacent` (7.16 Å measured → PASS)
- Ca W691 must land in `drug_contact` (3.66 Å measured → PASS)
- Ca R1354 must land in `pocket_adjacent` (7.53 Å measured → PASS)

These are *independent validations* that the alignment is correct: the resistance panel was
defined from biochemical/clinical data; the cryo-EM structure was solved independently. The
panel residues landing in the drug-contact / adjacent zone is non-circular confirmation.

## Verification

No GCP required. Everything is deterministic:
1. `sha256sum data/structures/fks1_contact_map.json` → must match `6a2033ac…`
2. `python scripts/derive_fks1_contact_map.py --pdb /path/to/9WZU.cif` → regenerates
   the contact map from the pinned source; SHA of 9WZU.cif must match `1de4c46b…`.
3. `python -m pytest tests/test_fks1_structural.py` → 33 tests, all must pass.
4. `python scripts/repro.py` → drift gate includes this prereg's SHA.

## What stays for a future Rung 3

- A cert benchmark: a held-out set of novel FKS1 variants with known phenotype, measuring
  whether `drug_contact` / `pocket_adjacent` tiers predict resistance better than sequence
  proximity alone. This requires GCP (reading + calling a new cohort).
- Populating the `structural` field with a full fit-consequence verdict (requires a docking
  protocol for FKS1 — no heme-iron analogue exists for echinocandins).
