# CYP51 orthologs — provenance (issue #76)

The pinned CYP51 (sterol 14α-demethylase / ERG11) protein sequences the
cross-species pocket-conservation map (`openafr/conservation.py`,
`scripts/pocket_conservation.py`) aligns to the modeled *C. albicans* pocket.
A conservation verdict is a statement *relative to named, pinned sequences*; like
`protocol.yaml` and the ERG11 re-caller reference, they must be hashable so the map
regenerates byte-for-byte with no network at run time.

## The sequences

| file | accession | organism / gene | length | role |
|---|---|---|---:|---|
| `P10613.fasta` | UniProt **P10613** (CP51_CANAL, SV=2) | *Candida albicans* SC5314 ERG11 | 528 aa | **reference** — the 5TZ1 structure frame |
| `Q4WNT5.fasta` | UniProt **Q4WNT5** (CP51A_ASPFU) | *Aspergillus fumigatus* Af293 **cyp51A** | 515 aa | verdict ortholog (mold azole target + resistance locus) |
| `Q96W81.fasta` | UniProt **Q96W81** (cyp51B) | *Aspergillus fumigatus* **cyp51B** | 524 aa | supplementary (reported, **not** scored) |

The fourth ortholog in the verdict set — *C. auris* (Candidozyma auris) ERG11 — is
**not** duplicated here. It is translated at run time from the already-pinned reference
CDS `data/earlywarning/erg11_reference/erg11_cds.fasta` (NCBI RefSeq `XM_085597798.1`,
strain B8441), the single source of truth the ERG11 re-caller also uses. Its 524-aa
protein carries V125/F126/Y132/K143 at the C. albicans-frame positions (the canonical
azole-resistance residues), confirming the two Candida-clade sequences share numbering.

## SHA-256 (raw file bytes — what the pre-registration pins)

```
82621a5bfe871807a9cebeb424afd6ad1e59c2d53449645b8ecab65af20911a3  P10613.fasta
2d8a8e15e4eae09ca390c86bb4e0a94c21239cdcc911d615b5bb856e20b77eb5  Q4WNT5.fasta
87e76a70c6d9ee98b21500030fcc1c6c1535d8fe46a6571ea132695bf7f51810  Q96W81.fasta
```

## Rebuild + verify byte-identically

```bash
for acc in P10613 Q4WNT5 Q96W81; do
  curl -sS "https://rest.uniprot.org/uniprotkb/${acc}.fasta" -o data/cyp51_orthologs/${acc}.fasta
done
shasum -a 256 data/cyp51_orthologs/*.fasta   # must match the block above
```

Retrieved 2026-08-29 from the UniProt REST API (`rest.uniprot.org/uniprotkb/<acc>.fasta`).

## Why these accessions

- **P10613** is *the* reviewed *C. albicans* ERG11 entry, and its numbering is the frame
  the 5TZ1 crystal structure (the docking receptor) uses — verified: structure residues
  V125/F126/Y132/K143 match P10613 at those positions (`assert_reference_frame`). Anchoring
  the map on this sequence keeps the pocket residues in the same coordinate frame as every
  gate and the frozen protocol.
- **Q4WNT5 (CYP51A)** is the clinically dominant azole target and resistance locus in
  *A. fumigatus* (TR34/L98H, TR46/Y121F sit on this paralog), so it — not CYP51B — carries
  the drug-target verdict for the broad-spectrum (mold) direction.
- **Q96W81 (CYP51B)** is the second *A. fumigatus* paralog; it is reported as a supplementary
  column for completeness but does not enter the conserved/divergent verdict.
