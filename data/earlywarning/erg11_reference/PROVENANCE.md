# C. auris ERG11 reference — provenance

The reference the ERG11 re-caller (`openafr/recaller.py`, `scripts/recall_erg11.py`)
translates isolate consensus sequences against. A resistance call is a statement
*relative to a reference*; that reference must be pinned and hashable, exactly like
`protocol.yaml` and the snapshot `INDEX.tsv` hashes. A call you cannot reproduce
against a named reference is not a call.

## The sequence

- **File:** `erg11_cds.fasta`
- **Accession:** `XM_085597798.1` (NCBI RefSeq mRNA)
- **Organism:** *Candidozyma auris* (formerly *Candida auris*), strain **B8441**
  — the clade I reference-genome strain (Pakistan).
- **Feature:** `ERG11` / `B9J08_03698`, sterol 14α-demethylase (CYP51),
  protein `XP_085452763.1`.
- **CDS:** 1575 nt, location 1..1575 → 524-aa protein + stop.
- **CDS SHA-256:** `764e6d1a429dcfda48e8904d8b438425aee42833d5fe1a354743e07bf629ce11`
  (the bare `ACGT` string, no header/newlines — what the re-caller hashes at load).

Rebuild + verify byte-identically:

```bash
python3 - <<'PY'
import urllib.request, hashlib
u="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
u+="?db=nuccore&id=XM_085597798.1&rettype=fasta_cds_na&retmode=text"
fa=urllib.request.urlopen(u, timeout=60).read().decode()
seq="".join(l.strip() for l in fa.splitlines() if not l.startswith(">"))
assert hashlib.sha256(seq.encode()).hexdigest().startswith("764e6d1a")
print("OK", len(seq))
PY
```

## Why this reference, and why the numbering is trustworthy

The canonical *C. auris* azole-resistance substitutions are reported in the
literature against the B8441 Erg11 coordinate frame (Lockhart et al. 2017;
Chowdhary et al.). The re-caller must emit tokens (`Y132F`, `K143R`, `F126L`, …)
in that same frame or the downstream mapping (`openafr.mapping`, which asserts the
wild-type letter) will reject them.

We do not assert that frame — we **verify** it. Translating this CDS gives the
wild-type residues the panel names its mutations *from*:

| position | this reference | panel mutation |
|---:|:--:|:--|
| 125 | **V** | V125A (part of the VF125AL clade-III haplotype) |
| 126 | **F** | F126L |
| 132 | **Y** | Y132F |
| 143 | **K** | K143R |

All four match. `tests/test_recaller.py` re-checks these on every run, so a
reference swap that breaks the numbering fails loudly rather than silently
mis-calling.

Protein N-terminus (sanity): `MALKDCIVDVVDRFSALPVP…`.

## Scope / honesty

- This is the **clade I** reference. Isolates from other clades carry
  clade-defining ERG11 polymorphisms that are *not* resistance mutations; the
  re-caller reports every non-synonymous difference from this reference, and the
  emergence detector (`openafr.emergence`) — not the re-caller — decides which are
  novel/rising. The re-caller tags only the known-resistance panel positions;
  everything else is emitted verbatim, never suppressed and never asserted to be
  resistance.
- Tandem-duplication / promoter calls (TR-type) are **not** representable as point
  substitutions and are out of scope for this CDS-level re-caller — the same honest
  boundary `openafr.mapping` already draws for non-substitution tokens.
