# C. auris FKS1 (GSC1) reference — provenance

The reference the FKS1 re-caller (`openafr/fks1_caller.py`, `scripts/recall_fks1.py`)
translates isolate consensus sequences against, for the echinocandin/FKS1 early-warning
track (v2). A resistance call is a statement *relative to a reference*; that reference must
be pinned and hashable, exactly like `protocol.yaml`, the ERG11 reference, and the snapshot
`INDEX.tsv` hashes. A call you cannot reproduce against a named reference is not a call.

## The sequence

- **File:** `fks1_cds.fasta`
- **Accession:** `XM_085597048.1` (NCBI RefSeq mRNA), protein `XP_085453772.1`.
- **Organism:** *Candidozyma auris* (formerly *Candida auris*), strain **B8441** — the
  clade I reference-genome strain (Pakistan); same strain + RefSeq annotation batch as the
  pinned ERG11 reference, so the two calls share a consistent coordinate frame.
- **Feature:** gene **`GSC1`** (the standard FKS1 / glucan-synthase alias),
  locus_tag **`B9J08_02922`**, GeneID 312328390, 1,3-β-D-glucan synthase catalytic subunit
  (the echinocandin target). Domains: PFAM PF02364 (1,3-β-glucan synthase) + PF14288
  (FKS1 domain).
- **CDS:** 5667 nt → **1888-aa** protein + stop.
- **CDS SHA-256:** `d50ececa0ff8977b61fce9422d2fc255a5e1d97144dfe6fcb24bb6501fbf15ee`
  (the bare `ACGT` string, no header/newlines — what the re-caller hashes at load).

Rebuild + verify byte-identically:

```bash
python3 - <<'PY'
import urllib.request, hashlib
u="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
u+="?db=nuccore&id=XM_085597048.1&rettype=fasta_cds_na&retmode=text"
fa=urllib.request.urlopen(u, timeout=90).read().decode()
seq="".join(l.strip() for l in fa.splitlines() if not l.startswith(">")).upper()
assert hashlib.sha256(seq.encode()).hexdigest() == \
    "d50ececa0ff8977b61fce9422d2fc255a5e1d97144dfe6fcb24bb6501fbf15ee"
print("OK", len(seq))
PY
```

## Why this reference, and why the numbering is trustworthy

Echinocandin resistance in *C. auris* is caused by point substitutions in two short,
conserved FKS1 **hot-spot** regions. The literature reports the dominant HS1 mutation as
**S639F / S639P / S639Y** and an HS2 residue at **R1354**. "S639" and "R1354" are *claims
about the reference frame* — exactly as "Y132" was for ERG11 — so we do **not** assert them,
we **verify** them by translating this CDS:

| hot-spot | residues | motif (translated) | anchor | panel |
|---|---|---|---|---|
| **HS1** | 635–643 | `FLTLSLRDP` | **S639** (residue 639 = Ser ✓) | S639F / S639P / S639Y |
| **HS2** | 1350–1358 | `DWIRRYTLS` | **R1354** (residue 1354 = Arg ✓) | *reported, not yet tagged* |

Both literature positions land on their stated wild-type letters. `load_reference()` and
`tests/test_fks1_caller.py` re-check both anchors on every run, so a reference swap that
breaks the numbering fails loudly rather than silently mis-calling.

Protein N-terminus (sanity): `MSYDNNHNYYDPNQQQGGQPGGEYYQQGAY…` — the canonical Fks1
cytoplasmic N-terminus (no signal peptide).

### Identity was verified, not assumed — a wrong-gene guess was caught here

An early draft guessed FKS1 = locus `B9J08_004112` (`B9J08_04866`, `XM_085598933.1`). That
locus carries a glucan-synthase InterPro domain but is **not Fks1**: it has a signal-peptide
N-terminus (`MKLASLAAVSLAVAGAVAK`), a Ser/Thr/Pro low-complexity "hot-spot" region, residue
639 = Ile, and **no `…SLRDP` HS1 motif**. Pinning by locus tag alone would have anchored the
caller to the wrong gene. The correct FKS1 was confirmed by its N-terminus, its FLTLSLRDP HS1
motif, and the S639/R1354 numbering check above.

## Scope / honesty

- This is the **clade I** reference. Isolates from other clades carry clade-defining FKS1
  polymorphisms that are *not* resistance mutations; the re-caller reports every
  non-synonymous difference within a hot-spot window, and the emergence detector — not the
  re-caller — decides which are novel/rising. Only the known-resistance panel (S639F/P/Y) is
  tagged.
- **HS2 is reported but not yet panel-tagged.** The HS2 window anchor (R1354) is
  numbering-verified, but the specific HS2 *mutant* letters are not yet confidently pinned
  from the C. auris literature, so HS2 substitutions are emitted verbatim/untagged until a
  literature pin fixes the mutant set. Widening the panel is then a one-line data edit in
  `FKS1_WINDOWS` — no downstream code change. This is a documented gap, not a silent drop.
- **Detection only.** Echinocandins do not coordinate a metal; the CYP51/heme-iron geometry
  "so-what" (`openafr/structural.py`) does not transfer to FKS1, and no structural verdict is
  claimed for an FKS1 call. See `.context/fks1_detection_scope.md` §6.
- Non-substitution echinocandin resistance (FKS1 promoter/expression changes, non-FKS1 loci)
  is out of scope for this CDS hot-spot re-caller — the same honest boundary the ERG11
  re-caller draws for TR-type calls.
