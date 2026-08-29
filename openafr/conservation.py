"""Cross-species CYP51 pocket conservation (issue #76).

The whole shortlist is docked against ONE species' CYP51 -- the *C. albicans*
5TZ1 chain-A stand-in (`work/receptor_A.pdb`), the same scaffold the frozen
protocol, the Y132F resistance port and every gate share. But the mission targets
*C. auris* and (for the broad-spectrum claim) *Aspergillus fumigatus*. If a residue
that lines the modeled pocket is DIFFERENT in those species, an albicans-optimized
geometry may not transfer -- and a candidate whose binding leans on that contact is
a weaker cross-species hypothesis than one anchored only on conserved residues.

This module answers, structurally and from pinned sequences, the question a docking
score against one species can never answer:

  For each residue that lines the modeled azole pocket, is it CONSERVED across the
  pathogens we care about, or is it species-specific?

It is deliberately built on the SAME pocket definition as `openafr.mapping` (issue
#23): a residue "lines the pocket" if any of its heavy atoms comes within
`CONTACT_CUTOFF` (4.5 A) of any atom of the co-crystallized azole (VT1 from 5TZ1,
`work/ref_VT1.pdb`). So "pocket residue" means the same thing here as in the
mutation->pocket mapper -- computed from the structure, never hand-picked.

Honesty constraints (mirroring openafr.mapping / openafr.recaller)
------------------------------------------------------------------
1. WE ASSERT THE NUMBERING, WE DO NOT ASSUME IT. Every pocket residue's identity in
   the *structure* must match the pinned *C. albicans* reference sequence at that
   position (the 5TZ1 frame IS the P10613 frame -- verified for V125/F126/Y132/K143).
   A single mismatch means the structure and the sequence are in different frames and
   the whole cross-species read is void; `assert_reference_frame` raises rather than
   silently aligning the wrong columns.

2. THE ALIGNMENT IS DETERMINISTIC AND PINNED. Conservation is read off a
   Needleman-Wunsch global alignment (BLOSUM62, affine gaps) of each ortholog to the
   *C. albicans* reference -- no external MSA binary, no network at run time. The
   sequences are committed with provenance + SHA-256 (`data/cyp51_orthologs/`,
   `data/earlywarning/erg11_reference/`), so the map regenerates byte-for-byte.

3. THE VERDICT SET IS THE DRUG TARGETS, STATED UP FRONT. The conserved/divergent
   verdict is computed over the three azole DRUG-TARGET orthologs -- *C. albicans*
   ERG11, *C. auris* ERG11, *A. fumigatus* CYP51A. The second *A. fumigatus* paralog
   CYP51B is reported as a supplementary column but does NOT enter the verdict:
   CYP51A is the clinically dominant azole target and resistance locus in the mold.

Not a gate, and not a per-candidate axis. It produces a shared pocket-conservation
MAP; which specific residues a given candidate's pose contacts needs the per-candidate
pose files (gitignored, regenerated off-box), so the map flags the pocket residues and
-- decisively -- the heme-coordination shell that EVERY candidate's N->Fe anchor shares.
"""
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import mapping
from openafr.pdbqt import iron_position

CONTACT_CUTOFF = mapping.CONTACT_CUTOFF   # 4.5 A -- same pocket definition as issue #23
COORD_SHELL_A = 8.0                       # <= this to the heme iron = coordination shell

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

# Standard genetic code (translation of the pinned C. auris ERG11 CDS).
_CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# BLOSUM62 (the 20 standard residues). Symmetric; queried both orders.
_B62_ORDER = "ARNDCQEGHILKMFPSTWYV"
_B62_ROWS = [
    [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],
    [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
    [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],
    [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],
    [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],
    [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],
    [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],
    [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],
    [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],
    [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],
    [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],
    [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],
    [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],
    [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2],
    [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],
    [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],
    [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],
    [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],
]
_B62 = {}
for _i, _a in enumerate(_B62_ORDER):
    for _j, _b in enumerate(_B62_ORDER):
        _B62[(_a, _b)] = _B62_ROWS[_i][_j]

GAP_OPEN = -11     # BLAST/EMBOSS BLOSUM62 defaults
GAP_EXTEND = -1


def blosum62(a, b):
    """BLOSUM62 score for two one-letter residues; unknown letters score -4."""
    return _B62.get((a, b), -4)


def translate(cds):
    """Translate a CDS (multiple of 3) to protein, dropping a single trailing stop."""
    aa = []
    for i in range(0, len(cds) - 2, 3):
        aa.append(_CODON[cds[i:i + 3].upper()])
    prot = "".join(aa)
    return prot[:-1] if prot.endswith("*") else prot


def read_fasta(path):
    """Return (header, sequence) of the first record in a FASTA file."""
    lines = pathlib.Path(path).read_text().splitlines()
    header = lines[0].lstrip(">").strip() if lines and lines[0].startswith(">") else ""
    seq = "".join(l.strip() for l in lines if l and not l.startswith(">"))
    return header, seq.upper()


def needleman_wunsch(a, b, gap_open=GAP_OPEN, gap_extend=GAP_EXTEND):
    """Global alignment (BLOSUM62, affine gaps). Returns (aln_a, aln_b) gap-padded.

    Deterministic: ties break toward diagonal, then up, then left, so the same two
    sequences always align identically -- a requirement for a reproducible map.
    """
    n, m = len(a), len(b)
    NEG = float("-inf")
    # Gotoh three-matrix affine-gap DP.
    M = [[NEG] * (m + 1) for _ in range(n + 1)]   # match/mismatch ending
    Ix = [[NEG] * (m + 1) for _ in range(n + 1)]  # gap in b (consume a)
    Iy = [[NEG] * (m + 1) for _ in range(n + 1)]  # gap in a (consume b)
    M[0][0] = 0
    for i in range(1, n + 1):
        Ix[i][0] = gap_open + (i - 1) * gap_extend
    for j in range(1, m + 1):
        Iy[0][j] = gap_open + (j - 1) * gap_extend
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = blosum62(a[i - 1], b[j - 1])
            M[i][j] = s + max(M[i - 1][j - 1], Ix[i - 1][j - 1], Iy[i - 1][j - 1])
            Ix[i][j] = max(M[i - 1][j] + gap_open, Ix[i - 1][j] + gap_extend)
            Iy[i][j] = max(M[i][j - 1] + gap_open, Iy[i][j - 1] + gap_extend)
    # Traceback from the best-scoring end state.
    i, j = n, m
    aln_a, aln_b = [], []
    state = max((M[n][m], "M"), (Ix[n][m], "Ix"), (Iy[n][m], "Iy"))[1]
    while i > 0 or j > 0:
        if state == "M" and i > 0 and j > 0:
            aln_a.append(a[i - 1]); aln_b.append(b[j - 1])
            prev = max((M[i - 1][j - 1], "M"), (Ix[i - 1][j - 1], "Ix"),
                       (Iy[i - 1][j - 1], "Iy"))[1]
            i -= 1; j -= 1; state = prev
        elif state == "Ix" and i > 0:
            aln_a.append(a[i - 1]); aln_b.append("-")
            cont = Ix[i][j] == Ix[i - 1][j] + gap_extend and i - 1 > 0
            i -= 1; state = "Ix" if cont else "M"
        elif state == "Iy" and j > 0:
            aln_a.append("-"); aln_b.append(b[j - 1])
            cont = Iy[i][j] == Iy[i][j - 1] + gap_extend and j - 1 > 0
            j -= 1; state = "Iy" if cont else "M"
        elif i > 0:
            aln_a.append(a[i - 1]); aln_b.append("-"); i -= 1
        else:
            aln_a.append("-"); aln_b.append(b[j - 1]); j -= 1
    return "".join(reversed(aln_a)), "".join(reversed(aln_b))


def reference_column_map(ref, other):
    """Map each 1-based ref position to the aligned residue in `other` ('-' = gap).

    Aligns `other` to `ref` and walks the alignment, so the returned dict is keyed by
    the reference (structure) numbering the pocket residues use.
    """
    aln_ref, aln_other = needleman_wunsch(ref, other)
    col = {}
    rpos = 0
    for cr, co in zip(aln_ref, aln_other):
        if cr != "-":
            rpos += 1
            col[rpos] = co
    return col


def percent_identity(ref, other):
    """Percent identity of `other` to `ref` over aligned (non-gap-in-ref) columns."""
    col = reference_column_map(ref, other)
    if not col:
        return 0.0
    same = sum(1 for p, c in col.items() if c == ref[p - 1])
    return 100.0 * same / len(col)


def pocket_residues(receptor_pdb, ligand_pdb, contact_cutoff=CONTACT_CUTOFF):
    """Pocket-lining residues on the modeled receptor, the openafr.mapping definition.

    Returns a list of dicts sorted by distance to the heme iron:
      {pos, resname, aa1, dist_to_ligand, dist_to_iron, coord_shell}
    A residue lines the pocket iff a heavy atom is within `contact_cutoff` of the
    co-crystal azole; `coord_shell` flags the inner shell (<= COORD_SHELL_A to iron).
    """
    residues = mapping.read_residues(receptor_pdb)
    ligand = mapping.read_ligand_atoms(ligand_pdb)
    fe = iron_position(receptor_pdb)
    out = []
    for pos, (resname, atoms) in residues.items():
        d_lig = mapping._min_dist(atoms, ligand)
        if d_lig is None or d_lig > contact_cutoff:
            continue
        d_iron = mapping._min_dist(atoms, [fe])
        out.append({
            "pos": pos,
            "resname": resname,
            "aa1": AA3TO1.get(resname, "X"),
            "dist_to_ligand": round(d_lig, 3),
            "dist_to_iron": round(d_iron, 3) if d_iron is not None else None,
            "coord_shell": d_iron is not None and d_iron <= COORD_SHELL_A,
        })
    out.sort(key=lambda r: (r["dist_to_iron"] if r["dist_to_iron"] is not None else 1e9))
    return out


def assert_reference_frame(pocket, ref_seq, ref_name="C. albicans ERG11"):
    """Every pocket residue's structure identity must match the reference sequence.

    Guards the numbering: the 5TZ1 structure frame must equal the reference-sequence
    frame or the cross-species columns are read at the wrong positions. Raises
    ValueError listing every mismatch; returns None on success.
    """
    bad = []
    for r in pocket:
        p, aa = r["pos"], r["aa1"]
        seq_aa = ref_seq[p - 1] if 1 <= p <= len(ref_seq) else "?"
        if aa != seq_aa:
            bad.append(f"pos {p}: structure {r['resname']}({aa}) != {ref_name} {seq_aa}")
    if bad:
        raise ValueError(
            "structure/reference numbering mismatch -- cross-species map is void:\n  "
            + "\n  ".join(bad))


def classify(ref_aa, target_aas):
    """Conservation verdict for one pocket residue over the drug-target orthologs.

    `target_aas` is the dict {species_label: aligned_one_letter (or '-')} for the
    verdict set (NOT including the reference itself, NOT including the supplementary
    CYP51B). The rule is frozen in the pre-registration:

      IDENTICAL     -- every target carries exactly `ref_aa` (no gaps).
      CONSERVATIVE  -- not identical, no gaps, and every substitution is positively
                       scoring under BLOSUM62 (score >= 1).
      DIVERGENT     -- any target has a gap/deletion, or a non-positive BLOSUM62
                       substitution (score <= 0).
    """
    if any(a == "-" for a in target_aas.values()):
        return "DIVERGENT"
    if all(a == ref_aa for a in target_aas.values()):
        return "IDENTICAL"
    subs = [a for a in target_aas.values() if a != ref_aa]
    if all(blosum62(ref_aa, a) >= 1 for a in subs):
        return "CONSERVATIVE"
    return "DIVERGENT"


def build_map(receptor_pdb, ligand_pdb, ref_seq, verdict_orthologs, supplementary=None,
              contact_cutoff=CONTACT_CUTOFF):
    """Assemble the cross-species pocket-conservation map.

    `ref_seq`             -- the C. albicans reference protein (structure frame).
    `verdict_orthologs`   -- {label: sequence} that DEFINE the conserved/divergent
                             verdict (C. auris ERG11, A. fumigatus CYP51A).
    `supplementary`       -- {label: sequence} reported but NOT scored (CYP51B).

    Returns {pocket, summary}. Raises if the structure/reference frame disagrees.
    """
    supplementary = supplementary or {}
    pocket = pocket_residues(receptor_pdb, ligand_pdb, contact_cutoff)
    assert_reference_frame(pocket, ref_seq)

    cols = {lab: reference_column_map(ref_seq, seq)
            for lab, seq in {**verdict_orthologs, **supplementary}.items()}

    rows = []
    for r in pocket:
        p, ref_aa = r["pos"], r["aa1"]
        verdict_aas = {lab: cols[lab].get(p, "-") for lab in verdict_orthologs}
        suppl_aas = {lab: cols[lab].get(p, "-") for lab in supplementary}
        verdict = classify(ref_aa, verdict_aas)
        rows.append({
            **r,
            "ref_aa": ref_aa,
            "orthologs": verdict_aas,
            "supplementary": suppl_aas,
            "verdict": verdict,
            "transferability_risk": verdict == "DIVERGENT",
        })

    shell = [r for r in rows if r["coord_shell"]]
    summary = {
        "n_pocket": len(rows),
        "n_identical": sum(1 for r in rows if r["verdict"] == "IDENTICAL"),
        "n_conservative": sum(1 for r in rows if r["verdict"] == "CONSERVATIVE"),
        "n_divergent": sum(1 for r in rows if r["verdict"] == "DIVERGENT"),
        "n_coord_shell": len(shell),
        "coord_shell_divergent": sorted(
            r["pos"] for r in shell if r["verdict"] == "DIVERGENT"),
        "divergent_positions": sorted(r["pos"] for r in rows if r["verdict"] == "DIVERGENT"),
        "identity_to_albicans": {
            lab: round(percent_identity(ref_seq, seq), 1)
            for lab, seq in {**verdict_orthologs, **supplementary}.items()
        },
    }
    return {"pocket": rows, "summary": summary}
