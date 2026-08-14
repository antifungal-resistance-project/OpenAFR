"""Mutation -> pocket mapping for the C. auris azole early-warning track (issue #23).

Pipeline stage: the *structural* half of the early-warning loop. This is where the
moat starts -- no other genomic surveillance feed does it.

Why this exists
---------------
The emergence detector (openafr.emergence, #22) turns a snapshot diff into a list
of flagged ERG11 (CYP51) substitution tokens -- "Y132F is rising", "F126L is new".
A token is not yet a *structural* statement. This module answers the question that
makes the alert actionable and that a metadata feed can never answer:

  Given a flagged substitution, *where* is that residue in the modeled azole
  pocket, does it touch the drug, and can we build the mutant and dock it?

It reuses the exact modeled receptor the whole project docks against -- the
5TZ1 chain-A *C. albicans* CYP51 scaffold (`work/receptor_A.pdb`), same heme, same
iron, same coordinate frame as the frozen protocol and the *C. auris* Y132F port
(#13) -- so a mapping verdict is in the same structural universe as every gate run.

How the pocket is defined -- data-driven, not a hand-picked radius
-----------------------------------------------------------------
A residue "lines the pocket" if any of its heavy atoms comes within `contact_cutoff`
(default 4.5 A, a standard vdW/contact distance) of any atom of the co-crystallized
azole in `work/ref_VT1.pdb` (the VT1 ligand from 5TZ1). This is the textbook
definition of a binding-site residue -- computed from the structure at map time, not
asserted -- so "Y132 contacts the drug" is measured, not claimed. Each verdict also
reports the residue's closest approach to the catalytic heme iron, the atom azoles
coordinate.

The three honesty constraints (mirrors emergence.py / mutate_receptor.py)
------------------------------------------------------------------------
1. WE ASSERT THE WILD TYPE, WE DO NOT ASSUME IT. A token's wt letter must match the
   residue actually present at that position in the modeled structure. If residue 132
   is not TYR, the mapping is refused with the observed residue reported -- a
   numbering or structure mismatch is surfaced, never papered over.

2. NON-RESIDUE TOKENS ARE REFUSED, NOT FORCED. Tokens that are not canonical point
   substitutions -- promoter / tandem-repeat calls (TR34), hotspot names, anything
   that is not wt-position-mut -- cannot map to a single residue, and say so. We do
   not pick a residue for them.

3. BUILDABILITY IS HONEST ABOUT THE ENVIRONMENT. The pinned env has no side-chain
   repacker, so only pure atom-DELETION mutations can be built deterministically
   (Y132F, V125A) -- the same whitelist mutate_receptor.py enforces. A mutation that
   ADDS atoms (K143R, F126L) is mapped (we know where it is and whether it touches
   the drug) but reported as needing a repacker to build the mutant side chain, not
   silently approximated. Mapping a residue and building its mutant side chain are
   two different claims; this module keeps them separate.

Stdlib only. Operates on `work/receptor_A.pdb` + `work/ref_VT1.pdb` and on the
substitution tokens emerging from openafr.emergence.detect_emergence.
"""
import math
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_RECEPTOR = os.path.join(ROOT, "work", "receptor_A.pdb")
DEFAULT_LIGAND = os.path.join(ROOT, "work", "ref_VT1.pdb")

# One-letter -> three-letter for all 20 standard amino acids, so ANY substitution
# token parses (not just the buildable whitelist) -- an unmappable mutation should
# still be recognised as a substitution and reported, not rejected as a typo.
AA1TO3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "E": "GLU",
    "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}

# Deterministic atom-DELETION mutations: the mutant side chain is a strict subset of
# the wild-type atoms in identical positions, so it is unique with no rotamer choice.
# This MUST stay in sync with scripts/mutate_receptor.py SUPPORTED -- test_mapping
# asserts they agree, so a change to one that is not mirrored in the other fails.
DELETION_MUTATIONS = {("TYR", "PHE"), ("VAL", "ALA")}

CONTACT_CUTOFF = 4.5   # A; heavy-atom contact distance defining a pocket-lining residue


def parse_substitution(token):
    """'Y132F' -> ('TYR', 132, 'PHE'), or None if `token` is not a point substitution.

    Case-insensitive. A stop codon (Y132*) parses with mut3 == '*'. Anything that is
    not <one-letter aa><integer><one-letter aa | *> returns None -- that is the
    signal to the caller that the token cannot be mapped to a single residue.
    """
    s = (token or "").strip().upper()
    m = re.fullmatch(r"([A-Z])(\d+)([A-Z*])", s)
    if not m:
        return None
    wt1, pos, mut1 = m.group(1), int(m.group(2)), m.group(3)
    if wt1 not in AA1TO3:
        return None
    mut3 = "*" if mut1 == "*" else AA1TO3.get(mut1)
    if mut3 is None:
        return None
    return AA1TO3[wt1], pos, mut3


def read_residues(receptor_pdb):
    """{position(int): (resname, [(atom_name, x, y, z), ...])} from an ATOM-record PDB.

    Heavy atoms only (hydrogens are skipped, matching openafr.pdbqt). Chain A only is
    assumed -- work/receptor_A.pdb is already a single-chain receptor.
    """
    residues = {}
    for l in open(receptor_pdb):
        if l[:6].strip() != "ATOM":
            continue
        name = l[12:16].strip()
        element = (l[76:78].strip() or name[0]).upper()
        if element == "H":
            continue
        pos = int(l[22:26])
        resname = l[17:20].strip()
        xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        entry = residues.setdefault(pos, (resname, []))
        entry[1].append((name, xyz))
    return residues


def read_ligand_atoms(ligand_pdb):
    """Heavy-atom coordinates [(x, y, z), ...] of the reference ligand."""
    atoms = []
    for l in open(ligand_pdb):
        if l[:6].strip() not in ("ATOM", "HETATM"):
            continue
        name = l[12:16].strip()
        element = (l[76:78].strip() or name[0]).upper()
        if element == "H":
            continue
        atoms.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return atoms


def _min_dist(residue_atoms, points):
    """Closest approach of any residue heavy atom to any point in `points` (or None)."""
    if not points:
        return None
    return min(math.dist(a[1], p) for a in residue_atoms for p in points)


def build_class(wt3, mut3):
    """'deletion' (deterministically buildable), or 'needs-repacker' (adds/repacks)."""
    return "deletion" if (wt3, mut3) in DELETION_MUTATIONS else "needs-repacker"


def map_mutation(token, residues, ligand_atoms, fe, contact_cutoff=CONTACT_CUTOFF):
    """Map one substitution token to its structural position in the modeled pocket.

    `residues`     -- read_residues(receptor)
    `ligand_atoms` -- read_ligand_atoms(reference ligand)
    `fe`           -- heme iron (x, y, z), from openafr.pdbqt.iron_position

    Returns a self-describing dict. `mappable` is False (with `reason`) when the token
    is not a point substitution, the residue is absent from the model, or the observed
    wild type does not match the token -- each an honest "cannot map", never a guess.
    A mappable verdict adds the structural readout: distance to the drug and to the
    iron, whether the residue lines the pocket, and whether the mutant can be built.
    """
    out = {"token": (token or "").strip().upper()}
    parsed = parse_substitution(token)
    if parsed is None:
        out.update(canonical=False, mappable=False,
                   reason="not a point substitution (promoter/tandem-repeat or "
                          "non-standard token) -- cannot map to a single residue")
        return out
    wt3, pos, mut3 = parsed
    out.update(canonical=True, wt=wt3, position=pos, mut=mut3,
               build_class=build_class(wt3, mut3))
    out["buildable"] = out["build_class"] == "deletion"

    if pos not in residues:
        out.update(mappable=False,
                   reason=f"residue {pos} is absent from the modeled structure "
                          f"(model covers {min(residues)}-{max(residues)})")
        return out
    observed, atoms = residues[pos]
    out["observed_wt"] = observed
    if observed != wt3:
        out.update(mappable=False, wt_matches=False,
                   reason=f"position {pos} is {observed} in the model, token expects "
                          f"{wt3} -- numbering or structure mismatch, refusing to map")
        return out

    out["wt_matches"] = True
    d_iron = _min_dist(atoms, [fe])
    d_lig = _min_dist(atoms, ligand_atoms)
    contact = d_lig is not None and d_lig <= contact_cutoff
    out.update(mappable=True, dist_to_iron=d_iron, dist_to_ligand=d_lig,
               contact_cutoff=contact_cutoff, pocket_contact=contact)
    out["structural_verdict"] = _verdict(out)
    return out


def _verdict(m):
    """One-line human summary of a mappable verdict."""
    where = (f"directly lines the azole pocket ({m['dist_to_ligand']:.2f} A from the "
             f"co-crystal drug)" if m["pocket_contact"]
             else f"is second-shell -- {m['dist_to_ligand']:.2f} A from the drug, "
                  f"does not contact it")
    build = ("mutant is deterministically buildable (atom deletion) -- can dock"
             if m["buildable"]
             else "mutant side chain needs a repacker to build -- mapped, not docked")
    return (f"{m['wt']}{m['position']}{m['mut']}: residue {where}, "
            f"{m['dist_to_iron']:.2f} A from the catalytic iron; {build}")


def map_signals(signals, receptor=DEFAULT_RECEPTOR, ligand=DEFAULT_LIGAND,
                contact_cutoff=CONTACT_CUTOFF):
    """Map a list of emergence signals (or bare tokens) to structural verdicts.

    Each element may be an emergence signal dict (uses its 'mutation' key) or a bare
    token string. Returns a list of map_mutation verdicts in the input order.
    """
    from openafr.pdbqt import iron_position
    residues = read_residues(receptor)
    ligand_atoms = read_ligand_atoms(ligand)
    fe = iron_position(receptor)
    tokens = [s["mutation"] if isinstance(s, dict) else s for s in signals]
    return [map_mutation(t, residues, ligand_atoms, fe, contact_cutoff)
            for t in tokens]
