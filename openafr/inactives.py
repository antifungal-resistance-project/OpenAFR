"""Verified-inactive decoy construction — measured non-inhibitors, not presumed ones.

Why this module exists
----------------------
The published ceiling (work/PREPRINT_geometry_ceiling.md §4.1) rests on decoys that are
*presumed* inactive: property-matched ChEMBL molecules that were never tested. That is the
paper's single largest caveat (Limitation 2), and §4.1 names the only clean exit —
"experimentally confirmed inactives — compounds measured not to inhibit". The §5 stopping
rule forbids another feature on the existing poses but explicitly allows "a new decoy
construction". This module builds that construction.

The one design rule: SINGLE-VARIABLE SWAP
-----------------------------------------
Everything that made the presumed-inactive set honest is held fixed — property matching
(MW/cLogP/rotatable bonds/HBD/HBA), the >= 1 aromatic nitrogen mechanism requirement, the
Tanimoto < 0.35 novelty floor, per-active pairing. The ONLY thing that changes is the
provenance of "inactive": measured, not assumed. If the matching changed too, a difference
in the gate could not be attributed to the measurement.

What counts as measured-inactive (the evidence rule)
----------------------------------------------------
Whole-cell MIC data against *C. albicans* from ChEMBL. Per record:

    inactive   MIC  >  X ug/mL with X >= 32   (tested to a real ceiling, nothing seen), or
               MIC  =  X ug/mL with X >= 64   (a measured, non-inhibitory MIC), or
               an explicit "Not Active"/"Inactive" activity comment
    active     MIC <= 8 ug/mL (or <= 10 uM), an "Active" comment, or pChEMBL >= 5
    ambiguous  measured, but between those bands (8-64 ug/mL), or censored below 32
    ignore     no usable quantitative or textual verdict

A compound is a **verified inactive** only if it has at least one `inactive` record and
ZERO `active` and ZERO `ambiguous` records — consistent inactivity across every assay it
appears in. One weak-but-real MIC anywhere disqualifies it. This is deliberately strict:
a false inactive in the decoy set depresses measured performance, and the whole point of
the exercise is to stop guessing about that.

The confound this cannot fully remove, stated plainly: a whole-cell MIC is not a CYP51
assay. A compound can inhibit the enzyme and still fail the cell assay (permeability,
efflux), which would make it a false inactive here. Where ChEMBL has direct CYP51 enzyme
data (target CHEMBL1780) we use it to disqualify such compounds; beyond that the residual
error is conservative — it depresses, not inflates, measured separation.
"""
import re

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

# Evidence thresholds (ug/mL unless noted). Frozen here, mirrored into the
# pre-registration; changing one is a construction change, not a tuning knob.
ACTIVE_UGML_MAX = 8.0        # <= this = active (the CLSI fluconazole-resistant breakpoint)
INACTIVE_UGML_MIN = 64.0     # an equality at or above this = measured non-inhibitory
CENSORED_MIN = 32.0          # "> X" only counts if the assay reached at least this far
ACTIVE_NM_MAX = 10000.0      # 10 uM, for records reported in nM
ACTIVE_PCHEMBL_MIN = 5.0     # ChEMBL's own potency flag

# Property matching + mechanism + novelty — IDENTICAL to scripts/make_decoys.py.
MW_TOL, LOGP_TOL, ROTB_TOL, HBD_TOL, HBA_TOL = 25.0, 1.5, 2, 1, 2
MAX_TANIMOTO = 0.35
MAX_HEAVY = 45               # the declared applicability domain (preprint §6.4)

_ACTIVE_COMMENT = re.compile(r"^\s*(active|potent)", re.I)
_INACTIVE_COMMENT = re.compile(r"(not\s+active|inactive|no\s+activity)", re.I)
_AROMATIC_N = Chem.MolFromSmarts("[n]")
_fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def classify_record(rec):
    """One ChEMBL activity record -> 'active' | 'inactive' | 'ambiguous' | 'ignore'.

    Textual verdicts are checked first because ChEMBL uses them exactly where it has no
    number to give. Quantitative rules only fire on units we can reason about.
    """
    comment = rec.get("activity_comment") or ""
    if _ACTIVE_COMMENT.search(comment):
        return "active"
    if _INACTIVE_COMMENT.search(comment):
        return "inactive"

    if _num(rec.get("pchembl_value")) is not None and _num(rec["pchembl_value"]) >= ACTIVE_PCHEMBL_MIN:
        return "active"

    val, rel = _num(rec.get("standard_value")), (rec.get("standard_relation") or "").strip()
    units = (rec.get("standard_units") or "").strip()
    if val is None:
        return "ignore"

    if units == "nM":
        # nM appears on enzyme/potency endpoints; only ever read as active evidence,
        # because a large nM value on a non-MIC endpoint is not a non-inhibitory MIC.
        return "active" if (rel in ("=", "<", "<=", "~") and val <= ACTIVE_NM_MAX) else "ignore"
    if units != "ug.mL-1":
        return "ignore"

    if rel in ("<", "<="):
        return "active" if val <= ACTIVE_UGML_MAX else "ambiguous"
    if rel in (">", ">="):
        return "inactive" if val >= CENSORED_MIN else "ambiguous"
    if rel in ("=", "~", ""):
        if val <= ACTIVE_UGML_MAX:
            return "active"
        if val >= INACTIVE_UGML_MIN:
            return "inactive"
        return "ambiguous"
    return "ignore"


def compound_verdict(records):
    """All of one compound's records -> (verdict, tally).

    'verified-inactive' requires consistent inactivity: >= 1 inactive record and no active
    or ambiguous ones. Any contrary evidence anywhere wins, because the cost of a false
    inactive in a benchmark is silent and asymmetric.
    """
    tally = {"active": 0, "inactive": 0, "ambiguous": 0, "ignore": 0}
    for r in records:
        tally[classify_record(r)] += 1
    if tally["active"]:
        return "has-active-evidence", tally
    if tally["ambiguous"]:
        return "ambiguous", tally
    if tally["inactive"]:
        return "verified-inactive", tally
    return "no-usable-evidence", tally


def props(mol):
    """The five matched properties + heavy-atom count (same set as make_decoys.py)."""
    return {
        "mw": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "rotb": Descriptors.NumRotatableBonds(mol),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "heavy": mol.GetNumHeavyAtoms(),
    }


def property_match(p, a):
    """Does candidate property dict p match active property dict a within the frozen tolerances?"""
    return (abs(p["mw"] - a["mw"]) <= MW_TOL
            and abs(p["logp"] - a["logp"]) <= LOGP_TOL
            and abs(p["rotb"] - a["rotb"]) <= ROTB_TOL
            and abs(p["hbd"] - a["hbd"]) <= HBD_TOL
            and abs(p["hba"] - a["hba"]) <= HBA_TOL)


def has_aromatic_nitrogen(mol):
    """The mechanism requirement: a plausible iron coordinator. Without it the criterion
    would only be asking 'does this molecule contain a nitrogen', which is circular."""
    return bool(mol.GetSubstructMatches(_AROMATIC_N))


def fingerprint(mol):
    return _fpgen.GetFingerprint(mol)


def load_actives(path):
    """Read a tab-separated .smi into [{name, mol, fp, ...props}]."""
    out = []
    for line in open(path):
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 2:
            continue
        smi, name = parts
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"unparseable active SMILES for {name}")
        out.append({"name": name, "mol": mol, "fp": fingerprint(mol), **props(mol)})
    return out


def select_inactives(candidates, actives, active_fps, per_active, max_heavy=MAX_HEAVY):
    """Pair verified-inactive compounds to actives under the frozen filters.

    candidates: iterable of (chembl_id, smiles), consumed in the order given (the caller
    fixes that order deterministically). Each compound is assigned to at most one active —
    the first, in the given active order, it property-matches — and each active takes at
    most `per_active`. Returns (chosen, rejects) where chosen is a list of dicts and
    rejects is a Counter-style dict of why-not counts.
    """
    quota = {a["name"]: 0 for a in actives}
    chosen, seen, rejects = [], set(), {
        "unparseable": 0, "duplicate": 0, "no-aromatic-n": 0,
        "over-domain": 0, "too-similar": 0, "no-property-match": 0, "quota-full": 0}
    for cid, smi in candidates:
        if smi in seen:
            rejects["duplicate"] += 1
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            rejects["unparseable"] += 1
            continue
        if not has_aromatic_nitrogen(mol):
            rejects["no-aromatic-n"] += 1
            continue
        p = props(mol)
        if p["heavy"] > max_heavy:
            rejects["over-domain"] += 1
            continue
        fp = fingerprint(mol)
        if max(DataStructs.BulkTanimotoSimilarity(fp, active_fps)) >= MAX_TANIMOTO:
            rejects["too-similar"] += 1
            continue
        partner = next((a for a in actives if property_match(p, a)), None)
        if partner is None:
            rejects["no-property-match"] += 1
            continue
        if quota[partner["name"]] >= per_active:
            rejects["quota-full"] += 1
            continue
        quota[partner["name"]] += 1
        seen.add(smi)
        chosen.append({"chembl_id": cid, "smiles": smi, "paired": partner["name"], **p})
    return chosen, rejects
