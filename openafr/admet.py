"""ADMET / physicochemical + structural-alert profiling of a candidate.

Why this module exists
----------------------
A molecule can coordinate the heme iron beautifully (the geometry criterion the whole
pipeline is built on) and still be un-druggable: too greasy to dissolve, too big to
permeate, carrying a reactive warhead or a known toxicophore. The geometry rank says
nothing about any of this. This module adds the orthogonal chemistry-realism axis the
shortlist was missing (#75): before a molecule is proposed for the bench, profile its
drug-likeness and flag the substructures medicinal chemists have learned to distrust.

What it reports (all informational, none a hard gate)
-----------------------------------------------------
Physicochemical profile, per molecule:
  MW, cLogP (Crippen), H-bond donors/acceptors, TPSA, rotatable bonds, aromatic-ring
  count, and an ESOL aqueous-solubility estimate (Delaney's linear model). These are the
  inputs to two classic, deliberately-simple drug-likeness rulesets:
    Lipinski Ro5   MW<=500, cLogP<=5, HBD<=5, HBA<=10  -> count violations
    Veber          rotatable bonds<=10 AND TPSA<=140   -> oral-bioavailability proxy
  We report violation *counts*, never a pass/fail verdict: the issue is explicit that
  these are informational, not gates, and every marketed azole below already bends a rule
  (ketoconazole is a well-known Ro5 offender and still a drug).

Structural alerts (RDKit FilterCatalog, offline, no network):
  PAINS (A/B/C)  pan-assay-interference substructures — frequent false positives in
                 screens; a hit here says "this may be reading the assay, not the target".
  BRENK          unwanted/reactive/putative-tox fragments (Brenk et al. 2008).
  NIH            NIH/MLSMR reactive-and-unwanted filters.
  Each catalog hit carries the matched fragment's description so it can be read, not just
  counted.

hERG cardiotoxicity (coarse, honestly labelled):
  We do NOT ship a predictive hERG model. We flag the classic hERG pharmacophore *risk
  factors* only — a protonatable (basic) nitrogen together with lipophilicity (cLogP high)
  and aromatic bulk — the recurring shape of hERG blockers (Aronov 2005). This is a
  "go read the cardiac liability" hint, explicitly not a prediction. See HERG_* below.

What this cannot tell you (stated plainly)
------------------------------------------
Every number here is a computed estimate, not a measurement. Ro5/Veber are 30-year-old
heuristics with famous exceptions; PAINS/BRENK flag *possible* liabilities that context
can excuse; the hERG flag is a pharmacophore hint, not an IC50. For candidates that are
already approved drugs, the authoritative source is the drug label, not these estimates —
this module annotates a rank, it never overrules the wet lab or the clinic.

Network boundary
----------------
None. Everything here is a local RDKit computation over a SMILES string; the module makes
no network calls and is fully deterministic.
"""
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

RDLogger.DisableLog("rdApp.*")

# --- drug-likeness thresholds (classic, frozen) -----------------------------------------
RO5_MW = 500.0
RO5_LOGP = 5.0
RO5_HBD = 5
RO5_HBA = 10
VEBER_ROTB = 10
VEBER_TPSA = 140.0

# hERG pharmacophore risk-factor thresholds (coarse hint, NOT a predictive model)
HERG_LOGP = 3.7        # lipophilicity above which hERG risk rises (Aronov 2005, roughly)
HERG_AROM_RINGS = 2    # aromatic bulk


def _esol_logs(mol, logp):
    """Delaney (2004) ESOL estimate of aqueous log-solubility (mol/L).

    logS = 0.16 - 0.63*cLogP - 0.0062*MW + 0.066*RotB - 0.74*AP
    where AP is the aromatic proportion (aromatic heavy atoms / heavy atoms).
    """
    mw = Descriptors.MolWt(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    heavy = mol.GetNumHeavyAtoms()
    arom = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
    ap = (arom / heavy) if heavy else 0.0
    return 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotb - 0.74 * ap


def _has_basic_nitrogen(mol):
    """A protonatable (basic) nitrogen: an sp3 amine or an amidine/guanidine-like N.

    Deliberately coarse — excludes amide, nitro, aromatic-pyridine-type and nitrile N,
    which are not appreciably basic at physiological pH. Matches the hERG-pharmacophore
    "basic center" and the general medchem notion of a protonatable site.
    """
    # aliphatic primary/secondary/tertiary amine (not amide, not adjacent to C=O)
    amine = Chem.MolFromSmarts("[NX3;!$(NC=[O,S,N]);!$(N=*);!$([N+])]")
    amidine = Chem.MolFromSmarts("[NX3][CX3]=[NX2]")  # amidine / guanidine
    return mol.HasSubstructMatch(amine) or mol.HasSubstructMatch(amidine)


def _alert_catalog():
    """FilterCatalog loaded with PAINS (A/B/C), BRENK, and NIH rule sets."""
    params = FilterCatalogParams()
    for cat in (
        FilterCatalogParams.FilterCatalogs.PAINS_A,
        FilterCatalogParams.FilterCatalogs.PAINS_B,
        FilterCatalogParams.FilterCatalogs.PAINS_C,
        FilterCatalogParams.FilterCatalogs.BRENK,
        FilterCatalogParams.FilterCatalogs.NIH,
    ):
        params.AddCatalog(cat)
    return FilterCatalog(params)


# module-level singleton — building the catalog is not free
_CATALOG = None


def _catalog():
    global _CATALOG
    if _CATALOG is None:
        _CATALOG = _alert_catalog()
    return _CATALOG


def structural_alerts(mol):
    """List of (catalog, description) for every PAINS/BRENK/NIH match. Empty == clean."""
    hits = []
    for entry in _catalog().GetMatches(mol):
        # e.g. "PAINS_A", "Brenk", "NIH" -> the family lives in the FilterSet prop
        props = entry.GetPropList()
        family = entry.GetProp("FilterSet") if "FilterSet" in props else ""
        hits.append((family, entry.GetDescription()))
    return hits


def profile(smiles):
    """Full ADMET/physicochemical profile of one SMILES.

    Returns a dict; raises ValueError on an unparseable SMILES so a caller can decide
    whether to skip or fail (the pipeline elsewhere refuses rather than guesses).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("unparseable SMILES: %r" % (smiles,))

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    heavy = mol.GetNumHeavyAtoms()
    logs = _esol_logs(mol, logp)

    ro5_violations = sum(
        (mw > RO5_MW, logp > RO5_LOGP, hbd > RO5_HBD, hba > RO5_HBA)
    )
    veber_ok = rotb <= VEBER_ROTB and tpsa <= VEBER_TPSA

    alerts = structural_alerts(mol)
    pains = [d for fam, d in alerts if fam.startswith("PAINS")]
    other_alerts = [(fam, d) for fam, d in alerts if not fam.startswith("PAINS")]

    basic_n = _has_basic_nitrogen(mol)
    herg_risk = basic_n and logp >= HERG_LOGP and arom_rings >= HERG_AROM_RINGS

    return {
        "smiles": smiles,
        "mw": round(mw, 1),
        "clogp": round(logp, 2),
        "hbd": hbd,
        "hba": hba,
        "tpsa": round(tpsa, 1),
        "rotatable_bonds": rotb,
        "aromatic_rings": arom_rings,
        "heavy_atoms": heavy,
        "esol_logs": round(logs, 2),      # log10 solubility (mol/L); more negative = less soluble
        "ro5_violations": ro5_violations,
        "veber_ok": veber_ok,
        "pains": pains,                   # list of PAINS descriptions (empty == clean)
        "structural_alerts": other_alerts,  # (family, description) BRENK/NIH
        "n_alerts": len(alerts),
        "basic_nitrogen": basic_n,
        "herg_risk_factors": herg_risk,   # coarse pharmacophore hint, NOT a prediction
    }
