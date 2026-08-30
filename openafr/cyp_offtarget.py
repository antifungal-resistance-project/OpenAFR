"""Human drug-metabolizing-CYP off-target liability of a candidate (issue #74).

Why this module exists
----------------------
The whole pipeline selects molecules by ONE feature: an aromatic nitrogen that reaches
and coordinates the CYP51 heme iron. That is also, precisely, the mechanism by which the
marketed azoles poison the *human* cytochrome P450s that metabolize most drugs. Type-II
binding — an sp2 ring-nitrogen lone pair ligating the ferric heme iron as the sixth axial
ligand — does not know it is standing in a fungal pocket rather than CYP3A4. Ketoconazole
and itraconazole are textbook potent CYP3A4 inhibitors for exactly this reason, which is
why they carry black-box drug-interaction and hepatotoxicity warnings. So the feature our
geometry criterion REWARDS on the target is, unchanged, the feature a medicinal chemist
DREADS on the human panel. A geometry rank that never looks at the human CYPs is silent on
the single most predictable liability an azole-shaped candidate carries.

This module adds that counter-screen. For the five human CYPs that dominate xenobiotic
metabolism — 3A4, 2C9, 2D6, 2C19, 1A2 — it reports a per-isoform inhibition-liability
level from local, deterministic structural reasoning. It is a per-candidate CONFIDENCE
AXIS, not a gate: it annotates a rank, it never promotes and never hard-rejects.

What it reports (all coarse hints, none an IC50)
------------------------------------------------
Two, deliberately separated, signal types:

  1. TYPE-II HEME LIGATION (mechanistic, pan-CYP). An accessible azole warhead — an
     imidazole / 1,2,4-triazole / tetrazole / pyrazole ring carrying a pyridine-type
     (two-coordinate, unprotonated, no-H) ring nitrogen — is the antifungal pharmacophore
     itself, and it can axially ligate the heme iron of ANY cytochrome P450. An
     imidazole / triazole / tetrazole warhead (the strong, documented antifungal-class
     type-II ligand) is flagged HIGH against all five isoforms at once. A weaker ligator —
     a pyrazole (1,2-diazole) or a bare non-azole azine nitrogen (pyridine / pyrimidine /
     …) — is a real but softer type-II center and is flagged MODERATE pan-CYP.

  2. PER-ISOFORM SUBSTRATE PHARMACOPHORE (coarse, literature-grounded hints). Each human
     CYP has a well-known, recurring substrate shape. These are flagged MODERATE on the
     matching isoform only, and are explicitly hints, not predictions:
       3A4  — large + lipophilic (MW>=400, cLogP>=3): the promiscuous big-lipophilic sink.
       2D6  — a basic (protonatable) nitrogen near an aromatic ring: the classic 2D6 motif.
       2C9  — an anionic/acidic group (carboxylic acid / tetrazole) on a lipophilic body.
       2C19 — a moderately lipophilic aromatic amide/basic scaffold (the loosest hint).
       1A2  — a small, planar, poly-aromatic body (the flat-PAH-like substrate class).

Per isoform the liability is the STRONGER of the two signals: HIGH (azole type-II) >
MODERATE (azine type-II or a pharmacophore match) > LOW. The `basis` list records every
reason that fired, so a reader sees WHY, not just a level.

What this is NOT
----------------
Not a docking result and not a QSAR IC50. Docking these five large, famously flexible
pockets rigidly would manufacture false precision (the issue says as much); a trained
predictor would need a model file and a provenance we do not have offline. So this axis
stays honest about its resolution: it is a mechanistic + pharmacophore *liability flag*,
the same tier of claim as the hERG risk-factor hint in openafr.admet. For an approved-drug
candidate the authoritative source is the drug label (e.g. ketoconazole's is unambiguous),
not this estimate. The value here is that it fires on the NOVEL candidates too, before any
label exists, from the one feature they were selected for.

Network boundary
----------------
None. Everything is a local RDKit computation over a SMILES string; fully deterministic.
"""
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

# Reuse the SAME basic-nitrogen definition the ADMET axis uses, so "basic N" means one
# thing across the whole confidence stack.
from openafr.admet import _has_basic_nitrogen

RDLogger.DisableLog("rdApp.*")

# The human panel, in reporting order (metabolic dominance / clinical DDI weight).
PANEL = ("CYP3A4", "CYP2C9", "CYP2D6", "CYP2C19", "CYP1A2")

# --- per-isoform pharmacophore thresholds (coarse, frozen) ------------------------------
A34_MW = 400.0        # 3A4 favours large ...
A34_LOGP = 3.0        # ... lipophilic substrates
C9_LOGP = 2.0         # 2C9 favours lipophilic acids
C19_LOGP_LO = 2.0     # 2C19: a moderate-lipophilicity band
C19_LOGP_HI = 5.0
A2_MW = 300.0         # 1A2 favours small ...
A2_AROM_RINGS = 2     # ... planar poly-aromatic bodies
A2_AROM_FRACTION = 0.45  # heavy-atom fraction that is aromatic (flatness proxy)

# --- SMARTS (built once) ----------------------------------------------------------------
# A pyridine-type ring nitrogen: aromatic, exactly two ring connections, no H, not
# protonated -> an in-plane lone pair free to ligate the heme iron (the type-II center).
_LIGATING_N = Chem.MolFromSmarts("[nX2H0;!$([n+])]")
# Azole warheads: a 5-membered aromatic ring with >=2 ring nitrogens (imidazole,
# 1,2,4- and 1,2,3-triazole, tetrazole, pyrazole). The antifungal pharmacophore.
# Strong, documented antifungal-class type-II warheads -> HIGH pan-CYP.
_STRONG_WARHEADS = {
    "imidazole": Chem.MolFromSmarts("c1cnc[nX3,nX2]1"),
    "1,2,4-triazole": Chem.MolFromSmarts("c1n[nX2,nX3]cn1"),
    "1,2,3-triazole": Chem.MolFromSmarts("c1cn[nX2,nX3]n1"),
    "tetrazole": Chem.MolFromSmarts("c1nnn[nX2,nX3]1"),
}
# Weaker azole ligator (1,2-diazole) -> MODERATE pan-CYP.
_WEAK_WARHEADS = {
    "pyrazole": Chem.MolFromSmarts("c1cc[nX2,nX3][nX3,nX2]1"),
}
# 6-membered azine (pyridine / pyrimidine / pyrazine / triazine) carrying a ligating N.
_AZINE6 = Chem.MolFromSmarts("[nX2H0]1[c,n][c,n][c,n][c,n][c,n]1")
# An anionic / acidic warhead (2C9's substrate signature).
_CARBOXYLIC_ACID = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
_TETRAZOLE_ACID = Chem.MolFromSmarts("c1nnn[nH]1")
_ACYL_SULFONAMIDE = Chem.MolFromSmarts("[CX3](=O)[NX3][SX4](=O)(=O)")
# A (non-aromatic) amide, a common 2C19-substrate scaffold feature.
_AMIDE = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")


def _azole_warhead(mol):
    """(name, strength) of the azole warhead present, or (None, None).

    strength is "strong" (imidazole/triazole/tetrazole -> HIGH) or "weak" (pyrazole ->
    MODERATE). Requires the ring to also carry an accessible ligating nitrogen, so an
    azole whose every nitrogen is substituted/protonated (no lone pair free) is not counted.
    """
    if not mol.HasSubstructMatch(_LIGATING_N):
        return None, None
    for name, patt in _STRONG_WARHEADS.items():
        if mol.HasSubstructMatch(patt):
            return name, "strong"
    for name, patt in _WEAK_WARHEADS.items():
        if mol.HasSubstructMatch(patt):
            return name, "weak"
    return None, None


def _acidic_group(mol):
    return (mol.HasSubstructMatch(_CARBOXYLIC_ACID)
            or mol.HasSubstructMatch(_TETRAZOLE_ACID)
            or mol.HasSubstructMatch(_ACYL_SULFONAMIDE))


def _worst(a, b):
    """Return the stronger of two liability levels."""
    order = {"low": 0, "moderate": 1, "high": 2}
    return a if order[a] >= order[b] else b


def profile(smiles):
    """Full human-CYP off-target liability profile of one SMILES.

    Returns a dict; raises ValueError on an unparseable SMILES (callers decide skip/fail).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("unparseable SMILES: %r" % (smiles,))

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    heavy = mol.GetNumHeavyAtoms()
    arom_atoms = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
    arom_fraction = (arom_atoms / heavy) if heavy else 0.0
    basic_n = _has_basic_nitrogen(mol)
    acidic = _acidic_group(mol)

    warhead, strength = _azole_warhead(mol)
    # A bare (non-azole) azine nitrogen: has a ligating N but no azole ring -> weaker
    # type-II ligand. (Azole rings already imply a ligating N, so exclude them here.)
    azine_n = warhead is None and mol.HasSubstructMatch(_AZINE6)

    # --- assemble the per-isoform verdicts ---------------------------------------------
    panel = {cyp: {"liability": "low", "basis": []} for cyp in PANEL}

    def add(cyp, level, reason):
        panel[cyp]["liability"] = _worst(panel[cyp]["liability"], level)
        panel[cyp]["basis"].append(reason)

    # 1. mechanistic, pan-CYP type-II heme ligation
    if strength == "strong":
        for cyp in PANEL:
            add(cyp, "high", "type-II heme ligation (%s warhead)" % warhead)
    elif strength == "weak":
        for cyp in PANEL:
            add(cyp, "moderate", "type-II heme ligation (%s, weaker)" % warhead)
    elif azine_n:
        for cyp in PANEL:
            add(cyp, "moderate", "type-II heme ligation (azine N, weaker)")

    # 2. per-isoform substrate pharmacophore hints (coarse)
    if mw >= A34_MW and logp >= A34_LOGP:
        add("CYP3A4", "moderate", "large + lipophilic (MW>=%.0f, cLogP>=%.0f)" % (A34_MW, A34_LOGP))
    if basic_n and arom_rings >= 1:
        add("CYP2D6", "moderate", "basic N + aromatic (2D6 pharmacophore)")
    if acidic and logp >= C9_LOGP:
        add("CYP2C9", "moderate", "acidic group + lipophilic (2C9-preferred)")
    if (C19_LOGP_LO <= logp <= C19_LOGP_HI and arom_rings >= 1
            and (mol.HasSubstructMatch(_AMIDE) or basic_n)):
        add("CYP2C19", "moderate", "moderately lipophilic aromatic amide/basic (2C19 hint)")
    if (arom_rings >= A2_AROM_RINGS and mw <= A2_MW
            and arom_fraction >= A2_AROM_FRACTION):
        add("CYP1A2", "moderate", "small planar poly-aromatic (1A2-preferred)")

    n_high = sum(v["liability"] == "high" for v in panel.values())
    n_moderate = sum(v["liability"] == "moderate" for v in panel.values())
    flagged = [cyp for cyp in PANEL if panel[cyp]["liability"] != "low"]

    return {
        "smiles": smiles,
        "mw": round(mw, 1),
        "clogp": round(logp, 2),
        "azole_warhead": warhead,        # name or None (the antifungal pharmacophore)
        "warhead_strength": strength,    # "strong" (HIGH) / "weak" (MODERATE) / None
        "azine_nitrogen": azine_n,       # bare non-azole ligating azine N
        "type_ii_ligator": warhead is not None or azine_n,
        "panel": panel,                  # per-isoform {liability, basis}
        "n_high": n_high,
        "n_moderate": n_moderate,
        "flagged_isoforms": flagged,     # isoforms with liability != low
    }


def flag_summary(p):
    """One-glance flag string: what a reader should go check, worst-first.

    'type-II(imidazole):HIGHx5' == a strong azole warhead flags all five isoforms HIGH.
    A weaker ligator or a pure pharmacophore profile lists the matching isoforms. 'clean'
    == nothing fired.
    """
    if p["warhead_strength"] == "strong":
        return "type-II(%s):HIGHx%d" % (p["azole_warhead"], p["n_high"])
    if not p["flagged_isoforms"]:
        return "clean"
    short = [c.replace("CYP", "") for c in p["flagged_isoforms"]]
    if p["azole_warhead"]:
        lead = "%s:" % p["azole_warhead"]
    elif p["azine_nitrogen"]:
        lead = "azine:"
    else:
        lead = ""
    return lead + "+".join(short)
