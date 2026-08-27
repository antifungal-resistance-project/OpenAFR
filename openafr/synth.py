"""Synthesizability / availability — can a wet lab actually obtain this candidate?

Why this module exists
----------------------
A ranked hypothesis is only useful if someone can put the molecule in a tube. A geometry
rank says nothing about sourcing: a compound may be un-buyable and a nightmare to make.
This module adds the practicality axis the hand-off needs (#86): a synthetic-accessibility
score (local, deterministic) and an optional catalogue-presence probe.

Two signals, of very different strength
---------------------------------------
  SAscore   Ertl & Schuffenhauer's synthetic-accessibility score (2009), the RDKit-Contrib
            implementation. A 1 (trivial) .. 10 (very hard) heuristic from fragment
            contributions + a complexity penalty (ring fusion, stereo, size, macrocycles).
            Fully LOCAL and deterministic — no network. This is the substantive signal.
  presence  An OPTIONAL, best-effort PubChem lookup: does a compound record even exist for
            this name/SMILES? A CID is a *weak* proxy for "buyable" — presence in PubChem is
            not a vendor quote — so it is off by default and clearly labelled. It lives
            behind an injectable `opener` (default: urllib), exactly like openafr.precedent,
            so the parsing is tested offline with canned JSON and a live call can never take
            down a run.

Verdict bands (SAscore, informational)
--------------------------------------
  easy        SAscore <  SA_EASY  (<= ~3.5): routine, most drug-like molecules land here
  moderate    SA_EASY .. SA_HARD : synthesizable with effort
  hard        SAscore >= SA_HARD  (>= ~6): a real synthetic project — sourcing risk

These are soft annotations on a rank, never a gate. The current shortlist is repurposed
*approved/clinical* drugs, which are by definition buyable — the check matters most for any
future non-drug library, and is wired now so it is ready then.

What this cannot tell you (stated plainly)
------------------------------------------
SAscore is a heuristic, not a route; a low score is not a guarantee and a high one is not a
veto (nature makes hard molecules). PubChem presence is not a purchase order — a true
buyable check needs a vendor catalogue API (ZINC/eMolecules), out of scope here. For an
approved drug the honest answer is "buyable, see the drug's catalogue", which this module
records as `approved_drug` when told so, rather than pretending to price it.
"""
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request

from rdkit import Chem, RDLogger
from rdkit.Chem import RDConfig

RDLogger.DisableLog("rdApp.*")

# SAscore bands (Ertl scale 1..10)
SA_EASY = 3.5
SA_HARD = 6.0

# --- SAscore: import the RDKit-Contrib implementation once -------------------------------
_SASCORER = None


def _sascorer():
    """Lazily import RDKit-Contrib SA_Score/sascorer (loads fpscores next to itself)."""
    global _SASCORER
    if _SASCORER is None:
        sa_dir = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if sa_dir not in sys.path:
            sys.path.append(sa_dir)
        import sascorer  # noqa: E402  (path-appended contrib module)
        _SASCORER = sascorer
    return _SASCORER


def sa_score(mol):
    """Ertl SAscore (1 easy .. 10 hard) for an RDKit mol."""
    return _sascorer().calculateScore(mol)


def sa_band(score):
    if score < SA_EASY:
        return "easy"
    if score < SA_HARD:
        return "moderate"
    return "hard"


# --- optional PubChem presence probe (weak "does a record exist" signal) -----------------
_PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def _get_json(url, opener):
    try:
        with opener(url, timeout=20) as r:
            return json.loads(r.read().decode("utf-8"))
    except (urllib.error.URLError, urllib.error.HTTPError, ValueError, TimeoutError):
        return None


def pubchem_cid(name, opener=urllib.request.urlopen):
    """First PubChem CID for a compound name, or None. Best-effort, never raises.

    A found CID is a WEAK proxy for availability — presence in PubChem is not a vendor
    quote. Behind an injectable opener so parsing is tested offline.
    """
    url = "%s/compound/name/%s/cids/JSON" % (_PUBCHEM, urllib.parse.quote(name))
    data = _get_json(url, opener)
    try:
        cids = data["IdentifierList"]["CID"]
        return cids[0] if cids else None
    except (TypeError, KeyError, IndexError):
        return None


def profile(smiles, name=None, catalogued=False, opener=None):
    """Synthesizability profile for one SMILES.

    - SAscore + band are always computed (local, deterministic).
    - availability: 'catalogued-drug' when catalogued=True (the molecule comes from a
      curated buyable drug library, e.g. the Broad Drug Repurposing Hub — approved,
      clinical or catalogued tool compounds, all vendor-listed); else, if an opener is
      given, 'pubchem-record' (weak presence proxy) or 'no-pubchem-record'; else
      'not-checked'.
    Raises ValueError on an unparseable SMILES (the pipeline refuses rather than guesses).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("unparseable SMILES: %r" % (smiles,))
    score = round(sa_score(mol), 2)

    if catalogued:
        availability, cid = "catalogued-drug", None
    elif opener is not None and name:
        cid = pubchem_cid(name, opener=opener)
        availability = "pubchem-record" if cid else "no-pubchem-record"
    else:
        availability, cid = "not-checked", None

    return {
        "name": name,
        "smiles": smiles,
        "sa_score": score,
        "sa_band": sa_band(score),
        "availability": availability,
        "pubchem_cid": cid,
    }
