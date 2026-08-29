"""Assay precedent — has this candidate already been measured against the target?

Why this module exists
----------------------
The triage pipeline ranks *novel* molecules by heme-iron geometry (rank_candidates.py),
but "novel to us" is not "never tested". A wet lab may already have assayed a molecule
against fungal CYP51 and found it inactive — and, because of publication bias, that result
often never reaches a paper (it sits in a file drawer). Ranking such a compound highly and
proposing it for the bench re-runs an experiment someone already lost. This module asks the
one question the ranker cannot: *is there a deposited assay record for this molecule against
the target, and what did it say?*

Three sources, because inactives (and known actives) hide in all of them
------------------------------------------------------------------------
  ChEMBL    curated literature + HTS activities, quantitative (MIC/IC50). Verdict comes from
            the measured value via openafr.inactives (CLSI-anchored bands).
  PubChem   BioAssay — the largest store of HTS *inactives*, each assay pre-adjudicated by the
            depositor as Active / Inactive / Inconclusive / Unspecified. Verdict comes from
            that outcome directly; 'Unspecified' is read as no verdict (ignored), never as a
            non-inhibition claim.
  BindingDB curated binding affinities (IC50/Ki/Kd, nM) with a literature citation (PMID/DOI)
            per record — a store of measured *actives* that the other two often miss, because
            BindingDB harvests medicinal-chemistry SAR tables ChEMBL may not have curated (#78).
            Queried target-first (every ligand ever measured against the CYP51 enzyme) and then
            matched to a candidate locally by InChIKey, the inverse of the compound-first ChEMBL/
            PubChem probes. Its nM affinities are read through the SAME nM rule the ChEMBL path
            uses — only ever as *active* evidence, never as a non-inhibitory verdict.
The sources are merged at the tally level, so a candidate flagged inactive in any is caught,
and one active record in any flips the verdict — the same 'contrary evidence wins' rule
inactives.verdict_from_tally encodes, applied once across all sources.

What it can and cannot tell you (stated plainly)
------------------------------------------------
A hit here is strong: the molecule was measured and the record exists. But a *miss* is weak:

    NO record  !=  never tested.

It may mean only that a failing measurement was never deposited. This module reports
"no-precedent", never "untested" — the absence is not evidence.

On-target vs. phenotypic vs. off-target
---------------------------------------
Records are sorted into three buckets, and the sort differs by source because the two label
their targets differently:
  on_target   the CYP51/ERG11 enzyme itself — the only direct verdict. ChEMBL: target id in
              CYP51_TARGETS. PubChem: assay name names the enzyme (cyp51/erg11/sterol
              14-demethylase), since PubChem has no single CYP51 target id to filter on.
  phenotypic  a whole-cell MIC against a FUNGUS — ChEMBL: a non-enzyme µg/mL record whose
              organism (target_organism/assay_organism) is a fungus; PubChem: an assay named
              antifungal/Candida/Aspergillus/etc. Weaker (permeability/efflux confound) and a
              coarse hint to go read the assays, never a fungal verdict on its own. It is now
              organism-filtered on the ChEMBL side (#79): an antibacterial or mammalian-cell µg/mL
              record is NOT antifungal evidence and drops to off_target, so a `tested-active`
              here is a fungal cell result, not (e.g.) an antibacterial one.
  off_target  everything else — not a resistance signal, a *novelty* signal: a compound already
              characterised against other targets is a known molecule, not a fresh chemotype.

Near-neighbour precedent
------------------------
Even when the exact molecule is absent, a very close analogue measured inactive is a soft
negative. nearest_tested_neighbour() reports the most similar previously-tested compound
(Tanimoto over the same Morgan fingerprint the pipeline uses elsewhere). It annotates a
geometry rank, never overrides it.

Network boundary
----------------
Every ChEMBL/PubChem/BindingDB call lives behind an injectable `opener` (default: urllib) so the
pure bucket/verdict logic is exercised offline with canned JSON. The network functions make no
promises about availability and fail per-molecule, never taking down a whole run. BindingDB's
target-first index is built once (one call per CYP51 UniProt) and injected into lookup(), so the
per-candidate path stays a local InChIKey match with no extra network per molecule.
"""
import json
import re
import urllib.error
import urllib.parse
import urllib.request

from rdkit import Chem, DataStructs, RDLogger

from openafr import inactives

RDLogger.DisableLog("rdApp.*")

# ChEMBL target(s) that ARE the enzyme we triage against. C. albicans sterol 14-alpha
# demethylase (CYP51 / ERG11); the same target inactives.py names for disqualifying
# enzyme-level actives. Frozen here so a precedent verdict is auditable; extend only with a
# target that is genuinely the same enzyme (e.g. an orthologue) and say so in the run notes.
CYP51_TARGETS = frozenset({"CHEMBL1780"})

# BindingDB UniProt accession(s) that ARE the enzyme we triage against — the same target the
# ChEMBL CYP51_TARGETS id names and the docked 5TZ1 structure is: C. albicans sterol 14-alpha
# demethylase (ERG11 / CYP51), UniProt P10613. Frozen for auditability; extend only with a
# genuine CYP51 orthologue and say so in the run notes. BindingDB is queried target-first, so a
# candidate's OFF-target BindingDB records are deliberately not collected — only the direct,
# on-target enzyme verdict this source uniquely adds beyond ChEMBL + PubChem (#78).
CYP51_UNIPROTS = frozenset({"P10613"})

# PubChem has no single CYP51 target id, so its buckets are read from the assay name. Kept
# deliberately specific: ENZYME must name the sterol-demethylase itself (a bare "demethylase"
# would sweep in unrelated histone/DNA demethylases), PHENOTYPIC names a whole-organism assay.
_ENZYME_ASSAY = re.compile(
    r"cyp51|erg11|sterol\s*14|14\W*alpha[\s\-]*demethyl|lanosterol[\s\-]*14", re.I)
_PHENOTYPIC_ASSAY = re.compile(
    r"antifungal|candida|aspergillus|cryptococc|trichophyton|dermatophyt|\bfungal\b", re.I)

# A whole-cell µg/mL MIC bears on ANTIFUNGAL activity only if it was run against a fungus. The
# ChEMBL side of group_by_target used to route EVERY non-enzyme µg/mL record into the phenotypic
# bucket regardless of organism, so an antibacterial MIC (flucloxacillin/tazobactam vs
# Staphylococcus/E. coli) or a mammalian-cell-line assay counted as an antifungal signal — the
# exact noise issue #79 removes. This matches an organism string (ChEMBL target_organism /
# assay_organism) against the clinically-relevant fungal genera plus the generic terms, and is
# kept aligned in spirit with _PHENOTYPIC_ASSAY (the PubChem side is already fungal-name-gated,
# so PubChem phenotypic records were never contaminated — only the ChEMBL side needed this).
_FUNGAL_ORGANISM = re.compile(
    r"\bfung|\bmould\b|\bmold\b|\byeast|candida|nakaseomyces|clavispora|meyerozyma|"
    r"pichia|debaryomyces|kluyveromyces|saccharomyces|schizosaccharomyces|zygosaccharomyces|"
    r"aspergillus|penicillium|talaromyces|paecilomyces|purpureocillium|scopulariopsis|"
    r"acremonium|fusarium|scedosporium|lomentospora|cryptococc|trichosporon|rhodotorula|"
    r"malassezia|trichophyton|microsporum|nannizzia|epidermophyton|dermatophyt|"
    r"histoplasma|blastomyces|coccidioides|paracoccidioides|sporothrix|"
    r"rhizopus|mucor|rhizomucor|lichtheimia|cunninghamella|absidia|apophysomyces|"
    r"exophiala|cladosporium|alternaria|curvularia|geotrichum|magnaporthe|"
    r"pneumocystis|phialophora|fonsecaea", re.I)

# PubChem's own adjudication -> our shared vocabulary. 'Unspecified' is NOT a non-inhibition
# claim (the depositor simply did not call it), so it is ignored, never read as inactive.
_PUBCHEM_LABEL = {"Active": "active", "Inactive": "inactive", "Inconclusive": "ambiguous"}

# A near neighbour this similar to a previously-tested compound inherits that compound's
# result as a *soft* signal. Set well above the novelty floor make_decoys.py uses (0.35): we
# only speak up when the analogy is strong enough that a chemist would not dismiss it.
NEIGHBOUR_TANIMOTO = 0.85

_CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
_PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_BINDINGDB_REST = "https://bindingdb.org/rest"

# candidate-facing label <- inactives verdict label. The two are the same finding: a molecule
# that would qualify as a verified-inactive decoy is one a wet lab measured inactive.
_VERDICT_LABEL = {
    "has-active-evidence": "tested-active",
    "verified-inactive": "tested-inactive",
    "ambiguous": "tested-ambiguous",
    "no-usable-evidence": "record-inconclusive",
}


def inchikey(smiles):
    """Canonical InChIKey for a SMILES, or None if RDKit cannot parse/standardise it.

    The InChIKey is the join key against both databases: it collapses tautomer/charge notation
    that differs between our library and a deposited record, so two writings of the same
    molecule still match. A None here is a candidate we cannot look up, reported as such —
    never a silent skip.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    key = Chem.MolToInchiKey(mol)
    return key or None


# --- the source-agnostic bucket/tally core --------------------------------------------

def _empty_buckets():
    """Per-bucket tallies of classified labels. Merging two sources is a field-wise sum."""
    return {b: {"active": 0, "inactive": 0, "ambiguous": 0, "ignore": 0}
            for b in ("on_target", "phenotypic", "off_target")}


def _add(tally, label):
    tally[label] += 1


def merge_buckets(a, b):
    """Field-wise sum of two bucket dicts — how ChEMBL and PubChem evidence is combined."""
    out = _empty_buckets()
    for bucket in out:
        for label in out[bucket]:
            out[bucket][label] = a[bucket][label] + b[bucket][label]
    return out


def _bucket_verdict(tally):
    """One bucket's tally -> {'verdict', 'tally'}, or None when the bucket is empty."""
    if sum(tally.values()) == 0:
        return None
    return {"verdict": _VERDICT_LABEL[inactives.verdict_from_tally(tally)], "tally": tally}


def render(buckets):
    """Bucketed tallies -> the candidate-facing summary dict.

    Keys: verdict (on-target, candidate-facing; 'no-precedent' when the enzyme bucket is
    empty), tally, phenotypic (secondary whole-cell verdict or None), n_off_target (novelty
    signal), n_records. This is the pure renderer both sources feed; it touches no network.
    """
    on = _bucket_verdict(buckets["on_target"])
    pheno = _bucket_verdict(buckets["phenotypic"])
    n_off = sum(buckets["off_target"].values())
    n_records = sum(sum(t.values()) for t in buckets.values())
    return {
        "verdict": on["verdict"] if on else "no-precedent",
        "tally": on["tally"] if on else None,
        "phenotypic": pheno,
        "n_off_target": n_off,
        "n_records": n_records,
    }


def fingerprint(smiles):
    """Morgan fingerprint for a SMILES (the pipeline's shared one), or None if unparseable."""
    mol = Chem.MolFromSmiles(smiles)
    return inactives.fingerprint(mol) if mol is not None else None


def load_neighbours(path):
    """Read a previously-tested set for the near-neighbour check.

    Format per line: `smiles<TAB>id[<TAB>verdict]`. The verified-inactives set this repo already
    builds (scripts/make_verified_inactives.py) is the intended input — every entry there is a
    measured inactive, so an omitted third column defaults to 'tested-inactive'. Unparseable
    SMILES are skipped (they cannot seed a similarity), never silently miscounted.
    """
    out = []
    for line in open(path):
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        smi, ident = parts[0], parts[1]
        verdict = parts[2] if len(parts) > 2 else "tested-inactive"
        fp = fingerprint(smi)
        if fp is not None:
            out.append({"id": ident, "fp": fp, "verdict": verdict})
    return out


def nearest_tested_neighbour(query_fp, neighbours):
    """Most similar previously-tested compound to a candidate, above NEIGHBOUR_TANIMOTO.

    query_fp   : the candidate's Morgan fingerprint (inactives.fingerprint).
    neighbours : [{"id", "fp", "verdict"}, ...] — compounds with a known precedent verdict.
    Returns the single best {"id", "verdict", "similarity"} at or above the threshold, or None
    when nothing is close enough. A soft signal only: it annotates a rank, never changes it.
    """
    best = None
    for n in neighbours:
        s = DataStructs.TanimotoSimilarity(query_fp, n["fp"])
        if s >= NEIGHBOUR_TANIMOTO and (best is None or s > best["similarity"]):
            best = {"id": n["id"], "verdict": n["verdict"], "similarity": s}
    return best


# --- ChEMBL adapter -------------------------------------------------------------------

def record_organism(rec):
    """The organism a ChEMBL activity record was measured against ('' if none is recorded).

    target_organism is the assay's biological subject — for a whole-cell MIC, the organism whose
    cells were dosed; assay_organism is the fallback. Either is what tells a µg/mL record apart as
    fungal (relevant to antifungal activity) or not. Only these two structured fields are read;
    the free-text assay_description is deliberately NOT mined, to avoid a stray "vs Candida
    comparator" mention flipping a non-fungal assay to fungal.
    """
    return (rec.get("target_organism") or rec.get("assay_organism") or "").strip()


def is_fungal_organism(name):
    """True iff `name` denotes a fungus — so a whole-cell MIC against it is antifungal evidence.

    An empty/unrecognised organism is False: a µg/mL record we cannot confirm is fungal must not
    masquerade as an antifungal signal (it becomes an off-target novelty count instead, #79).
    """
    return bool(name) and _FUNGAL_ORGANISM.search(name) is not None


def group_by_target(records):
    """Split ChEMBL activity records into (on_target, phenotypic, off_target).

    on_target  : records against a CYP51_TARGETS enzyme target — the direct verdict.
    phenotypic : a whole-cell µg/mL MIC run against a **fungus** — the only cell assay that bears
                 on antifungal activity. Still weaker than the enzyme verdict (permeability/efflux
                 confound), reported separately, never an enzyme verdict on its own — but now
                 organism-filtered (#79), so an antibacterial or mammalian-cell MIC no longer
                 pollutes it.
    off_target : everything else — a 'this is a known compound' novelty signal, kept as a count.
                 This now ALSO catches a non-enzyme µg/mL record whose organism is a non-fungus
                 OR is unrecorded: it is still a measured record (a novelty signal) but makes no
                 antifungal claim. A record with no target id likewise cannot support any claim.
    """
    on_target, phenotypic, off_target = [], [], []
    for r in records:
        tid = r.get("target_chembl_id")
        if tid in CYP51_TARGETS:
            on_target.append(r)
        elif (r.get("standard_units") or "").strip() == "ug.mL-1" \
                and is_fungal_organism(record_organism(r)):
            phenotypic.append(r)
        else:
            off_target.append(r)
    return on_target, phenotypic, off_target


def phenotypic_organism_breakdown(records):
    """From a molecule's ChEMBL records, the organisms behind its whole-cell µg/mL MICs.

    Returns {'fungal': [...], 'nonfungal': [...]} — sorted, distinct organism names — so a reader
    can see *which* organisms a phenotypic verdict rests on, and which non-fungal MICs were
    reclassified to off_target rather than counted as antifungal evidence (#79). An unrecorded
    organism is listed under 'nonfungal' as '(unspecified)', since it too was excluded from the
    fungal verdict. PubChem records carry no structured organism field and never appear here.
    """
    fungal, nonfungal = set(), set()
    for r in records:
        if r.get("target_chembl_id") in CYP51_TARGETS:
            continue
        if (r.get("standard_units") or "").strip() != "ug.mL-1":
            continue
        org = record_organism(r)
        if is_fungal_organism(org):
            fungal.add(org)
        else:
            nonfungal.add(org or "(unspecified)")
    return {"fungal": sorted(fungal), "nonfungal": sorted(nonfungal)}


def chembl_buckets(records):
    """ChEMBL activity records -> bucketed tallies, via target grouping + inactives.classify_record."""
    on, pheno, off = group_by_target(records)
    b = _empty_buckets()
    for bucket, recs in (("on_target", on), ("phenotypic", pheno), ("off_target", off)):
        for r in recs:
            _add(b[bucket], inactives.classify_record(r))
    return b


def candidate_verdict(on_target_records):
    """On-target ChEMBL records -> a candidate-facing precedent verdict {'verdict', 'tally'}.

    With no records the verdict is 'no-precedent' — deliberately NOT 'untested', because a
    missing record does not prove the molecule was never measured (see the module docstring).
    """
    if not on_target_records:
        return {"verdict": "no-precedent", "tally": None}
    raw, tally = inactives.compound_verdict(on_target_records)
    return {"verdict": _VERDICT_LABEL[raw], "tally": tally}


def summarise(records):
    """All of a molecule's ChEMBL activity records -> the full precedent summary (ChEMBL only).

    Adds `phenotypic_organisms` (the fungal/non-fungal organism breakdown, #79) to render()'s
    keys, so the whole-cell verdict carries the organisms it rests on.
    """
    summary = render(chembl_buckets(records))
    summary["phenotypic_organisms"] = phenotypic_organism_breakdown(records)
    return summary


# --- PubChem adapter ------------------------------------------------------------------

def classify_pubchem_outcome(outcome):
    """PubChem 'Activity Outcome' -> our shared label. 'Unspecified'/unknown -> 'ignore'."""
    return _PUBCHEM_LABEL.get(outcome, "ignore")


def pubchem_bucket(assay_name):
    """Route a PubChem assay to a bucket by its name (PubChem has no CYP51 target id)."""
    name = assay_name or ""
    if _ENZYME_ASSAY.search(name):
        return "on_target"
    if _PHENOTYPIC_ASSAY.search(name):
        return "phenotypic"
    return "off_target"


def pubchem_buckets(rows):
    """PubChem assay rows [{'outcome', 'assay_name'}, ...] -> bucketed tallies."""
    b = _empty_buckets()
    for r in rows:
        _add(b[pubchem_bucket(r.get("assay_name"))], classify_pubchem_outcome(r.get("outcome")))
    return b


# --- BindingDB adapter ----------------------------------------------------------------

# Inequality prefixes a BindingDB affinity may carry; longest-first so '<=' is matched before '<'.
_BDB_RELATIONS = ("<=", ">=", "<", ">", "=", "~")


def parse_bindingdb_affinity(raw):
    """A BindingDB affinity string ('250', '>10000', '<0.5') -> (relation, value_nM), or None.

    BindingDB reports IC50/Ki/Kd in nM, optionally prefixed with an inequality. A value we cannot
    read as a number is dropped (None), never guessed — one unreadable affinity must not fabricate
    a verdict. A bare number carries the '=' relation.
    """
    s = (raw or "").strip()
    rel = "="
    for r in _BDB_RELATIONS:
        if s.startswith(r):
            rel, s = r, s[len(r):].strip()
            break
    try:
        return rel, float(s)
    except ValueError:
        return None


def bindingdb_record_to_activity(rec):
    """A normalised BindingDB record -> a ChEMBL-shaped activity dict for inactives.classify_record.

    The affinity is an nM enzyme potency (IC50/Ki/Kd), so it flows through the SAME nM rule the
    ChEMBL path uses: only ever read as *active* evidence, because a large nM value on a binding
    endpoint is not a non-inhibitory MIC (openafr.inactives). Returns None when the affinity is
    unparseable — the record then supports no verdict rather than a guessed one.
    """
    parsed = parse_bindingdb_affinity(rec.get("affinity"))
    if parsed is None:
        return None
    rel, val = parsed
    return {"standard_value": val, "standard_relation": rel, "standard_units": "nM",
            "activity_comment": None, "pchembl_value": None}


def bindingdb_buckets(records):
    """Normalised BindingDB records -> bucketed tallies (all on-target by construction).

    Every record here was fetched target-first against a CYP51_UNIPROTS enzyme, so it lands in the
    on_target bucket — BindingDB adds a direct enzyme verdict, never a phenotypic or off-target
    count. An affinity that will not parse is ignored, not tallied.
    """
    b = _empty_buckets()
    for r in records:
        act = bindingdb_record_to_activity(r)
        if act is not None:
            _add(b["on_target"], inactives.classify_record(act))
    return b


# --- network layer (thin, injectable, not exercised against the real services) --------

def _get_json(url, opener):
    with opener(url, timeout=30) as resp:
        return json.load(resp)


def _get_json_or_none_on_404(url, opener):
    """PubChem returns HTTP 404 for 'compound/assay not found' — a normal miss, not an error."""
    try:
        return _get_json(url, opener)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        raise


def chembl_id_for_inchikey(ik, opener=urllib.request.urlopen):
    """First ChEMBL molecule id whose standard InChIKey matches `ik`, or None if unknown."""
    q = urllib.parse.urlencode({"molecule_structures__standard_inchi_key": ik})
    data = _get_json(f"{_CHEMBL_BASE}/molecule.json?{q}", opener)
    mols = data.get("molecules") or []
    return mols[0]["molecule_chembl_id"] if mols else None


def chembl_activities(chembl_id, opener=urllib.request.urlopen, page_limit=1000):
    """All ChEMBL activity records for a molecule id, following pagination to the end."""
    q = urllib.parse.urlencode({"molecule_chembl_id": chembl_id, "limit": page_limit})
    url = f"{_CHEMBL_BASE}/activity.json?{q}"
    out = []
    while url:
        data = _get_json(url, opener)
        out.extend(data.get("activities") or [])
        nxt = (data.get("page_meta") or {}).get("next")
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
    return out


def pubchem_cid_for_inchikey(ik, opener=urllib.request.urlopen):
    """First PubChem CID whose InChIKey matches `ik`, or None if PubChem has never seen it."""
    data = _get_json_or_none_on_404(
        f"{_PUBCHEM_BASE}/compound/inchikey/{urllib.parse.quote(ik)}/cids/JSON", opener)
    cids = ((data or {}).get("IdentifierList") or {}).get("CID") or []
    return cids[0] if cids else None


def pubchem_assays(cid, opener=urllib.request.urlopen):
    """PubChem assay summary for a CID -> [{'outcome', 'assay_name'}, ...] (empty if none).

    The assaysummary endpoint returns a column-oriented Table; this reads the two columns the
    bucket/verdict logic needs by name, so PubChem re-ordering or adding columns cannot shift
    the wrong field into place.
    """
    data = _get_json_or_none_on_404(
        f"{_PUBCHEM_BASE}/compound/cid/{cid}/assaysummary/JSON", opener)
    table = (data or {}).get("Table")
    if not table:
        return []
    cols = table["Columns"]["Column"]
    oi, ni = cols.index("Activity Outcome"), cols.index("Assay Name")
    return [{"outcome": row["Cell"][oi], "assay_name": row["Cell"][ni]}
            for row in table.get("Row", [])]


def _clean_bindingdb_smiles(smile):
    """BindingDB SMILES carry a trailing CXSMILES tail (' |r|', ' |c:..|'); the leading token
    is the plain structure RDKit needs to derive an InChIKey."""
    return (smile or "").split()[0] if smile else ""


def bindingdb_target_ligands(uniprot, opener=urllib.request.urlopen, cutoff=10_000_000):
    """Every ligand BindingDB has measured against `uniprot`, normalised with its provenance.

    Returns [{'smiles','monomerid','affinity_type','affinity','pmid','doi'}, ...]. `cutoff` is
    BindingDB's server-side affinity ceiling in nM; kept high so nothing is filtered out there —
    our own nM rule (bindingdb_record_to_activity) decides usability. The response's top-level key
    is BindingDB's own (mis-spelled) wrapper, so it is read by position, not by name, and the typo
    cannot break the parse.
    """
    q = urllib.parse.urlencode({"uniprot": uniprot, "cutoff": cutoff,
                                "response": "application/json"})
    data = _get_json(f"{_BINDINGDB_REST}/getLigandsByUniprots?{q}", opener)
    wrapper = next(iter(data.values())) if isinstance(data, dict) and data else {}
    out = []
    for a in (wrapper.get("affinities") or []):
        out.append({"smiles": _clean_bindingdb_smiles(a.get("smile")),
                    "monomerid": a.get("monomerid"), "affinity_type": a.get("affinity_type"),
                    "affinity": a.get("affinity"), "pmid": a.get("pmid"), "doi": a.get("doi")})
    return out


def build_bindingdb_index(uniprots=CYP51_UNIPROTS, opener=urllib.request.urlopen):
    """{InChIKey -> [BindingDB records]} for the CYP51 target(s) — the fetch-once join table.

    BindingDB is queried target-first, so this whole index is built with one call per UniProt and
    then matched against every candidate locally by InChIKey (the same join key the compound-first
    probes use). A record whose SMILES RDKit cannot standardise is skipped — it cannot be matched
    to a candidate. Returns {} when given no targets. Pass the result to lookup(bdb_index=...).
    """
    index = {}
    for up in uniprots:
        for r in bindingdb_target_ligands(up, opener):
            ik = inchikey(r["smiles"])
            if ik is not None:
                index.setdefault(ik, []).append(r)
    return index


def lookup(smiles, sources=("chembl", "pubchem"), opener=urllib.request.urlopen, bdb_index=None):
    """Full precedent check for one SMILES across the requested sources.

    Canonicalises to an InChIKey, resolves it in each source, fetches records, and renders one
    merged summary. Returns summarise()'s keys plus 'inchikey', 'chembl_id', 'pubchem_cid', and
    'bindingdb' (the matched BindingDB records, with PMID/DOI provenance — empty when none).
    Distinguishes 'unparseable' (RDKit rejected the SMILES) from 'no-precedent' (a real molecule
    with no enzyme record found). BindingDB is target-first, so it is not fetched per molecule
    here: pass a pre-built `bdb_index` (build_bindingdb_index) and this matches the InChIKey into
    it locally. Network errors propagate to the caller, which isolates them per-molecule.
    """
    ik = inchikey(smiles)
    if ik is None:
        return {"inchikey": None, "chembl_id": None, "pubchem_cid": None,
                "verdict": "unparseable", "tally": None, "phenotypic": None,
                "phenotypic_organisms": {"fungal": [], "nonfungal": []},
                "bindingdb": [], "n_off_target": 0, "n_records": 0}

    buckets = _empty_buckets()
    chembl_id = pubchem_cid = None
    chembl_records = []                          # kept so the phenotypic organism breakdown (#79)
    bindingdb_hits = []                          # matched BindingDB records (provenance for #78)
    if "chembl" in sources:                      # can name the organisms behind the verdict
        chembl_id = chembl_id_for_inchikey(ik, opener)
        if chembl_id:
            chembl_records = chembl_activities(chembl_id, opener)
            buckets = merge_buckets(buckets, chembl_buckets(chembl_records))
    if "pubchem" in sources:
        pubchem_cid = pubchem_cid_for_inchikey(ik, opener)
        if pubchem_cid:
            buckets = merge_buckets(buckets, pubchem_buckets(pubchem_assays(pubchem_cid, opener)))
    if "bindingdb" in sources and bdb_index:
        bindingdb_hits = bdb_index.get(ik, [])
        if bindingdb_hits:
            buckets = merge_buckets(buckets, bindingdb_buckets(bindingdb_hits))

    summary = render(buckets)
    summary.update({"inchikey": ik, "chembl_id": chembl_id, "pubchem_cid": pubchem_cid,
                    "phenotypic_organisms": phenotypic_organism_breakdown(chembl_records),
                    "bindingdb": bindingdb_hits})
    return summary
