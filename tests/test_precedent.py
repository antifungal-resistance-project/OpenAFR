"""Tests for the assay-precedent check (openafr/precedent.py).

What these lock:
  * the candidate-facing verdict is inactives.compound_verdict's judgement, only re-labelled --
    so the two cannot drift, and 'contrary evidence wins' still holds end to end,
  * a missing record is 'no-precedent', NEVER 'untested' -- absence is not evidence,
  * records are grouped by target: only the CYP51 enzyme gives an on-target verdict, whole-cell
    MICs are demoted to a phenotypic verdict, unrelated targets are a novelty count only,
  * the near-neighbour signal fires only above the (high) similarity floor and picks the best,
  * the network layer's parsing/pagination is exercised through an injected fake opener, so the
    whole path is tested with NO network.
"""
import io
import json
import urllib.error

from rdkit import Chem

from openafr import precedent as P


def rec(**kw):
    """A ChEMBL activity record with only the fields the rules read."""
    base = {"standard_value": None, "standard_relation": None, "standard_units": None,
            "activity_comment": None, "pchembl_value": None, "target_chembl_id": None}
    base.update(kw)
    return base


def enzyme(**kw):
    return rec(target_chembl_id=next(iter(P.CYP51_TARGETS)), **kw)


# --- the verdict is a re-label of the shared inactives judgement -----------------------

def test_no_records_is_no_precedent_not_untested():
    """A missing record must never read as 'untested' -- a failing assay may just be undeposited."""
    v = P.candidate_verdict([])
    assert v["verdict"] == "no-precedent" and v["tally"] is None


def test_measured_inactive_enzyme_records_are_tested_inactive():
    v = P.candidate_verdict([enzyme(standard_value="64", standard_relation="=",
                                    standard_units="ug.mL-1"),
                             enzyme(activity_comment="Not Active")])
    assert v["verdict"] == "tested-inactive" and v["tally"]["inactive"] == 2


def test_one_active_record_flips_to_tested_active():
    """Contrary evidence still wins, inherited straight from inactives.compound_verdict."""
    v = P.candidate_verdict([enzyme(standard_value="128", standard_relation="=",
                                    standard_units="ug.mL-1"),
                             enzyme(pchembl_value="6.5")])
    assert v["verdict"] == "tested-active"


def test_weak_measurement_is_ambiguous():
    v = P.candidate_verdict([enzyme(standard_value="16", standard_relation="=",
                                    standard_units="ug.mL-1")])
    assert v["verdict"] == "tested-ambiguous"


# --- grouping by target ---------------------------------------------------------------

def test_grouping_separates_enzyme_phenotypic_and_offtarget():
    on, pheno, off = P.group_by_target([
        enzyme(standard_value="64", standard_relation="=", standard_units="ug.mL-1"),
        rec(target_chembl_id="CHEMBL999", standard_value="128", standard_relation="=",
            standard_units="ug.mL-1", target_organism="Candida albicans"),  # fungal whole-cell MIC
        rec(target_chembl_id="CHEMBL555", standard_value="5", standard_units="nM"),  # off-target
        rec(),                                               # no target -> off-target
    ])
    assert len(on) == 1 and len(pheno) == 1 and len(off) == 2


def test_summary_reports_phenotypic_and_offtarget_separately():
    s = P.summarise([
        rec(target_chembl_id="CHEMBL999", standard_value="128", standard_relation="=",
            standard_units="ug.mL-1", target_organism="Candida albicans"),
        rec(target_chembl_id="CHEMBL555", standard_value="5", standard_units="nM"),
    ])
    assert s["verdict"] == "no-precedent"                    # nothing on the enzyme itself
    assert s["phenotypic"]["verdict"] == "tested-inactive"   # the fungal whole-cell MIC, demoted
    assert s["n_off_target"] == 1 and s["n_records"] == 2


# --- organism-aware phenotypic bucketing (#79) ----------------------------------------

def test_is_fungal_organism_recognises_fungi_and_rejects_others():
    assert P.is_fungal_organism("Candida albicans")
    assert P.is_fungal_organism("Aspergillus fumigatus")
    assert P.is_fungal_organism("Saccharomyces cerevisiae")
    assert P.is_fungal_organism("Cryptococcus neoformans")
    # a bacterium and a human cell line are NOT antifungal evidence
    assert not P.is_fungal_organism("Staphylococcus aureus")
    assert not P.is_fungal_organism("Escherichia coli")
    assert not P.is_fungal_organism("Homo sapiens")
    # an unrecorded organism cannot be confirmed fungal -> False (becomes off-target, not fungal)
    assert not P.is_fungal_organism("")


def test_record_organism_prefers_target_then_assay_organism():
    assert P.record_organism(rec(target_organism="Candida glabrata")) == "Candida glabrata"
    assert P.record_organism(rec(assay_organism="Aspergillus niger")) == "Aspergillus niger"
    assert P.record_organism(rec()) == ""


def test_nonfungal_and_unspecified_mic_go_to_offtarget_not_phenotypic():
    """The #79 fix: an antibacterial or organism-less µg/mL MIC is NOT an antifungal signal."""
    on, pheno, off = P.group_by_target([
        rec(target_chembl_id="CHEMBL111", standard_value="4", standard_relation="=",
            standard_units="ug.mL-1", target_organism="Staphylococcus aureus"),   # antibacterial
        rec(target_chembl_id="CHEMBL222", standard_value="8", standard_relation="=",
            standard_units="ug.mL-1"),                                            # no organism
        rec(target_chembl_id="CHEMBL333", standard_value="2", standard_relation="=",
            standard_units="ug.mL-1", target_organism="Candida krusei"),          # fungal
    ])
    assert len(on) == 0 and len(pheno) == 1 and len(off) == 2
    assert pheno[0]["target_organism"] == "Candida krusei"


def test_antibacterial_mic_no_longer_reads_as_phenotypic_active():
    """Flucloxacillin's real bug: an 'active' µg/mL MIC vs S. aureus once flagged phenotypic."""
    s = P.summarise([
        rec(target_chembl_id="CHEMBL111", standard_value="1", standard_relation="=",
            standard_units="ug.mL-1", target_organism="Staphylococcus aureus"),
    ])
    assert s["phenotypic"] is None                           # not an antifungal signal at all
    assert s["n_off_target"] == 1
    assert s["phenotypic_organisms"] == {"fungal": [], "nonfungal": ["Staphylococcus aureus"]}


def test_summary_carries_the_fungal_organisms_behind_the_verdict():
    s = P.summarise([
        rec(target_chembl_id="CHEMBL1", standard_value="64", standard_relation="=",
            standard_units="ug.mL-1", target_organism="Candida albicans"),
        rec(target_chembl_id="CHEMBL2", standard_value="64", standard_relation="=",
            standard_units="ug.mL-1", assay_organism="Aspergillus fumigatus"),
        rec(target_chembl_id="CHEMBL3", standard_value="4", standard_relation="=",
            standard_units="ug.mL-1", target_organism="Homo sapiens"),           # reclassified out
    ])
    assert s["phenotypic"]["verdict"] == "tested-inactive"
    assert s["phenotypic_organisms"]["fungal"] == ["Aspergillus fumigatus", "Candida albicans"]
    assert s["phenotypic_organisms"]["nonfungal"] == ["Homo sapiens"]


# --- near-neighbour soft signal -------------------------------------------------------

def _fp(smiles):
    return inactives_fp(smiles)


def inactives_fp(smiles):
    from openafr import inactives
    return inactives.fingerprint(Chem.MolFromSmiles(smiles))


def test_nearest_neighbour_fires_only_above_threshold_and_picks_best():
    fluconazole = "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F"
    query = inactives_fp(fluconazole)
    neighbours = [
        {"id": "SELF", "fp": inactives_fp(fluconazole), "verdict": "tested-inactive"},   # ~1.0
        {"id": "FAR", "fp": inactives_fp("c1ccccc1"), "verdict": "tested-active"},        # far
    ]
    best = P.nearest_tested_neighbour(query, neighbours)
    assert best["id"] == "SELF" and best["similarity"] >= P.NEIGHBOUR_TANIMOTO


def test_nearest_neighbour_returns_none_when_nothing_is_close():
    query = inactives_fp("OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F")
    far = [{"id": "FAR", "fp": inactives_fp("CCO"), "verdict": "tested-inactive"}]
    assert P.nearest_tested_neighbour(query, far) is None


def test_load_neighbours_defaults_verdict_and_skips_unparseable(tmp_path):
    p = tmp_path / "tested.smi"
    p.write_text("OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F\tFLU\n"      # no verdict -> default
                 "c1ccccc1\tBEN\ttested-active\n"                  # explicit verdict
                 "not-a-smiles\tBAD\n")                            # unparseable -> skipped
    ns = P.load_neighbours(str(p))
    assert [n["id"] for n in ns] == ["FLU", "BEN"]
    assert ns[0]["verdict"] == "tested-inactive"                   # the default
    assert ns[1]["verdict"] == "tested-active"


def test_fingerprint_is_none_for_unparseable():
    assert P.fingerprint("not-a-smiles") is None
    assert P.fingerprint("c1ccccc1") is not None


# --- canonicalisation -----------------------------------------------------------------

def test_inchikey_is_none_for_unparseable_smiles():
    assert P.inchikey("not a smiles") is None


def test_inchikey_is_tautomer_stable_join_key():
    """Two writings of the same molecule must yield the same key, or the ChEMBL join fails."""
    assert P.inchikey("c1ccccc1") == P.inchikey("C1=CC=CC=C1")


# --- network layer through an injected fake opener (no real network) ------------------

class _FakeOpener:
    """Maps URL substrings to canned JSON payloads, mimicking urllib.request.urlopen."""
    def __init__(self, routes):
        self.routes = routes
        self.calls = []

    def __call__(self, url, timeout=None):
        self.calls.append(url)
        for needle, payload in self.routes.items():
            if needle in url:
                if isinstance(payload, Exception):
                    raise payload
                return io.BytesIO(json.dumps(payload).encode())
        raise AssertionError(f"unexpected URL: {url}")


def _http_404(url="x"):
    return urllib.error.HTTPError(url, 404, "Not Found", {}, None)


def _pubchem_table(rows):
    """A PubChem assaysummary Table with the two columns the reader needs, plus decoys."""
    cols = ["AID", "Activity Outcome", "Target GeneID", "Assay Name", "Assay Type"]
    return {"Table": {"Columns": {"Column": cols},
                      "Row": [{"Cell": ["1", oc, "", name, "Confirmatory"]}
                              for oc, name in rows]}}


def test_lookup_end_to_end_with_fake_chembl():
    opener = _FakeOpener({
        "molecule.json": {"molecules": [{"molecule_chembl_id": "CHEMBL42"}]},
        "activity.json": {"activities": [
            enzyme(standard_value="100", standard_relation=">", standard_units="ug.mL-1")],
            "page_meta": {"next": None}},
    })
    out = P.lookup("OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F", sources=("chembl",), opener=opener)
    assert out["chembl_id"] == "CHEMBL42"
    assert out["verdict"] == "tested-inactive"


def test_lookup_reports_no_precedent_when_chembl_has_never_seen_it():
    opener = _FakeOpener({"molecule.json": {"molecules": []}})
    out = P.lookup("c1ccccc1", sources=("chembl",), opener=opener)
    assert out["chembl_id"] is None and out["verdict"] == "no-precedent"


def test_lookup_reports_unparseable_without_hitting_network():
    opener = _FakeOpener({})                                  # any call would AssertionError
    out = P.lookup("not a smiles", opener=opener)
    assert out["verdict"] == "unparseable" and opener.calls == []


# --- PubChem adapter: outcome mapping, name-based buckets, Unspecified ignored ---------

def test_pubchem_outcome_mapping():
    assert P.classify_pubchem_outcome("Active") == "active"
    assert P.classify_pubchem_outcome("Inactive") == "inactive"
    assert P.classify_pubchem_outcome("Inconclusive") == "ambiguous"
    # 'Unspecified' is NOT read as non-inhibition -- the depositor simply did not call it.
    assert P.classify_pubchem_outcome("Unspecified") == "ignore"
    assert P.classify_pubchem_outcome("") == "ignore"


def test_pubchem_bucket_by_assay_name():
    assert P.pubchem_bucket("CYP51 sterol 14-alpha demethylase inhibition") == "on_target"
    assert P.pubchem_bucket("Antifungal activity against Candida albicans") == "phenotypic"
    assert P.pubchem_bucket("qHTS for inhibitors of RGS12 GoLoco Motif") == "off_target"
    # a bare 'demethylase' must NOT sweep in unrelated (e.g. histone) demethylases
    assert P.pubchem_bucket("Histone lysine demethylase KDM assay") == "off_target"


def test_pubchem_buckets_verdict_and_ignored_unspecified():
    b = P.pubchem_buckets([
        {"outcome": "Inactive", "assay_name": "CYP51 sterol 14-demethylase inhibition"},
        {"outcome": "Unspecified", "assay_name": "In vitro antifungal activity vs Candida"},
        {"outcome": "Active", "assay_name": "qHTS for some unrelated target"},
    ])
    s = P.render(b)
    assert s["verdict"] == "tested-inactive"          # the on-target enzyme outcome
    assert s["phenotypic"]["verdict"] == "record-inconclusive"  # Unspecified -> ignored
    assert s["n_off_target"] == 1


def test_pubchem_assays_reads_columns_by_name():
    opener = _FakeOpener({"assaysummary/JSON": _pubchem_table(
        [("Inactive", "CYP51 inhibition"), ("Active", "other")])})
    rows = P.pubchem_assays(1234, opener=opener)
    assert rows == [{"outcome": "Inactive", "assay_name": "CYP51 inhibition"},
                    {"outcome": "Active", "assay_name": "other"}]


def test_pubchem_cid_404_is_a_miss_not_an_error():
    opener = _FakeOpener({"cids/JSON": _http_404()})
    assert P.pubchem_cid_for_inchikey("AAAAAAAAAAAAAA-AAAAAAAAAA-A", opener=opener) is None


def test_pubchem_assays_404_returns_empty():
    opener = _FakeOpener({"assaysummary/JSON": _http_404()})
    assert P.pubchem_assays(1234, opener=opener) == []


# --- the two sources merge at the tally level ------------------------------------------

def test_merge_contrary_evidence_wins_across_sources():
    """ChEMBL says inactive, PubChem says active -> the merged on-target verdict is active."""
    chembl = P.chembl_buckets([enzyme(standard_value="128", standard_relation="=",
                                      standard_units="ug.mL-1")])   # inactive
    pubchem = P.pubchem_buckets([{"outcome": "Active", "assay_name": "CYP51 inhibition"}])
    merged = P.render(P.merge_buckets(chembl, pubchem))
    assert merged["verdict"] == "tested-active"


def test_lookup_merges_both_sources():
    opener = _FakeOpener({
        "molecule.json": {"molecules": [{"molecule_chembl_id": "CHEMBL42"}]},
        "activity.json": {"activities": [], "page_meta": {"next": None}},   # ChEMBL: nothing
        "cids/JSON": {"IdentifierList": {"CID": [3365]}},
        "assaysummary/JSON": _pubchem_table([("Inactive", "CYP51 demethylase inhibition")]),
    })
    out = P.lookup("OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F", opener=opener)
    assert out["chembl_id"] == "CHEMBL42" and out["pubchem_cid"] == 3365
    assert out["verdict"] == "tested-inactive"        # supplied entirely by PubChem


def test_activities_follow_pagination():
    # Seed URL carries the molecule id + limit; the follow-up URL carries the offset from
    # page_meta.next. Routing on those two distinct needles proves the loop chases 'next'.
    opener = _FakeOpener({
        "molecule_chembl_id=CHEMBL42": {
            "activities": [enzyme(activity_comment="Not Active")],
            "page_meta": {"next": "/chembl/api/data/activity.json?offset=1"}},
        "offset=1": {"activities": [enzyme(activity_comment="Active")],
                     "page_meta": {"next": None}},
    })
    acts = P.chembl_activities("CHEMBL42", opener=opener)
    assert len(acts) == 2                                     # both pages collected


# --- BindingDB adapter: nM affinities, target-first index, merge into the verdict (#78) ----

def _bdb_response(recs):
    """A BindingDB getLigandsByUniprots JSON body, including its own mis-spelled wrapper key
    and the ' |r|' CXSMILES tail real records carry, so the parser is exercised as it ships."""
    return {"getLindsByUniprotsResponse": {"affinities": [
        {"query": "Lanosterol 14-alpha demethylase", "monomerid": m, "smile": s + " |r|",
         "affinity_type": t, "affinity": a, "pmid": p, "doi": "10.0/x"}
        for (s, m, t, a, p) in recs]}}


def test_bindingdb_affinity_parsing_relations_and_junk():
    assert P.parse_bindingdb_affinity("250") == ("=", 250.0)
    assert P.parse_bindingdb_affinity(">10000") == (">", 10000.0)
    assert P.parse_bindingdb_affinity("<=0.5") == ("<=", 0.5)
    # unreadable affinity is dropped, never guessed
    assert P.parse_bindingdb_affinity("n/a") is None
    assert P.parse_bindingdb_affinity(None) is None


def test_bindingdb_nM_is_active_only_never_inactive():
    """An nM binding value is active evidence when potent, else ignored — never a fake inactive."""
    potent = P.bindingdb_buckets([{"affinity": "250", "affinity_type": "IC50"}])
    assert potent["on_target"]["active"] == 1
    weak = P.bindingdb_buckets([{"affinity": "50000000", "affinity_type": "Kd"}])   # > 10 uM
    assert weak["on_target"] == {"active": 0, "inactive": 0, "ambiguous": 0, "ignore": 1}
    junk = P.bindingdb_buckets([{"affinity": "n/a", "affinity_type": "IC50"}])
    assert sum(junk["on_target"].values()) == 0                  # unparseable is not tallied


def test_bindingdb_target_ligands_parse_strips_cxsmiles_tail():
    opener = _FakeOpener({"getLigandsByUniprots": _bdb_response(
        [("c1ccccc1", "111", "IC50", "42", "999")])})
    recs = P.bindingdb_target_ligands("P10613", opener=opener)
    assert recs == [{"smiles": "c1ccccc1", "monomerid": "111", "affinity_type": "IC50",
                     "affinity": "42", "pmid": "999", "doi": "10.0/x"}]


def test_build_bindingdb_index_keys_by_inchikey():
    opener = _FakeOpener({"getLigandsByUniprots": _bdb_response(
        [("c1ccccc1", "111", "IC50", "42", "999")])})
    index = P.build_bindingdb_index(uniprots=("P10613",), opener=opener)
    ik = P.inchikey("c1ccccc1")
    assert list(index) == [ik] and index[ik][0]["monomerid"] == "111"


def test_lookup_bindingdb_flips_no_precedent_to_tested_active():
    """A candidate ChEMBL/PubChem never saw, but BindingDB measured potent, becomes tested-active
    with its PMID provenance carried through — the file-drawer catch #78 exists to add."""
    smi = "c1ccccc1"
    index = P.build_bindingdb_index(uniprots=("P10613",), opener=_FakeOpener(
        {"getLigandsByUniprots": _bdb_response([(smi, "111", "IC50", "42", "999")])}))
    opener = _FakeOpener({
        "molecule.json": {"molecules": []},                     # ChEMBL: nothing
        "cids/JSON": _http_404(),                               # PubChem: nothing
    })
    out = P.lookup(smi, sources=("chembl", "pubchem", "bindingdb"), opener=opener, bdb_index=index)
    assert out["verdict"] == "tested-active"
    assert [h["pmid"] for h in out["bindingdb"]] == ["999"]
    assert out["n_records"] == 1


def test_lookup_ignores_bindingdb_when_source_not_requested():
    smi = "c1ccccc1"
    index = P.build_bindingdb_index(uniprots=("P10613",), opener=_FakeOpener(
        {"getLigandsByUniprots": _bdb_response([(smi, "111", "IC50", "42", "999")])}))
    opener = _FakeOpener({"molecule.json": {"molecules": []}})
    out = P.lookup(smi, sources=("chembl",), opener=opener, bdb_index=index)
    assert out["verdict"] == "no-precedent" and out["bindingdb"] == []
