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
            standard_units="ug.mL-1"),                       # whole-cell MIC, not the enzyme
        rec(target_chembl_id="CHEMBL555", standard_value="5", standard_units="nM"),  # off-target
        rec(),                                               # no target -> off-target
    ])
    assert len(on) == 1 and len(pheno) == 1 and len(off) == 2


def test_summary_reports_phenotypic_and_offtarget_separately():
    s = P.summarise([
        rec(target_chembl_id="CHEMBL999", standard_value="128", standard_relation="=",
            standard_units="ug.mL-1"),
        rec(target_chembl_id="CHEMBL555", standard_value="5", standard_units="nM"),
    ])
    assert s["verdict"] == "no-precedent"                    # nothing on the enzyme itself
    assert s["phenotypic"]["verdict"] == "tested-inactive"   # the whole-cell MIC, demoted
    assert s["n_off_target"] == 1 and s["n_records"] == 2


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
