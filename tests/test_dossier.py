"""Tests for the candidate-dossier fusion + priority table (openafr/dossier.py).

What these lock:
  * build_one assembles every documented schema key from already-computed profile dicts,
    and the result JSON-round-trips (it is machine-readable by construction),
  * the A/B/C priority table fires per rule: tested-inactive / PAINS / unrankable / soft
    negative / profile error all force C; the full all-of-A conjunction earns A; a single
    failed A condition drops to B,
  * informational reasons TRAVEL regardless of tier — an azole warhead is always logged,
    a tested-active is tagged 'reference' (never automatic A),
  * build_all preserves geometry-rank order independent of priority, stamps provenance,
    and computes a valid priority with precedent omitted (the offline path),
  * the whole thing is deterministic (sort_keys JSON is stable across runs).
"""
import json

import pytest

from openafr import dossier


# --- minimal already-computed profile fixtures (the shapes the real profilers return) -----

def _admet(ro5=0, alerts=0, pains=None):
    return {"mw": 341.4, "clogp": 3.1, "tpsa": 62.0, "ro5_violations": ro5,
            "veber_ok": True, "esol_logs": -4.2, "pains": pains or [],
            "structural_alerts": [] if alerts == 0 else [("BRENK", "x")],
            "n_alerts": alerts, "herg_risk_factors": False}


def _cyp(warhead=None):
    flagged = list(dossier._cyp.PANEL) if warhead else []
    return {"azole_warhead": warhead, "warhead_strength": "strong" if warhead else None,
            "azine_nitrogen": False, "n_high": len(flagged) if warhead else 0,
            "n_moderate": 0, "flagged_isoforms": flagged}


def _synth(availability="catalogued-drug"):
    return {"sa_score": 3.1, "sa_band": "easy", "availability": availability, "pubchem_cid": 1}


def _geom(rank=1, best_fe=2.4, n_passing=7, n_poses=9):
    return {"rank": rank, "best_fe": best_fe, "n_passing": n_passing,
            "n_poses": n_poses, "best_score": -8.1}


# a real, parseable SMILES so inchikey() works without a network call (fluconazole)
SMI = "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F"


def _build(prec=None, near=None, error=None, top_cutoff=1, **kw):
    a = kw.pop("admet", _admet())
    c = kw.pop("cyp", _cyp())
    s = kw.pop("synth", _synth())
    g = kw.pop("geometry", _geom())
    return dossier.build_one("OAFR-1", SMI, g, admet=a, cyp=c, synth=s,
                             prec=prec, near=near, error=error, top_cutoff=top_cutoff, **kw)


# --- schema ---------------------------------------------------------------------------------

def test_schema_keys_and_json_roundtrip():
    d = _build(prec={"verdict": "no-precedent", "phenotypic": None, "n_off_target": 4,
                     "inchikey": "IK"})
    for k in ("name", "smiles", "inchikey", "geometry", "precedent", "developability",
              "selectivity", "synthesis", "priority", "priority_reasons",
              "selection_hypothesis"):
        assert k in d, "missing key: %s" % k
    assert "_near" not in d                 # transient priority input is stripped
    assert json.loads(json.dumps(d)) == d   # machine-readable, round-trips


# --- automatic C ---------------------------------------------------------------------------

def test_tested_inactive_forces_C():
    d = _build(prec={"verdict": "tested-inactive", "phenotypic": None, "n_off_target": 0})
    assert d["priority"] == "C"
    assert any("tested-inactive" in r for r in d["priority_reasons"])


def test_pains_forces_C():
    d = _build(admet=_admet(pains=["some PAINS motif"]))
    assert d["priority"] == "C"
    assert any("PAINS" in r for r in d["priority_reasons"])


def test_unrankable_forces_C():
    d = _build(geometry=_geom(rank=None, best_fe=None, n_passing=0, n_poses=0))
    assert d["priority"] == "C"
    assert any("unrankable" in r for r in d["priority_reasons"])


def test_soft_negative_neighbour_forces_C():
    near = {"tanimoto": 0.91, "name": "known-dud", "verdict": "tested-inactive"}
    d = _build(near=near)
    assert d["priority"] == "C"
    assert any("soft negative" in r for r in d["priority_reasons"])


def test_profile_error_forces_C():
    d = dossier.build_one("bad", "not-a-smiles", _geom(), admet=None, cyp=None, synth=None,
                          error="unparseable SMILES", top_cutoff=1)
    assert d["priority"] == "C"
    assert any("profile error" in r for r in d["priority_reasons"])


# --- automatic A and the B fallback --------------------------------------------------------

def test_full_conjunction_earns_A():
    d = _build(cyp=_cyp(warhead=None),  # a clean, non-azole candidate
               prec={"verdict": "no-precedent", "phenotypic": None, "n_off_target": 4})
    assert d["priority"] == "A"
    assert any(r.startswith("priority A") for r in d["priority_reasons"])


# --- coordinator-identity gating (RESULTS_coordinator.md) -----------------------------------

def _coord(verdict, kind="nitrile", aromatic=6.9):
    return {"reported_kind": kind, "in_mechanism": verdict == "in-mechanism",
            "reported_distance": 2.42, "aromatic_distance": aromatic, "verdict": verdict}


def test_coordinator_key_in_schema():
    d = _build(coord=_coord("in-mechanism", kind="aromatic-ring", aromatic=2.42))
    assert "coordinator" in d
    assert json.loads(json.dumps(d)) == d


def test_out_of_mechanism_blocks_A_not_C():
    # a non-aromatic coordinator blocks automatic-A and is flagged, but does NOT force C:
    # the active luliconazole also coordinates via a nitrile, so it is unvalidated, not disproven
    d = _build(coord=_coord("out-of-mechanism", kind="nitrile", aromatic=6.9))
    assert d["priority"] == "B"
    assert any("out-of-mechanism" in r for r in d["priority_reasons"])
    assert any("in-mechanism coordination" in r for r in d["priority_reasons"])


def test_dual_mode_blocks_A_drops_to_B():
    # rank inflated by a nitrile but the aromatic N reaches the band -> not C, but not A
    d = _build(coord=_coord("dual-mode", kind="nitrile", aromatic=2.79))
    assert d["priority"] == "B"
    assert any("dual-mode" in r for r in d["priority_reasons"])
    assert any("in-mechanism coordination" in r for r in d["priority_reasons"])


def test_in_mechanism_still_earns_A():
    d = _build(coord=_coord("in-mechanism", kind="aromatic-ring", aromatic=2.42))
    assert d["priority"] == "A"


def test_absent_coordinator_is_backward_compatible():
    # no coordinator supplied -> the in-mechanism A-condition passes, behaviour unchanged
    d = _build()
    assert d["priority"] == "A"
    assert d["coordinator"] is None


@pytest.mark.parametrize("break_it,expect_in", [
    (dict(geometry=_geom(rank=99)), "top-decile"),          # rank outside top decile
    (dict(admet=_admet(ro5=3)), "clean drug-likeness"),      # too many Ro5 violations
    (dict(synth=_synth(availability="no-pubchem-record")), "physically obtainable"),
    (dict(geometry=_geom(n_passing=1)), "convergent poses"), # not convergent
])
def test_single_failed_A_condition_drops_to_B(break_it, expect_in):
    d = _build(cyp=_cyp(warhead=None),
               prec={"verdict": "no-precedent", "phenotypic": None, "n_off_target": 4},
               top_cutoff=10, **break_it)
    assert d["priority"] == "B"
    assert any(expect_in in r for r in d["priority_reasons"])


# --- reasons travel regardless of tier -----------------------------------------------------

def test_azole_warhead_always_logged():
    d = _build(cyp=_cyp(warhead="triazole"),
               prec={"verdict": "no-precedent", "phenotypic": None, "n_off_target": 4})
    assert any("azole warhead" in r for r in d["priority_reasons"])


def test_known_active_tagged_reference_not_A():
    d = _build(cyp=_cyp(warhead=None),
               prec={"verdict": "tested-active", "phenotypic": None, "n_off_target": 4})
    assert d["priority"] != "A"          # a known active is real, but not a discovery
    assert any("reference" in r for r in d["priority_reasons"])


# --- batch: ordering, provenance, offline path ---------------------------------------------

def test_build_all_preserves_rank_order_and_stamps_provenance():
    recs = [
        {"name": "r1", "smiles": SMI, "geometry": _geom(rank=1, best_fe=2.1),
         "admet": _admet(), "cyp": _cyp(), "synth": _synth(),
         "prec": {"verdict": "no-precedent", "phenotypic": None, "n_off_target": 0}},
        {"name": "r2", "smiles": SMI, "geometry": _geom(rank=2, best_fe=2.6),
         "admet": _admet(), "cyp": _cyp(), "synth": _synth(),
         "prec": {"verdict": "tested-inactive", "phenotypic": None, "n_off_target": 0}},
    ]
    prov = {"protocol_sha256": "abc", "protocol_unmodified": True}
    out = dossier.build_all(recs, prov)
    assert [d["name"] for d in out] == ["r1", "r2"]      # input order preserved
    assert out[1]["priority"] == "C"                      # r2 was tested-inactive
    assert all(d["provenance"] is prov for d in out)


def test_offline_path_has_null_precedent_and_valid_priority():
    d = _build(prec=None, cyp=_cyp(warhead=None))
    assert d["precedent"] is None
    assert d["priority"] in ("A", "B", "C")   # 'not a known active' holds when prec is None


def test_deterministic_json():
    kw = dict(prec={"verdict": "no-precedent", "phenotypic": None, "n_off_target": 4})
    a = json.dumps(_build(**kw), sort_keys=True)
    b = json.dumps(_build(**kw), sort_keys=True)
    assert a == b
