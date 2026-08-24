"""Known-answer tests for the end-to-end front door (scripts/score_molecule.py).

score_molecule chains triage -> prep -> dock -> geometry into one command. The dock step
needs a VM, so the same discipline the recaller used applies: the *decision* and
*interpretation* logic is pure and tested here, leaving only the subprocess dock for the
real run. These lock:
  * plan(): only in-envelope molecules earn a docking run; nothing is dropped.
  * classify_geometry(): each band label vs. the validated actives' iron-bound band.
  * render_molecule(): the pending / not-docked / no-pose / scored branches + caveats.
"""
import pytest

rdkit = pytest.importorskip("rdkit")  # triage/plan need rdkit; skip cleanly without it

from scripts.applicability_triage import classify, load_panel
from scripts.score_molecule import (
    DOCKABLE,
    HEADER_CAVEATS,
    VALIDATED_FE_BAND,
    classify_geometry,
    molecule_record,
    plan,
    render_molecule,
)

# Reuse the same reference molecules the triage tests pin.
FLUCONAZOLE = "C1=CC(=C(C=C1F)F)C(CN2C=NC=N2)(CN3C=NC=N3)O"           # near-known (in panel)
POSACONAZOLE = ("CC[C@@H]([C@H](C)O)N1C(=O)N(C=N1)C2=CC=C(C=C2)N3CCN(CC3)C4=CC=C(C=C4)"
                "OC[C@H]5C[C@](OC5)(CN6C=NC=N6)C7=C(C=C(C=C7)F)F")     # out-of-domain (51 heavy)
BENZENE = "c1ccccc1"                                                    # out-of-scope
VT_1598 = ("C1=CC(=CC=C1COC2=CC=C(C=C2)C#CC3=CN=C(C=C3)C([C@](CN4C=NN=N4)"
           "(C5=C(C=C(C=C5)F)F)O)(F)F)C#N")                            # in-envelope-novel


@pytest.fixture(scope="module")
def panel():
    return load_panel()


# ---- plan(): the docking gate ------------------------------------------------------------

def test_plan_docks_only_in_envelope_molecules(panel):
    rows = [("vt1598", VT_1598), ("benzene", BENZENE), ("posaconazole", POSACONAZOLE),
            ("fluconazole", FLUCONAZOLE), ("junk", "not-a-molecule")]
    planned = plan(rows, panel)
    dock = {name: d for name, _s, _t, d in planned}
    assert dock["vt1598"] is True          # in-envelope-novel -> dock
    assert dock["fluconazole"] is True     # near-known is still inside the envelope -> dock
    assert dock["benzene"] is False        # out-of-scope
    assert dock["posaconazole"] is False   # out-of-domain
    assert dock["junk"] is False           # unparseable


def test_plan_preserves_order_and_drops_nothing(panel):
    rows = [("benzene", BENZENE), ("vt1598", VT_1598)]
    planned = plan(rows, panel)
    assert [p[0] for p in planned] == ["benzene", "vt1598"]  # order kept, nothing dropped


def test_dockable_set_matches_triage_verdicts():
    # Guard against the two rule sources drifting apart.
    assert DOCKABLE == {"in-envelope-novel", "near-known"}


# ---- classify_geometry(): interpretation vs. the validated band --------------------------

def test_geometry_within_validated_band_coordinates():
    lo, hi = VALIDATED_FE_BAND
    label, coord = classify_geometry(best_fe=(lo + hi) / 2, n_passing=1)
    assert coord is True
    assert "within the validated actives' band" in label


def test_geometry_tighter_than_band_still_coordinates():
    label, coord = classify_geometry(best_fe=VALIDATED_FE_BAND[0] - 0.2, n_passing=1)
    assert coord is True
    assert "tighter" in label


def test_geometry_loose_edge_is_iron_bound_but_flagged():
    # between the actives band and the coordination cutoff: iron-bound but noted as loose.
    label, coord = classify_geometry(best_fe=VALIDATED_FE_BAND[1] + 0.05, n_passing=1)
    assert coord is True
    assert "loose edge" in label


def test_geometry_beyond_cutoff_is_not_iron_bound():
    label, coord = classify_geometry(best_fe=3.5, n_passing=0)
    assert coord is False
    assert "NOT iron-bound" in label


def test_geometry_no_nitrogen_pose_is_uninterpretable():
    label, coord = classify_geometry(best_fe=None, n_passing=0)
    assert coord is False
    assert "nothing to interpret" in label


# ---- render_molecule(): the four branches ------------------------------------------------

def test_render_out_of_scope_is_not_docked(panel):
    t = classify(BENZENE, panel)
    out = render_molecule("benzene", t, geom=None, pending=False)
    assert "not docked" in out
    assert "makes no claim" in out


def test_render_in_envelope_pending_shows_pending(panel):
    t = classify(VT_1598, panel)
    out = render_molecule("vt1598", t, geom=None, pending=True)
    assert "PENDING" in out


def test_render_docked_no_pose_is_reported_not_pending(panel):
    t = classify(VT_1598, panel)
    out = render_molecule("vt1598", t, geom=None, pending=False)
    assert "NO usable pose" in out
    assert "PENDING" not in out


def test_render_scored_molecule_shows_geometry_and_number(panel):
    t = classify(VT_1598, panel)
    geom = {"best_fe": 2.60, "n_passing": 3, "n_poses": 9, "best_score": -8.4}
    out = render_molecule("vt1598", t, geom=geom, pending=False)
    assert "within the validated actives' band" in out
    assert "2.600 A" in out
    assert "3/9" in out
    assert "-8.4 kcal/mol" in out


# ---- molecule_record(): the JSON mirror of render_molecule -------------------------------

def test_record_out_of_scope_is_not_docked(panel):
    t = classify(BENZENE, panel)
    rec = molecule_record("benzene", BENZENE, t, geom=None, pending=False)
    assert rec["triage"] == "out-of-scope"
    assert rec["dockable"] is False
    assert rec["geometry_status"] == "not-docked"
    assert rec["geometry"] is None


def test_record_in_envelope_pending(panel):
    t = classify(VT_1598, panel)
    rec = molecule_record("vt1598", VT_1598, t, geom=None, pending=True)
    assert rec["dockable"] is True
    assert rec["geometry_status"] == "pending"
    assert rec["geometry"] is None


def test_record_docked_no_pose(panel):
    t = classify(VT_1598, panel)
    rec = molecule_record("vt1598", VT_1598, t, geom=None, pending=False)
    assert rec["geometry_status"] == "no-pose"


def test_record_scored_carries_geometry_and_interpretation(panel):
    t = classify(VT_1598, panel)
    geom = {"best_fe": 2.60, "n_passing": 3, "n_poses": 9, "best_score": -8.4}
    rec = molecule_record("vt1598", VT_1598, t, geom=geom, pending=False)
    assert rec["geometry_status"] == "scored"
    g = rec["geometry"]
    assert g["closest_n_fe"] == 2.60
    assert g["iron_bound_poses"] == 3
    assert g["n_poses"] == 9
    assert g["vina_score"] == -8.4
    assert g["coordinates_iron"] is True
    assert "within the validated actives' band" in g["interpretation"]


def test_record_is_json_serializable(panel):
    import json
    t = classify(VT_1598, panel)
    rec = molecule_record("vt1598", VT_1598, t, geom=None, pending=True)
    json.dumps(rec)  # must not raise


# ---- honesty caveats are carried, not just in docstrings ---------------------------------

def test_header_carries_the_standing_caveats():
    lower = HEADER_CAVEATS.lower()
    assert "hypothesis" in lower           # a rank is a hypothesis, not a hit
    assert "auc" in lower                  # the durable metric
    assert "not sufficient" in lower       # decoy ceiling: coordinating is necessary-not-sufficient
    assert "validate_gate2" in HEADER_CAVEATS
