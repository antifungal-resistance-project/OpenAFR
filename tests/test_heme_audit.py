"""Locks the heme/iron model audit (issue #72, scripts/audit_heme_model.py).

The audited numbers ARE the physical facts the whole criterion rests on. If a
receptor swap, a prep-step change, or a coordinate edit ever moves them, this
test fails loudly instead of letting the foundational measurement drift silently.
"""
import pathlib
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from scripts.audit_heme_model import (  # noqa: E402
    COORD_BAND,
    _verdicts,
    audit_coordination_sphere,
    audit_crystal_reproduction,
    run,
)


def test_all_audit_verdicts_pass():
    v = _verdicts(run())
    assert all(v.values()), f"heme audit regressed: {v}"


def test_coordination_sphere_is_cys_thiolate_ligated_heme():
    s = audit_coordination_sphere()
    # four equatorial porphyrin nitrogens, all in the porphyrin band
    assert len(s["porphyrin_nitrogens"]) == 4
    assert all(1.8 <= d <= 2.2 for d, _ in s["porphyrin_nitrogens"])
    # proximal axial ligand is the conserved Cys470 thiolate at ~2.35 A
    prox = s["proximal_cys_thiolate"]
    assert prox is not None
    assert prox[1] == "SG" and prox[2] == "470"
    assert prox[0] == pytest.approx(2.347, abs=0.01)


def test_distal_axial_site_is_open():
    # the sixth coordination position must be free for the incoming azole
    s = audit_coordination_sphere()
    assert s["distal_axial_ligands"] == []
    assert s["explicit_waters_in_model"] == 0


def test_crystal_azole_reproduces_canonical_coordination():
    c = audit_crystal_reproduction()
    assert c["crystal_min_nfe"] == pytest.approx(2.254, abs=0.01)
    assert COORD_BAND[0] <= c["crystal_min_nfe"] <= COORD_BAND[1]


def test_min_nitrogen_heuristic_selects_the_coordinating_nitrogen():
    c = audit_crystal_reproduction()
    # the closest nitrogen (NAF, the triazole N1) wins by a comfortable margin,
    # so the min-N heuristic is not latching a non-coordinating amine
    assert c["ligand_nitrogen_distances"][0][1] == "NAF"
    assert c["closest_minus_second_A"] >= 0.5
