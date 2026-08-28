"""Known-answer tests for the protonation/tautomer axis (issue #84).

Two layers are locked:
  - the deterministic enumeration (scripts/enumerate_protomers.py): a molecule with a single
    relevant microstate at pH 7.4 (protonation stable across 6.4-8.4 AND one tautomer) is
    `protomer_explicit`; one with a titratable group or multiple tautomers is not, and the
    reference (obabel -p 7.4) microspecies is always a member of the enumerated set;
  - the per-candidate classifier (scripts/candidate_protomer.py) frozen in
    work/PREREGISTRATION_candidate_protomer.md, BEFORE any per-microstate distance was read:
    ROBUST iff every microstate coordinates AND the spread is within STATE_TOL;
    MICROSTATE-DEPENDENT the moment a microstate fails to coordinate OR the spread exceeds it.
"""
from scripts.candidate_protomer import MIN_COORD_SEEDS, STATE_TOL, classify
from scripts.enumerate_protomers import audit_one


# ---- enumeration -------------------------------------------------------------

def test_rigid_nonionizable_molecule_is_protomer_explicit():
    # fluconazole: no group titrating across pH 6.4-8.4, tautomer-locked -> single microstate
    rec = audit_one("OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F", "fluconazole", 16)
    assert rec["protomer_explicit"] is True
    assert rec["n_states"] == 1
    assert rec["n_protonation_states"] == 1


def test_reference_microspecies_is_always_a_member():
    # whatever the verdict, the obabel -p 7.4 form (the docked scorecard form) must be present
    rec = audit_one("Cc1[nH]c(=O)c(cc1-c1ccc2nccn2c1)C#N", "olprinone", 16)
    assert rec["reference_smiles"] in [s["smiles"] for s in rec["states"]]
    assert sum(s["is_reference"] for s in rec["states"]) == 1


def test_multi_tautomer_molecule_is_not_explicit_and_enumerates():
    # olprinone (2-pyridone/2-hydroxypyridine + more): several tautomers -> multi-state
    rec = audit_one("Cc1[nH]c(=O)c(cc1-c1ccc2nccn2c1)C#N", "olprinone", 16)
    assert rec["protomer_explicit"] is False
    assert rec["n_states"] > 1


def test_state_cap_is_respected():
    # flucloxacillin enumerates past the cap; the audit must not exceed MAX_STATES
    rec = audit_one(
        "Cc1onc(c1C(=O)N[C@H]1[C@H]2SC(C)(C)[C@@H](N2C1=O)C(O)=O)-c1c(F)cccc1Cl",
        "flucloxacillin", 16)
    assert rec["n_states"] <= 16
    # reference survives the cap
    assert rec["reference_smiles"] in [s["smiles"] for s in rec["states"]]


# ---- classifier --------------------------------------------------------------

def test_explicit_candidates_are_not_docked_here():
    # a protomer-explicit molecule never reaches classify(); guard the contract instead
    rec = audit_one("OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F", "fluconazole", 16)
    assert rec["protomer_explicit"] is True  # -> builder emits EXPLICIT without any dock


def test_all_states_coordinate_and_tight_is_robust():
    # every microstate reaches the iron within the noise band -> the rank is microstate-robust
    verdict, spread = classify([("st01", 2.541, 5), ("st02", 2.563, 3)])
    assert verdict == "ROBUST"
    assert round(spread, 3) == 0.022


def test_a_state_that_fails_to_coordinate_is_microstate_dependent():
    # one tautomer never reaches the iron (0/5) -> the rank depends on the microstate
    verdict, _ = classify([("st01", 2.548, 3), ("st02", 3.342, 0), ("st03", 2.72, 4)])
    assert verdict == "MICROSTATE-DEPENDENT"


def test_wide_state_spread_is_microstate_dependent_even_if_all_coordinate():
    # every microstate coordinates but they disagree by > STATE_TOL on where the iron sits
    verdict, spread = classify([("st01", 2.20, 5), ("st02", 2.20 + STATE_TOL + 0.01, 5)])
    assert verdict == "MICROSTATE-DEPENDENT"
    assert spread > STATE_TOL


def test_seed_count_below_minimum_counts_as_non_coordinating():
    # a microstate at MIN_COORD_SEEDS-1 is treated as not reaching the iron
    verdict, _ = classify([("st01", 2.5, 5), ("st02", 2.5, MIN_COORD_SEEDS - 1)])
    assert verdict == "MICROSTATE-DEPENDENT"
