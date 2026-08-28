"""Known-answer tests for the stereochemistry axis (issue #85).

Two layers are locked:
  - the deterministic enumeration (scripts/enumerate_stereo.py): a molecule whose SMILES
    fixes every stereocenter is `stereo_explicit`; one that leaves a center or an E/Z bond
    unassigned is not, and enumerates to the right isomer count;
  - the per-candidate classifier (scripts/candidate_stereo.py) frozen in
    work/PREREGISTRATION_candidate_stereo.md, BEFORE any per-isomer distance was read:
    ROBUST iff every isomer coordinates AND the spread is within ISO_TOL; ISOMER-DEPENDENT
    the moment an isomer fails to coordinate OR the spread exceeds the tolerance.
"""
from scripts.candidate_stereo import ISO_TOL, MIN_COORD_SEEDS, classify
from scripts.enumerate_stereo import audit_one


# ---- enumeration -------------------------------------------------------------

def test_fully_specified_smiles_is_stereo_explicit():
    # voriconazole's real SMILES: both stereocenters carry @/@@ -> single defined isomer
    rec = audit_one("C[C@@H](c1ncncc1F)[C@](O)(Cn1cncn1)c1ccc(F)cc1F", "voriconazole", 16)
    assert rec["stereo_explicit"] is True
    assert rec["n_isomers"] == 1
    assert rec["n_unassigned"] == 0


def test_achiral_molecule_is_stereo_explicit():
    # fluconazole is achiral -> trivially stereo-explicit, one isomer
    rec = audit_one("OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F", "fluconazole", 16)
    assert rec["stereo_explicit"] is True
    assert rec["n_isomers"] == 1


def test_unspecified_center_is_ambiguous_and_enumerates():
    # ketoconazole supplied without either stereocenter -> not explicit, 4 isomers
    smi = "CC(=O)N1CCN(CC1)c1ccc(OCC2COC(Cn3ccnc3)(c3ccc(Cl)cc3Cl)O2)cc1"
    rec = audit_one(smi, "ketoconazole", 16)
    assert rec["stereo_explicit"] is False
    assert rec["n_unassigned"] == 2
    assert rec["n_isomers"] == 4


# ---- classifier --------------------------------------------------------------

def test_explicit_candidates_are_not_docked_here():
    # a stereo-explicit molecule never reaches classify(); guard the contract instead
    rec = audit_one("OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F", "fluconazole", 16)
    assert rec["stereo_explicit"] is True  # -> builder emits EXPLICIT without any isomer dock


def test_all_isomers_coordinate_and_tight_is_robust():
    # flutrimazole: both isomers coordinate (5/5, 3/5) within the noise band
    verdict, spread = classify([("iso01", 2.541, 5), ("iso02", 2.563, 3)])
    assert verdict == "ROBUST"
    assert round(spread, 3) == 0.022


def test_an_isomer_that_fails_to_coordinate_is_isomer_dependent():
    # ketoconazole: iso03 never reaches the iron (0/5) -> the rank depends on the isomer
    verdict, _ = classify([("iso01", 2.548, 3), ("iso02", 2.72, 2),
                           ("iso03", 3.342, 0), ("iso04", 2.841, 4)])
    assert verdict == "ISOMER-DEPENDENT"


def test_wide_isomer_spread_is_isomer_dependent_even_if_all_coordinate():
    # every isomer coordinates but they disagree by > ISO_TOL on where the iron sits
    verdict, spread = classify([("iso01", 2.20, 5), ("iso02", 2.20 + ISO_TOL + 0.01, 5)])
    assert verdict == "ISOMER-DEPENDENT"
    assert spread > ISO_TOL


def test_seed_count_below_minimum_counts_as_non_coordinating():
    # an isomer at MIN_COORD_SEEDS-1 is treated as not reaching the iron
    verdict, _ = classify([("iso01", 2.5, 5), ("iso02", 2.5, MIN_COORD_SEEDS - 1)])
    assert verdict == "ISOMER-DEPENDENT"
