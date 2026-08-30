"""Known-answer tests for the consolidated scorecard logic (scripts/scorecard.py, #88).

The scorecard's whole job is to be honest when the axes are seen together: it must
demote a candidate the moment ANY measured axis says so, must never let a control's
known-inhibitor status read as a strike, and must never let a card look complete by
omission. These lock exactly those rules on the pure `derive_concerns`, plus the
invariant that the not-yet-checked axis set is non-empty (so every card carries it).
"""
from scripts.scorecard import NOT_YET_CHECKED, derive_concerns

# A candidate clean on every measured axis (lersivirine's real shape).
CLEAN_CF = {"reaches_fungal_iron": True, "stable": True, "selectivity_margin": 2.07}
CLEAN_AM = {"flags": "clean"}
CLEAN_SY = {"sa_band": "easy"}
CLEAN_DM = {"flags": []}
CLEAN_PR = {"verdict": "no-precedent", "phenotypic": "—"}


def test_all_axes_clean_yields_no_concern():
    assert derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR, False) == []


def test_seed_unstable_is_flagged():
    # LDN-27219: coordinates in only 2/5 seeds -> stable False -> demoted
    cf = {**CLEAN_CF, "stable": False}
    concerns = derive_concerns(cf, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR, False)
    assert any("seed-unstable" in c for c in concerns)


def test_non_selective_hit_is_flagged():
    # nolatrexed: reaches the human iron as readily as the fungal one (margin ~ 0)
    cf = {**CLEAN_CF, "selectivity_margin": -0.06}
    concerns = derive_concerns(cf, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR, False)
    assert any("selectivity window" in c for c in concerns)


def test_positive_margin_is_not_flagged():
    cf = {**CLEAN_CF, "selectivity_margin": 0.01}  # strictly positive -> a window exists
    concerns = derive_concerns(cf, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR, False)
    assert not any("selectivity window" in c for c in concerns)


def test_admet_and_domain_flags_surface_verbatim():
    am = {"flags": "Ro5x2,insoluble"}
    dm = {"flags": ["coordinator-fragile"]}
    concerns = derive_concerns(CLEAN_CF, am, CLEAN_SY, dm, CLEAN_PR, False)
    assert "ADMET: Ro5x2,insoluble" in concerns
    assert "domain: coordinator-fragile" in concerns


def test_clean_admet_string_is_not_a_concern():
    # the literal "clean" sentinel must not be turned into a concern
    concerns = derive_concerns(CLEAN_CF, {"flags": "clean"}, CLEAN_SY, CLEAN_DM, CLEAN_PR, False)
    assert not any(c.startswith("ADMET:") for c in concerns)


def test_control_phenotypic_active_is_not_a_strike():
    # a known azole active in a cell assay is the point of a control, not a demotion
    pr = {"verdict": "no-precedent", "phenotypic": "tested-active"}
    ctrl = derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, pr, True)
    novel = derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, pr, False)
    assert not any("phenotypic" in c for c in ctrl)
    assert any("phenotypic" in c for c in novel)  # but a NOVEL hit must be flagged


def test_known_enzyme_inhibitor_flagged_even_for_control():
    # fluconazole/voriconazole: enzyme tested-active must always surface (it is the label)
    pr = {"verdict": "tested-active", "phenotypic": "—"}
    concerns = derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, pr, True)
    assert any("known CYP51 inhibitor" in c for c in concerns)


def test_not_reaching_iron_is_flagged():
    cf = {**CLEAN_CF, "reaches_fungal_iron": False}
    concerns = derive_concerns(cf, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR, False)
    assert any("does not reach fungal iron" in c for c in concerns)


def test_resistance_lost_is_flagged_for_novel_candidate():
    # vatalanib: 5/5 stable in WT but collapses to 2/5 in the C. auris Y132F pocket
    res = {"verdict": "LOST"}
    concerns = derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR, False, res)
    assert any("Y132F" in c for c in concerns)


def test_resistance_retained_is_not_a_concern():
    res = {"verdict": "RETAINED"}
    concerns = derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR, False, res)
    assert not any("Y132F" in c for c in concerns)


def test_resistance_lost_is_not_a_strike_against_a_control():
    # controls are the Y132F-defeated azoles; their mutant behavior is a robustness check
    res = {"verdict": "LOST"}
    concerns = derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR, True, res)
    assert not any("Y132F" in c for c in concerns)


def test_resistance_unchecked_adds_no_concern():
    # backward-compatible default: no resistance record -> no concern
    assert derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR, False) == []


def test_stereo_isomer_dependent_is_flagged_for_novel_candidate():
    # a rank that flips with an unspecified stereocenter is not trustworthy
    st = {"verdict": "ISOMER-DEPENDENT"}
    concerns = derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR,
                               False, None, st)
    assert any("stereochemistry" in c for c in concerns)


def test_stereo_explicit_and_robust_are_not_concerns():
    for verdict in ("EXPLICIT", "ROBUST"):
        st = {"verdict": verdict}
        concerns = derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR,
                                   False, None, st)
        assert not any("stereochemistry" in c for c in concerns)


def test_stereo_isomer_dependent_is_not_a_strike_against_a_control():
    # the ambiguous molecules are the racemic azole controls, docked to exercise the axis
    st = {"verdict": "ISOMER-DEPENDENT"}
    concerns = derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR,
                               True, None, st)
    assert not any("stereochemistry" in c for c in concerns)


def test_stereo_unchecked_adds_no_concern():
    assert derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR, False) == []


def test_protomer_microstate_dependent_is_flagged_for_novel_candidate():
    # a rank that collapses/shifts across pH-7.4 protonation/tautomer states is not trustworthy
    pt = {"verdict": "MICROSTATE-DEPENDENT"}
    concerns = derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR,
                               False, None, None, pt)
    assert any("protonation/tautomer" in c for c in concerns)


def test_protomer_explicit_and_robust_are_not_concerns():
    for verdict in ("EXPLICIT", "ROBUST"):
        pt = {"verdict": verdict}
        concerns = derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR,
                                   False, None, None, pt)
        assert not any("protonation/tautomer" in c for c in concerns)


def test_protomer_microstate_dependent_is_not_a_strike_against_a_control():
    # the multi-state molecules include the tautomeric azole controls, docked to exercise the axis
    pt = {"verdict": "MICROSTATE-DEPENDENT"}
    concerns = derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR,
                               True, None, None, pt)
    assert not any("protonation/tautomer" in c for c in concerns)


def test_protomer_unchecked_adds_no_concern():
    # backward-compatible default: no protomer record -> no concern
    assert derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR, False) == []


def test_cyp_type_ii_high_is_flagged_for_novel_candidate():
    # L-838417 / olprinone: an azole warhead -> HIGH type-II across the human CYP panel
    cyp = {"n_high": 5}
    concerns = derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR,
                               False, None, None, None, None, cyp)
    assert any("human-CYP" in c and "CYP3A4" in c for c in concerns)


def test_cyp_moderate_only_is_not_a_concern():
    # a weaker azine/pharmacophore profile (no HIGH isoform) is reported, not a demotion
    cyp = {"n_high": 0, "n_moderate": 5}
    concerns = derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR,
                               False, None, None, None, None, cyp)
    assert not any("human-CYP" in c for c in concerns)


def test_cyp_type_ii_is_not_a_strike_against_a_control():
    # ketoconazole IS the textbook CYP3A4 inhibitor; its HIGH type-II validates the axis
    cyp = {"n_high": 5}
    concerns = derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR,
                               True, None, None, None, None, cyp)
    assert not any("human-CYP" in c for c in concerns)


def test_cyp_unchecked_adds_no_concern():
    # backward-compatible default: no cyp record -> no concern
    assert derive_concerns(CLEAN_CF, CLEAN_AM, CLEAN_SY, CLEAN_DM, CLEAN_PR, False) == []


def test_unchecked_axis_set_is_nonempty_so_every_card_carries_it():
    # the honesty guard: there is always at least one axis a card must admit it skipped
    assert len(NOT_YET_CHECKED) >= 1
    assert all(isinstance(v, str) and v for v in NOT_YET_CHECKED.values())
