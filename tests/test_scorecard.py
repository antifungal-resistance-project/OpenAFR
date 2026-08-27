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


def test_unchecked_axis_set_is_nonempty_so_every_card_carries_it():
    # the honesty guard: there is always at least one axis a card must admit it skipped
    assert len(NOT_YET_CHECKED) >= 1
    assert all(isinstance(v, str) and v for v in NOT_YET_CHECKED.values())
