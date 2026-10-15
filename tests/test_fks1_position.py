"""Tests for the FKS1 novel-variant POSITION CONTEXT derivation (openafr/fks1_position.py).

Position context (#4 Rung 1) escalates a bare UNCHARACTERIZED_VARIANT with geometry-free
geography -- which hot-spot, offset from the anchor, nearest known marker -- WITHOUT making
any fit/phenotype/probability claim. These lock the frozen non-claims
(work/PREREGISTRATION_fks1_position_context.md):

  * a CORE panel marker returns None (it is already RESISTANCE_MARKER_DETECTED, not novel);
  * a known-but-LOW-PPV letter (D642Y, M690I) is its own honest tier, not "novel";
  * a genuinely novel letter at a panel residue (S639L, W691C) is the strongest tier;
  * a novel change elsewhere in a window is adjacent-to / background, by distance only;
  * anything outside every hot-spot window, or unparseable, returns None;
  * the derivation stays consistent with the pinned caller panel (no independent letters).

Pure: no reads, no tools, no network.
"""
from openafr import fks1_position as fp
from openafr import fks1_caller as F


def test_core_panel_marker_returns_none():
    # S639F is a canonical high-PPV marker: already detected, so no position context.
    assert fp.position_context("S639F") is None
    assert fp.position_context("W691L") is None


def test_low_ppv_letter_is_its_own_tier_not_novel():
    # D642Y / M690I ARE the low-PPV panel letters; the verdict abstains on them. They must
    # not be mislabelled "novel substitution" -- that would misread the token.
    for tok in ("D642Y", "M690I"):
        ctx = fp.position_context(tok)
        assert ctx is not None
        assert ctx["tier"] == fp.KNOWN_LOW_PPV_MARKER
        assert "low-PPV" in fp.describe(ctx)
        assert "novel" not in fp.describe(ctx)


def test_novel_letter_at_panel_residue_is_strongest_tier():
    ctx = fp.position_context("S639L")  # novel letter L at the S639 anchor residue
    assert ctx["tier"] == fp.NOVEL_AT_PANEL_RESIDUE
    assert ctx["window"] == "HS1"
    assert ctx["anchor"] == "S639"
    assert ctx["anchor_distance"] == 0
    assert ctx["at_panel_residue"] is True
    assert ctx["panel_letters_at_position"] == ["F", "P", "Y"]
    assert "F/P/Y" in fp.describe(ctx)


def test_benchmark_circular_variant_w691c_is_novel_at_anchor():
    # W691C is the anti-circularity case (reported only in the benchmark -> untagged by the
    # caller). Position context still locates it honestly: novel letter at the W691 anchor.
    ctx = fp.position_context("W691C")
    assert ctx["tier"] == fp.NOVEL_AT_PANEL_RESIDUE
    assert ctx["window"] == "HS3"
    assert ctx["anchor"] == "W691"
    assert ctx["anchor_distance"] == 0


def test_adjacent_vs_background_by_distance_only():
    near = fp.position_context("L636P")  # 1 aa from marker position 635
    assert near["tier"] == fp.ADJACENT_TO_PANEL
    assert near["nearest_core_marker"] == "635"
    assert near["distance_to_nearest"] == 1
    assert "adjacent" in fp.describe(near)
    # A window-background position far from any core marker (HS1 spans 635-643; 643 is >2 aa
    # from the nearest core position 639/642-excluded... 642 is low-PPV so excluded as core).
    bg = fp.position_context("P643A")
    assert bg["window"] == "HS1"
    assert bg["tier"] == fp.HOTSPOT_BACKGROUND
    assert "background" in fp.describe(bg)
    assert "after the S639 anchor" in fp.describe(bg)


def test_describe_of_none_is_empty():
    assert fp.describe(None) == ""
    assert fp.describe({}) == ""


def test_outside_every_window_returns_none():
    assert fp.position_context("A100V") is None       # not in any hot-spot
    assert fp.position_context("M1I") is None


def test_unparseable_token_returns_none():
    for bad in ("", "junk", "639F", "SX639F", None):
        assert fp.position_context(bad) is None


def test_low_ppv_positions_excluded_from_core_markers():
    # Distance-to-nearest must measure distance to a marker we would actually CALL, so the
    # low-PPV positions (642 D642Y, 690 M690I) are not core positions.
    assert 642 not in fp._CORE_POSITIONS
    assert 690 not in fp._CORE_POSITIONS
    assert {635, 639, 691, 1354} <= fp._CORE_POSITIONS


def test_panel_stays_consistent_with_caller():
    # The derivation reuses the caller's panel; a letter the caller tags core must never be
    # offered position context, and vice-versa -- the two can never drift.
    for pos, (wt, letters) in F.RESISTANCE_PANEL.items():
        for mut in letters:
            tok = f"{wt}{pos}{mut}"
            if F.is_core_panel_token(tok):
                assert fp.position_context(tok) is None
