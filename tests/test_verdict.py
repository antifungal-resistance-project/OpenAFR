"""Tests for the concordance-only diagnostic verdict contract (openafr/verdict.py).

These lock the four honesty rules the contract (docs/DIAGNOSTIC_VERDICT_CONTRACT.md) exists
to enforce:
  * the verdict is one of exactly four categorical values, with NO probability field;
  * NO_KNOWN_MARKER is emitted only when every panel residue was confidently read, and is
    never a susceptibility claim (the RUO scope string is always present);
  * an UNCHARACTERIZED_VARIANT (non-synonymous, non-panel) gets an explicit verdict, and its
    token is tagged 'uncharacterized', not dropped;
  * an uncalled PANEL residue forces UNRESOLVED, never NO_KNOWN_MARKER, but a detected panel
    hit still wins over an unrelated uncalled position;
  * FKS1 is detection-only -- its structural field is always null; ERG11 carries one.
Everything runs against synthetic caller outputs -- no reads, no tools, no network.
"""
from openafr import verdict as v


def _call(tokens=(), panel_hits=(), uncalled_panel=()):
    """A minimal caller-output dict with the keys both callers expose."""
    return {"tokens": list(tokens), "panel_hits": list(panel_hits),
            "uncalled_panel": list(uncalled_panel)}


def test_panel_hit_is_resistance_marker_detected():
    out = v.azole_verdict(_call(tokens=["Y132F"], panel_hits=["Y132F"]))
    assert out["verdict"] == v.RESISTANCE_MARKER_DETECTED
    assert out["resolution"] == "resolved"
    assert out["called_tokens"] == [{"token": "Y132F", "class": "known_resistance"}]


def test_nonpanel_change_is_uncharacterized_not_silence():
    out = v.azole_verdict(_call(tokens=["Y123H"]))
    assert out["verdict"] == v.UNCHARACTERIZED_VARIANT
    assert out["called_tokens"] == [{"token": "Y123H", "class": "uncharacterized"}]


def test_no_tokens_is_no_known_marker_not_susceptible():
    out = v.azole_verdict(_call())
    assert out["verdict"] == v.NO_KNOWN_MARKER
    assert out["called_tokens"] == []
    # NO_KNOWN_MARKER must never read as a clinical susceptibility call.
    assert "not a clinical determination" in out["scope"]
    assert "suscept" not in out["verdict"].lower()


def test_uncalled_panel_residue_forces_unresolved():
    # A key marker position could not be read and nothing positive was seen elsewhere:
    # we cannot rule resistance out -> UNRESOLVED, never NO_KNOWN_MARKER (rule 4).
    out = v.azole_verdict(_call(uncalled_panel=[132]))
    assert out["verdict"] == v.UNRESOLVED
    assert out["resolution"] == "unresolved"
    assert out["resolution_note"] == "panel_residue_uncalled:132"


def test_uncharacterized_token_with_uncalled_panel_is_unresolved_but_surfaces_token():
    # An uncharacterized change is visible AND a panel residue is unreadable: overall
    # UNRESOLVED (panel check incomplete), yet the observed token is still reported.
    out = v.azole_verdict(_call(tokens=["T123I"], uncalled_panel=[143]))
    assert out["verdict"] == v.UNRESOLVED
    assert out["called_tokens"] == [{"token": "T123I", "class": "uncharacterized"}]


def test_panel_hit_wins_over_unrelated_uncalled_position():
    # A confirmed marker is a confirmed marker even if another position was uncalled.
    out = v.azole_verdict(_call(tokens=["Y132F"], panel_hits=["Y132F"],
                                uncalled_panel=[143]))
    assert out["verdict"] == v.RESISTANCE_MARKER_DETECTED
    assert out["resolution"] == "resolved"


def test_hard_unresolved_reason_takes_precedence():
    out = v.azole_verdict(unresolved_reason="consensus length 1576 is not a multiple of 3")
    assert out["verdict"] == v.UNRESOLVED
    assert out["called_tokens"] == []
    assert "multiple of 3" in out["resolution_note"]


def test_no_probability_field_anywhere():
    for out in (
        v.azole_verdict(_call(tokens=["Y132F"], panel_hits=["Y132F"])),
        v.echinocandin_verdict(_call()),
        v.azole_verdict(unresolved_reason="coverage_gap"),
    ):
        keys = " ".join(out.keys()).lower()
        assert "prob" not in keys and "score" not in keys and "mic" not in keys


def test_fks1_is_detection_only_structural_always_null():
    out = v.echinocandin_verdict(_call(tokens=["S639F"], panel_hits=["S639F"]))
    assert out["drug_class"] == "echinocandin"
    assert out["structural"] is None
    # Even if a structural verdict is (wrongly) offered, FKS1 must drop it.
    forced = v.build_verdict("FKS1", _call(), structural={"pocket": "whatever"})
    assert forced["structural"] is None


def test_erg11_carries_structural_when_provided():
    pocket = {"reaches_iron": True}
    out = v.azole_verdict(_call(tokens=["Y132F"], panel_hits=["Y132F"]), structural=pocket)
    assert out["drug_class"] == "azole"
    assert out["structural"] == pocket


def test_scope_and_provenance_present():
    out = v.azole_verdict(_call(), provenance={"reference_accession": "XM_..."})
    assert out["scope"] == v.RUO_SCOPE
    assert out["provenance"]["reference_accession"] == "XM_..."


def test_unknown_gene_rejected():
    try:
        v.build_verdict("CDR1", _call())
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for unknown gene")


# --- PPV tier: FKS1 v3 (#137) -------------------------------------------------
# When the FKS1 caller reports `core_panel_hits`, only those high-PPV hits DETECT resistance;
# a panel hit that is present but not core (D642Y, M690I) abstains rather than over-calls.

def _fks1_call(tokens=(), panel_hits=(), core_panel_hits=(), uncalled_panel=()):
    return {"tokens": list(tokens), "panel_hits": list(panel_hits),
            "core_panel_hits": list(core_panel_hits), "uncalled_panel": list(uncalled_panel)}


def test_core_panel_hit_detects_resistance():
    out = v.echinocandin_verdict(
        _fks1_call(tokens=["S639F"], panel_hits=["S639F"], core_panel_hits=["S639F"]))
    assert out["verdict"] == v.RESISTANCE_MARKER_DETECTED
    assert out["called_tokens"] == [{"token": "S639F", "class": "known_resistance"}]


def test_low_ppv_only_hit_abstains_not_detects():
    # D642Y is a tagged panel hit but low-PPV -> UNCHARACTERIZED_VARIANT (honest abstain),
    # NOT resistance-detected and NOT a susceptibility claim. The token is still surfaced,
    # labelled 'known_low_ppv' so the reader sees it was recognised and deliberately not read
    # as resistance.
    out = v.echinocandin_verdict(
        _fks1_call(tokens=["D642Y"], panel_hits=["D642Y"], core_panel_hits=[]))
    assert out["verdict"] == v.UNCHARACTERIZED_VARIANT
    assert out["called_tokens"] == [{"token": "D642Y", "class": "known_low_ppv"}]


def test_core_plus_low_ppv_still_detects():
    # A carrier of both a core marker and a low-PPV position still detects resistance.
    out = v.echinocandin_verdict(
        _fks1_call(tokens=["S639F", "D642Y"], panel_hits=["S639F", "D642Y"],
                   core_panel_hits=["S639F"]))
    assert out["verdict"] == v.RESISTANCE_MARKER_DETECTED
    classes = {ct["token"]: ct["class"] for ct in out["called_tokens"]}
    assert classes == {"S639F": "known_resistance", "D642Y": "known_low_ppv"}


def test_erg11_has_no_tier_every_panel_hit_detects():
    # ERG11 exposes no core_panel_hits key -> backward-compatible: every panel hit detects.
    out = v.azole_verdict(_call(tokens=["Y132F"], panel_hits=["Y132F"]))
    assert out["verdict"] == v.RESISTANCE_MARKER_DETECTED
    assert out["called_tokens"] == [{"token": "Y132F", "class": "known_resistance"}]
