"""Tests for the one genotype -> verdict entrypoint (openafr/interpret.py, issue #134).

verdict.py (#133) is unit-tested against synthetic caller dicts; these instead drive the
FULL path -- raw sequence / variant-list input -> real re-caller -> #133 verdict object --
against the pinned references, so a break anywhere in the wiring shows up here. They lock:
  * both genes route to the right caller and the right drug-class;
  * a real substitution sequence yields the correct four-value verdict + panel tagging,
    identical to the equivalent variant-list input;
  * the callers' honest hard-failures (ERG11 frameshift, FKS1 un-anchorable CDS) become an
    UNRESOLVED verdict carrying the reason -- never a raised exception;
  * a missing/refused FKS1 hot-spot window forces UNRESOLVED at its panel residue rather than
    leaking out as NO_KNOWN_MARKER (contract rule 4);
  * a variant list is taken as covering the panel unless uncalled= declares a gap;
  * CYP51A and other unknown genes, and malformed input, are refused loudly;
  * provenance is auto-stamped from the pinned reference and structural passes through ERG11
    but is dropped for FKS1.
Runs against the pinned references on disk; no network, no external tools.
"""
import pytest

from openafr import fks1_caller
from openafr import interpret as I
from openafr import recaller
from openafr import verdict as v


@pytest.fixture(scope="module")
def erg11_ref():
    return recaller.load_reference()


@pytest.fixture(scope="module")
def fks1_ref():
    return fks1_caller.load_reference()


def _mutate_codon(cds, pos1, new_codon):
    i = (pos1 - 1) * 3
    assert len(new_codon) == 3
    return cds[:i] + new_codon + cds[i + 3:]


# --- ERG11 / azole: sequence path -------------------------------------------

def test_erg11_wildtype_cds_is_no_known_marker(erg11_ref):
    ref_cds, _ = erg11_ref
    out = I.interpret("ERG11", cds=ref_cds, reference=erg11_ref)
    assert out["drug_class"] == "azole"
    assert out["gene"] == "ERG11"
    assert out["verdict"] == v.NO_KNOWN_MARKER
    assert out["called_tokens"] == []
    # NO_KNOWN_MARKER is never a susceptibility call.
    assert "not a clinical determination" in out["scope"]


def test_erg11_panel_mutation_cds_is_resistance_marker(erg11_ref):
    ref_cds, _ = erg11_ref
    # Y132F: residue 132 Tyr(TAT/TAC) -> Phe(TTT).
    mutated = _mutate_codon(ref_cds, 132, "TTT")
    out = I.interpret("ERG11", cds=mutated, reference=erg11_ref)
    assert out["verdict"] == v.RESISTANCE_MARKER_DETECTED
    assert out["called_tokens"] == [{"token": "Y132F", "class": "known_resistance"}]


def test_erg11_nonpanel_mutation_cds_is_uncharacterized(erg11_ref):
    ref_cds, ref_prot = erg11_ref
    # A non-panel residue: change residue 200 to something non-synonymous.
    wt = ref_prot[199]
    # pick a codon translating to a different residue; Ala GCT unless wt is already Ala.
    new_codon = "GGT" if wt != "G" else "GCT"
    mutated = _mutate_codon(ref_cds, 200, new_codon)
    out = I.interpret("ERG11", cds=mutated, reference=erg11_ref)
    assert out["verdict"] == v.UNCHARACTERIZED_VARIANT
    assert out["called_tokens"][0]["class"] == "uncharacterized"


def test_erg11_frameshift_cds_becomes_unresolved_not_exception(erg11_ref):
    ref_cds, _ = erg11_ref
    out = I.interpret("ERG11", cds=ref_cds[:-1], reference=erg11_ref)  # drop 1 base
    assert out["verdict"] == v.UNRESOLVED
    assert out["resolution"] == "unresolved"
    assert "multiple of 3" in out["resolution_note"]
    assert out["called_tokens"] == []


def test_erg11_uncalled_panel_codon_cds_is_unresolved(erg11_ref):
    ref_cds, _ = erg11_ref
    # An ambiguous base inside the Y132 codon -> that panel residue is uncalled.
    mutated = _mutate_codon(ref_cds, 132, "TNT")
    out = I.interpret("ERG11", cds=mutated, reference=erg11_ref)
    assert out["verdict"] == v.UNRESOLVED
    assert out["resolution_note"] == "panel_residue_uncalled:132"


# --- ERG11 / azole: variant-list path (must match the sequence path) --------

def test_erg11_variant_list_panel_hit_matches_sequence():
    out = I.interpret("ERG11", variants=["Y132F"])
    assert out["verdict"] == v.RESISTANCE_MARKER_DETECTED
    assert out["called_tokens"] == [{"token": "Y132F", "class": "known_resistance"}]


def test_erg11_variant_list_empty_is_no_known_marker():
    out = I.interpret("ERG11", variants=[])
    assert out["verdict"] == v.NO_KNOWN_MARKER


def test_erg11_variant_list_nonpanel_is_uncharacterized():
    out = I.interpret("ERG11", variants=["t123i"])   # case-insensitive input
    assert out["verdict"] == v.UNCHARACTERIZED_VARIANT
    assert out["called_tokens"] == [{"token": "T123I", "class": "uncharacterized"}]


def test_erg11_variant_list_uncalled_declares_coverage_gap():
    # A targeted panel that could not read residue 132: not wild-type -> UNRESOLVED.
    out = I.interpret("ERG11", variants=[], uncalled=[132])
    assert out["verdict"] == v.UNRESOLVED
    assert out["resolution_note"] == "panel_residue_uncalled:132"


def test_erg11_variant_no_change_token_is_dropped():
    # 'Y132Y' is not a substitution; it must not become an uncharacterized variant.
    out = I.interpret("ERG11", variants=["Y132Y"])
    assert out["verdict"] == v.NO_KNOWN_MARKER
    assert out["called_tokens"] == []


def test_erg11_malformed_variant_token_is_refused():
    with pytest.raises(ValueError, match="unparseable"):
        I.interpret("ERG11", variants=["Y132"])


# --- FKS1 / echinocandin: window + variant paths ----------------------------

def test_fks1_wildtype_windows_is_no_known_marker(fks1_ref):
    ref_cds, _ = fks1_ref
    windows = fks1_caller.slice_windows_from_cds(ref_cds)
    out = I.interpret("FKS1", windows=windows, reference=fks1_ref)
    assert out["drug_class"] == "echinocandin"
    assert out["verdict"] == v.NO_KNOWN_MARKER
    assert out["structural"] is None            # detection-only, always null


def test_fks1_s639f_windows_is_resistance_marker(fks1_ref):
    ref_cds, _ = fks1_ref
    # S639F: residue 639 Ser -> Phe (TTT).
    mutated = _mutate_codon(ref_cds, 639, "TTT")
    windows = fks1_caller.slice_windows_from_cds(mutated)
    out = I.interpret("FKS1", windows=windows, reference=fks1_ref)
    assert out["verdict"] == v.RESISTANCE_MARKER_DETECTED
    assert out["called_tokens"] == [{"token": "S639F", "class": "known_resistance"}]


def test_fks1_missing_hs1_window_forces_unresolved_not_no_marker(fks1_ref):
    ref_cds, _ = fks1_ref
    windows = fks1_caller.slice_windows_from_cds(ref_cds)
    windows.pop("HS1")                            # no data for the S639 hot-spot
    out = I.interpret("FKS1", windows=windows, reference=fks1_ref)
    # Rule 4: an unread panel residue must NOT read as NO_KNOWN_MARKER.
    assert out["verdict"] == v.UNRESOLVED
    assert "639" in out["resolution_note"]


def test_fks1_refused_hs1_window_forces_unresolved(fks1_ref):
    ref_cds, _ = fks1_ref
    windows = fks1_caller.slice_windows_from_cds(ref_cds)
    windows["HS1"] = windows["HS1"][:-1]          # indel -> window refused
    out = I.interpret("FKS1", windows=windows, reference=fks1_ref)
    assert out["verdict"] == v.UNRESOLVED
    assert "639" in out["resolution_note"]


def test_fks1_missing_hs2_does_not_force_unresolved(fks1_ref):
    # HS2 carries no panel markers, so its absence creates no known-marker uncertainty.
    ref_cds, _ = fks1_ref
    windows = fks1_caller.slice_windows_from_cds(ref_cds)
    windows.pop("HS2")
    out = I.interpret("FKS1", windows=windows, reference=fks1_ref)
    assert out["verdict"] == v.NO_KNOWN_MARKER


def test_fks1_unanchorable_cds_becomes_unresolved(fks1_ref):
    out = I.interpret("FKS1", cds="ACGT" * 10, reference=fks1_ref)  # wrong length
    assert out["verdict"] == v.UNRESOLVED
    assert out["resolution"] == "unresolved"


def test_fks1_variant_list_s639f_is_resistance_marker():
    out = I.interpret("FKS1", variants=["S639F"])
    assert out["drug_class"] == "echinocandin"
    assert out["verdict"] == v.RESISTANCE_MARKER_DETECTED
    assert out["structural"] is None


# --- gene routing / provenance / passthrough / misuse -----------------------

def test_cyp51a_is_refused_with_pointer():
    with pytest.raises(ValueError, match="CYP51A"):
        I.interpret("CYP51A", variants=["Y121F"])


def test_unknown_gene_is_refused():
    with pytest.raises(ValueError, match="unknown gene"):
        I.interpret("CDR1", variants=["X1Y"])


def test_requires_exactly_one_input():
    with pytest.raises(ValueError, match="exactly one"):
        I.interpret("ERG11")
    with pytest.raises(ValueError, match="exactly one"):
        I.interpret("ERG11", cds="ATG", variants=["Y132F"])


def test_windows_input_rejected_for_erg11():
    with pytest.raises(ValueError, match="FKS1-only"):
        I.interpret("ERG11", windows={"HS1": "ATG"})


def test_provenance_is_autostamped_and_mergeable():
    out = I.interpret("ERG11", variants=["Y132F"],
                      provenance={"reference_accession": "XM_085597798.1"})
    prov = out["provenance"]
    assert prov["reference_sha256"] == recaller.REFERENCE_CDS_SHA256
    assert prov["reference_accession"] == "XM_085597798.1"
    assert 132 in prov["known_panel"]


def test_structural_passes_through_erg11_but_dropped_for_fks1():
    pocket = {"reaches_iron": False}
    erg = I.interpret("ERG11", variants=["Y132F"], structural=pocket)
    assert erg["structural"] == pocket
    fks = I.interpret("FKS1", variants=["S639F"], structural=pocket)
    assert fks["structural"] is None


def test_no_probability_field_on_the_full_path(erg11_ref):
    ref_cds, _ = erg11_ref
    out = I.interpret("ERG11", cds=ref_cds, reference=erg11_ref)
    keys = " ".join(out.keys()).lower()
    assert "prob" not in keys and "score" not in keys and "mic" not in keys
