"""Tests for the FKS1 re-caller deterministic core (openafr/fks1_caller.py).

These lock the same honesty guarantees the ERG11 re-caller has (tests/test_recaller.py),
adapted to the windowed FKS1 model -- every downstream stage fires off the tokens this core
emits, so a silent mis-call here is a silent mis-call in the echinocandin alert:
  * the pinned reference translates to the window anchor residues (S639, R1354) -- the
    numbering check that makes an emitted 'S639F' mean residue 639,
  * a wild-type window emits NO call and is distinguishable from uncalled / missing,
  * a real substitution in a window emits the right ABSOLUTE-position token; S639F/P/Y are
    panel-tagged, an HS2 substitution is emitted but (by design) NOT tagged yet,
  * an ambiguous codon at a hot-spot is UNCALLED, never silently "wild type",
  * a frameshift/indel inside ONE window is REFUSED for that window without discarding the
    other (the whole point of windowing a 5.6 kb gene),
  * a window with no consensus is MISSING, never assumed wild-type,
  * format_call maps results to the (fks1_call, resistance_source) snapshot pair in the
    exact token shape emergence.parse_call consumes.
Runs against the pinned reference on disk; no network, no external tools.
"""
import pathlib

import pytest

from openafr import fks1_caller as F
from openafr import emergence as E


@pytest.fixture(scope="module")
def ref():
    """The pinned (cds, protein) reference, loaded once."""
    return F.load_reference()


def _mutate_codon(cds, pos1, new_codon):
    """Replace the codon for 1-based residue `pos1` with `new_codon` (same length)."""
    i = (pos1 - 1) * 3
    assert len(new_codon) == 3
    return cds[:i] + new_codon + cds[i + 3:]


def _windows(cds):
    """The per-window consensus nt dict for a full-length CDS (mutated or wild-type)."""
    return F.slice_windows_from_cds(cds)


# --- reference integrity / numbering ---------------------------------------

def test_reference_loads_and_hashes(ref):
    cds, protein = ref
    assert len(cds) == 5667
    assert len(protein) == 1888                 # stop stripped
    assert protein.startswith("MSYDNNHNYYDP")


def test_reference_numbering_matches_window_anchors(ref):
    _cds, protein = ref
    # The whole point of the reference: S at 639 and R at 1354 are the wild-type letters
    # the panel/window number their mutations *from*. Drift = every token mis-numbered.
    assert protein[638] == "S"                  # residue 639
    assert protein[1353] == "R"                 # residue 1354
    # and the conserved motifs are where they should be
    assert protein[634:643] == "FLTLSLRDP"      # HS1 residues 635-643
    assert protein[1349:1358] == "DWIRRYTLS"    # HS2 residues 1350-1358


def test_reference_hash_mismatch_is_refused(tmp_path):
    bad = tmp_path / "fks1_cds.fasta"
    bad.write_text(">bad\nATGAAA\n")
    with pytest.raises(F.ReferenceError, match="hash mismatch"):
        F.load_reference(str(bad))


def test_missing_reference_file_is_refused(tmp_path):
    with pytest.raises(F.ReferenceError, match="reference not found"):
        F.load_reference(str(tmp_path / "does_not_exist.fasta"))


def test_reference_file_with_no_sequence_is_refused(tmp_path):
    empty = tmp_path / "fks1_cds.fasta"
    empty.write_text(">header only\n")
    with pytest.raises(F.ReferenceError, match="has no sequence"):
        F.load_reference(str(empty))


def test_reference_with_non_acgt_bases_is_refused(tmp_path):
    bad = tmp_path / "fks1_cds.fasta"
    bad.write_text(">bad\nATGNNNAAA\n")
    with pytest.raises(F.ReferenceError, match="non-ACGT bases"):
        F.load_reference(str(bad))


def test_reference_length_not_multiple_of_three_is_refused(tmp_path):
    bad = tmp_path / "fks1_cds.fasta"
    bad.write_text(">bad\nATGAA\n")
    with pytest.raises(F.ReferenceError, match="not a multiple of 3"):
        F.load_reference(str(bad), allow_hash_mismatch=True)


def test_reference_too_short_for_window_is_refused(tmp_path):
    short = tmp_path / "fks1_cds.fasta"
    short.write_text(">short\n" + "ATG" * 10 + "\n")   # 10 aa, HS1 needs residue 635+
    with pytest.raises(F.ReferenceError, match="too short"):
        F.load_reference(str(short), allow_hash_mismatch=True)


def test_reference_failing_numbering_check_is_refused(tmp_path):
    # Full-length, in-frame, all-ACGT, but residue 639 is no longer Ser. Even with the hash
    # waived, the numbering check -- which can never be waived -- must refuse it.
    real_cds = F.load_reference()[0]
    mangled = _mutate_codon(real_cds, 639, "GGT")       # S639 -> Gly
    bad = tmp_path / "fks1_cds.fasta"
    bad.write_text(">mangled\n" + mangled + "\n")
    with pytest.raises(F.ReferenceError, match="numbering check FAILED"):
        F.load_reference(str(bad), allow_hash_mismatch=True)


# --- wild type: no call, distinguishable from uncalled/missing --------------

def test_call_windows_loads_pinned_reference_when_omitted(ref):
    # With no reference passed, both call_windows and call_window must fall back to the
    # pinned on-disk reference and still emit the right token (exercises the default-load
    # path in both functions).
    mut = _mutate_codon(ref[0], 639, "TTT")             # S639F
    res = F.call_windows(_windows(mut))                  # no reference= kwarg
    assert res["tokens"] == ["S639F"]
    direct = F.call_window(F.FKS1_WINDOWS["HS1"], _windows(mut)["HS1"])
    assert direct["tokens"] == ["S639F"]


def test_wild_type_windows_emit_no_call(ref):
    res = F.call_windows(_windows(ref[0]), reference=ref)
    assert res["tokens"] == []
    assert res["windows"]["HS1"]["status"] == "wild-type"
    assert res["windows"]["HS2"]["status"] == "wild-type"
    call, source = F.format_call(res)
    assert call == ""
    assert source == "sra-fks1-recaller:HS1=wild-type,HS2=wild-type"


# --- real substitutions -----------------------------------------------------

def test_s639f_is_called_and_panel_tagged(ref):
    mut = _mutate_codon(ref[0], 639, "TTT")             # Ser -> Phe
    res = F.call_windows(_windows(mut), reference=ref)
    assert res["tokens"] == ["S639F"]
    assert res["panel_hits"] == ["S639F"]
    call, source = F.format_call(res)
    assert call == "S639F"
    assert source == "sra-fks1-recaller:HS1=called,HS2=wild-type"


@pytest.mark.parametrize("codon,letter", [("CCT", "P"), ("TAT", "Y")])
def test_s639_other_panel_letters_are_tagged(ref, codon, letter):
    mut = _mutate_codon(ref[0], 639, codon)
    res = F.call_windows(_windows(mut), reference=ref)
    assert res["tokens"] == [f"S639{letter}"]
    assert res["panel_hits"] == [f"S639{letter}"]


def test_s639_non_panel_letter_is_emitted_but_not_tagged(ref):
    mut = _mutate_codon(ref[0], 639, "ACT")             # Ser -> Thr (S639T: not in panel)
    res = F.call_windows(_windows(mut), reference=ref)
    assert res["tokens"] == ["S639T"]
    assert res["panel_hits"] == []                      # verbatim, never asserted resistance


def test_hs2_substitution_is_emitted_but_not_tagged(ref):
    # R1354 sits in HS2, whose mutant letters are not yet panel-encoded -- so a real HS2
    # substitution must be REPORTED (so emergence can see it) but NOT tagged as known.
    mut = _mutate_codon(ref[0], 1354, "AGT")            # Arg -> Ser
    res = F.call_windows(_windows(mut), reference=ref)
    assert res["tokens"] == ["R1354S"]
    assert res["panel_hits"] == []
    call, source = F.format_call(res)
    assert call == "R1354S"
    assert source == "sra-fks1-recaller:HS1=wild-type,HS2=called"


def test_non_panel_substitution_in_hs1_emitted_not_tagged(ref):
    mut = _mutate_codon(ref[0], 635, "GGT")             # F635 -> Gly
    res = F.call_windows(_windows(mut), reference=ref)
    assert res["tokens"] == ["F635G"]
    assert res["panel_hits"] == []


def test_substitutions_across_windows_sort_by_position(ref):
    mut = _mutate_codon(ref[0], 1354, "AGT")            # R1354S (HS2)
    mut = _mutate_codon(mut, 639, "TTT")                # S639F (HS1)
    res = F.call_windows(_windows(mut), reference=ref)
    assert res["tokens"] == ["S639F", "R1354S"]         # position order across windows
    assert res["panel_hits"] == ["S639F"]


# --- honesty: uncalled != wild type -----------------------------------------

def test_ambiguous_codon_at_hotspot_is_uncalled_not_wild_type(ref):
    mut = _mutate_codon(ref[0], 639, "TNT")             # N -> cannot be translated
    res = F.call_windows(_windows(mut), reference=ref)
    assert res["tokens"] == []                          # NOT called wild type
    assert 639 in res["windows"]["HS1"]["uncalled_positions"]
    assert res["uncalled_panel"] == [639]
    call, source = F.format_call(res)
    assert call == ""
    assert source == "sra-fks1-recaller:HS1=partial(uncalled:639),HS2=wild-type"


def test_gap_base_is_uncalled(ref):
    mut = _mutate_codon(ref[0], 640, "A-G")
    res = F.call_windows(_windows(mut), reference=ref)
    assert 640 in res["windows"]["HS1"]["uncalled_positions"]


# --- honesty: frame is asserted per window, not guessed ---------------------

def test_frameshift_in_one_window_refuses_only_that_window(ref):
    # Drop a base from HS1 (out of frame) but leave HS2 wild-type. HS1 must be refused with
    # a reason; HS2 must still be called wild-type -- one bad window never discards the other.
    wins = _windows(ref[0])
    wins["HS1"] = wins["HS1"][:-1]
    res = F.call_windows(wins, reference=ref)
    assert res["windows"]["HS1"]["status"] == "refused"
    assert "multiple of 3" in res["windows"]["HS1"]["refused"]
    assert res["windows"]["HS2"]["status"] == "wild-type"
    assert res["tokens"] == []
    _call, source = F.format_call(res)
    assert source.startswith("sra-fks1-recaller:HS1=refused(")


def test_in_frame_length_mismatch_is_refused_as_indel(ref):
    # A whole extra codon: multiple of 3, but longer than the window -> an in-window indel,
    # refused rather than mis-aligned.
    wins = _windows(ref[0])
    wins["HS1"] = wins["HS1"] + "GGG"
    res = F.call_windows(wins, reference=ref)
    assert res["windows"]["HS1"]["status"] == "refused"
    assert "indel in window" in res["windows"]["HS1"]["refused"]


def test_empty_window_is_refused(ref):
    wins = _windows(ref[0])
    wins["HS1"] = ""
    res = F.call_windows(wins, reference=ref)
    assert res["windows"]["HS1"]["status"] == "refused"
    assert "empty" in res["windows"]["HS1"]["refused"]


def test_missing_window_is_marked_missing_not_wild_type(ref):
    wins = _windows(ref[0])
    del wins["HS2"]
    res = F.call_windows(wins, reference=ref)
    assert res["windows"]["HS2"]["status"] == "missing"
    assert res["windows"]["HS2"]["tokens"] == []
    _call, source = F.format_call(res)
    assert "HS2=missing" in source


# --- slice bridge -----------------------------------------------------------

def test_slice_windows_wrong_length_is_refused(ref):
    with pytest.raises(F.WindowError, match="cannot anchor"):
        F.slice_windows_from_cds(ref[0][:-3])


def test_slice_windows_have_expected_lengths(ref):
    wins = _windows(ref[0])
    assert len(wins["HS1"]) == F.FKS1_WINDOWS["HS1"].nt_len == 27
    assert len(wins["HS2"]) == F.FKS1_WINDOWS["HS2"].nt_len == 27


# --- is_panel_token ---------------------------------------------------------

def test_is_panel_token_tags_the_known_panel():
    assert F.is_panel_token("S639F")
    assert F.is_panel_token("S639P")
    assert F.is_panel_token("S639Y")


def test_is_panel_token_rejects_non_panel_and_junk():
    assert not F.is_panel_token("S639T")   # panel residue, non-panel mutant letter
    assert not F.is_panel_token("R1354S")  # HS2 reported but not yet panel-tagged
    assert not F.is_panel_token("F635G")   # non-panel residue
    assert not F.is_panel_token("S639")    # not a substitution token
    assert not F.is_panel_token("")
    assert not F.is_panel_token(None)


def test_is_panel_token_matches_call_windows_panel_hits(ref):
    mut = _mutate_codon(ref[0], 639, "TTT")   # S639F
    res = F.call_windows(_windows(mut), reference=ref)
    from_string = {t for t in res["tokens"] if F.is_panel_token(t)}
    assert from_string == set(res["panel_hits"])


# --- downstream contract ----------------------------------------------------

def test_emitted_call_round_trips_through_emergence_parse(ref):
    mut = _mutate_codon(ref[0], 1354, "AGT")   # R1354S
    mut = _mutate_codon(mut, 639, "TTT")       # S639F
    call, _source = F.format_call(F.call_windows(_windows(mut), reference=ref))
    parsed = E.parse_call(call)
    assert parsed == {"S639F", "R1354S"}
    assert all(E.is_substitution(t) for t in parsed)
