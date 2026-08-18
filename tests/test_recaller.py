"""Tests for the ERG11 re-caller deterministic core (openafr/recaller.py).

These lock the honesty guarantees the whole early-warning claim rests on -- every
downstream stage fires off the tokens this core emits, so a silent mis-call here is a
silent mis-call in the alert a clinician reads:
  * the pinned reference translates to the panel's wild-type hotspot residues
    (the numbering check that makes an emitted 'Y132F' mean residue 132),
  * a wild-type consensus emits NO call (never invents a substitution) and is
    distinguishable from an uncalled one,
  * a real substitution codon emits the right token and is panel-tagged,
  * an ambiguous codon at a hotspot is UNCALLED, never silently "wild type",
  * a frameshift/indel or length mismatch is REFUSED, never frame-guessed,
  * format_call maps results to the (erg11_call, resistance_source) snapshot pair
    in the exact token shape emergence.parse_call consumes.
Runs against the pinned reference on disk; no network, no external tools.
"""
import importlib.util
import pathlib

import pytest

from openafr import recaller as R
from openafr import emergence as E

_SCRIPT = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "recall_erg11.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("recall_erg11", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def ref():
    """The pinned (cds, protein) reference, loaded once."""
    return R.load_reference()


def _mutate_codon(cds, pos1, new_codon):
    """Replace the codon for 1-based residue `pos1` with `new_codon`."""
    i = (pos1 - 1) * 3
    assert len(new_codon) == 3
    return cds[:i] + new_codon + cds[i + 3:]


# --- reference integrity / numbering ---------------------------------------

def test_reference_loads_and_hashes(ref):
    cds, protein = ref
    assert len(cds) == 1575
    assert len(protein) == 524          # stop stripped
    assert protein.startswith("MALKDCIVDVVDRF")


def test_reference_numbering_matches_panel(ref):
    _cds, protein = ref
    # The whole point of the reference: these are the wild-type letters the panel
    # names its mutations *from*. If any drift, every token would be mis-numbered.
    assert (protein[124], protein[125], protein[131], protein[142]) == ("V", "F", "Y", "K")


def test_reference_hash_mismatch_is_refused(tmp_path):
    bad = tmp_path / "erg11_cds.fasta"
    # A valid-looking but wrong-length/content CDS: 6 in-frame bases.
    bad.write_text(">bad\nATGAAA\n")
    with pytest.raises(R.ReferenceError):
        R.load_reference(str(bad))


# --- wild type: no call, and distinguishable from uncalled ------------------

def test_wild_type_emits_no_call(ref):
    cds, _ = ref
    res = R.call_substitutions(cds, reference=ref)
    assert res["tokens"] == []
    assert res["uncalled_positions"] == []
    assert res["n_covered"] == res["n_residues"] == 524
    call, source = R.format_call(res)
    assert call == ""
    assert source == "sra-recaller:wild-type"


# --- a real substitution ----------------------------------------------------

def test_y132f_is_called_and_panel_tagged(ref):
    cds, _ = ref
    mut = _mutate_codon(cds, 132, "TTT")     # Tyr (TAC) -> Phe (TTT)
    res = R.call_substitutions(mut, reference=ref)
    assert res["tokens"] == ["Y132F"]
    assert res["panel_hits"] == ["Y132F"]
    call, source = R.format_call(res)
    assert call == "Y132F"
    assert source == "sra-recaller:called"


def test_multiple_substitutions_sort_by_position(ref):
    cds, _ = ref
    mut = _mutate_codon(cds, 143, "AGG")     # Lys (AAG/AAA) -> Arg (AGG)
    mut = _mutate_codon(mut, 132, "TTT")     # + Y132F
    res = R.call_substitutions(mut, reference=ref)
    assert res["tokens"] == ["Y132F", "K143R"]     # position order, not input order
    assert res["panel_hits"] == ["Y132F", "K143R"]


def test_non_panel_substitution_is_emitted_but_not_tagged(ref):
    cds, protein = ref
    # Pick a residue that is NOT a panel position and force a change. Residue 1 is Met;
    # change codon 2 to something different from wild type without touching a hotspot.
    wt2 = protein[1]
    # Ala codon GCT -> whatever residue 2 is, force a distinct amino acid.
    mut = _mutate_codon(cds, 2, "GGT")       # -> Gly
    res = R.call_substitutions(mut, reference=ref)
    assert res["tokens"] == [f"{wt2}2G"]
    assert res["panel_hits"] == []           # emitted verbatim, never asserted resistance


def test_stop_codon_substitution_is_called(ref):
    cds, protein = ref
    mut = _mutate_codon(cds, 132, "TAA")     # Tyr -> stop
    res = R.call_substitutions(mut, reference=ref)
    assert res["tokens"] == ["Y132*"]


# --- honesty: uncalled != wild type -----------------------------------------

def test_ambiguous_codon_at_hotspot_is_uncalled_not_wild_type(ref):
    cds, _ = ref
    mut = _mutate_codon(cds, 143, "ANG")     # N -> codon cannot be confidently translated
    res = R.call_substitutions(mut, reference=ref)
    assert res["tokens"] == []               # NOT called wild type
    assert 143 in res["uncalled_positions"]
    assert res["uncalled_panel"] == [143]
    assert res["n_covered"] == 523           # one residue lost to ambiguity
    call, source = R.format_call(res)
    assert call == ""
    assert source == "sra-recaller:partial(panel-uncalled:143)"


def test_gap_base_is_uncalled(ref):
    cds, _ = ref
    mut = _mutate_codon(cds, 100, "A-G")
    res = R.call_substitutions(mut, reference=ref)
    assert 100 in res["uncalled_positions"]


# --- honesty: frame is asserted, not guessed --------------------------------

def test_length_not_multiple_of_three_is_refused(ref):
    cds, _ = ref
    with pytest.raises(R.ConsensusError):
        R.call_substitutions(cds[:-1], reference=ref)


def test_wrong_length_multiple_of_three_is_refused(ref):
    cds, _ = ref
    with pytest.raises(R.ConsensusError):
        # -6 drops the stop codon AND one real coding codon -> 523 aa vs 524 (a
        # -3 that only drops the stop codon is benign and correctly NOT refused).
        R.call_substitutions(cds[:-6], reference=ref)


def test_empty_consensus_is_refused(ref):
    with pytest.raises(R.ConsensusError):
        R.call_substitutions("   ", reference=ref)


# --- the downstream contract ------------------------------------------------

# --- orchestration honesty: a failed isolate is MARKED, never dropped/guessed -----

def test_recall_one_marks_tool_failure_without_guessing(ref, monkeypatch):
    """If reads->consensus blows up, the isolate must come back with an empty call and
    a 'failed' source -- not a silent wild-type and not an exception that drops it."""
    script = _load_script()

    def boom(*a, **k):
        raise RuntimeError("fasterq-dump: connection reset")

    monkeypatch.setattr(script, "reads_to_consensus", boom)
    call, source = script._recall_one("SRR000001", ref, R.DEFAULT_REFERENCE)
    assert call == ""
    assert source.startswith("sra-recaller:failed(")
    assert "wild-type" not in source          # a failure is never sequenced-wild-type


def test_recall_one_marks_frame_refusal(ref, monkeypatch):
    """A returned consensus that is out of frame is refused, marked, and does not raise."""
    script = _load_script()
    cds, _ = ref

    def short_consensus(run_acc, workdir, reference_fasta, min_depth=10):
        p = pathlib.Path(workdir) / "c.fasta"
        p.write_text(">c\n" + cds[:-6] + "\n")   # a codon short -> refused
        return str(p)

    monkeypatch.setattr(script, "reads_to_consensus", short_consensus)
    call, source = script._recall_one("SRR000002", ref, R.DEFAULT_REFERENCE)
    assert call == ""
    assert source.startswith("sra-recaller:refused(")


def test_emitted_call_round_trips_through_emergence_parse(ref):
    cds, _ = ref
    mut = _mutate_codon(cds, 143, "AGG")
    mut = _mutate_codon(mut, 132, "TTT")
    call, _source = R.format_call(R.call_substitutions(mut, reference=ref))
    # The whole point of the token format: emergence.parse_call must recover the tokens
    # the re-caller emitted, and each must be a canonical substitution.
    parsed = E.parse_call(call)
    assert parsed == {"Y132F", "K143R"}
    assert all(E.is_substitution(t) for t in parsed)


# ---- is_panel_token: rebuild panel membership from a stored token string ----

def test_is_panel_token_tags_the_known_panel():
    # Every panel mutant letter at its residue is a hit; the tokens the caller tags.
    assert R.is_panel_token("V125A")
    assert R.is_panel_token("F126L")
    assert R.is_panel_token("Y132F")
    assert R.is_panel_token("K143R")


def test_is_panel_token_rejects_non_panel_and_junk():
    assert not R.is_panel_token("Y132H")   # panel residue, non-panel mutant letter
    assert not R.is_panel_token("A61T")    # non-panel residue
    assert not R.is_panel_token("K143")    # not a substitution token
    assert not R.is_panel_token("")
    assert not R.is_panel_token(None)


def test_is_panel_token_matches_call_substitutions_panel_hits(ref):
    # The rebuilt-from-string test must agree with the caller's own panel_hits, so a
    # prevalence count never disagrees with what the caller tagged at call time.
    cds, _ = ref
    mut = _mutate_codon(cds, 132, "TTT")   # Y132F
    mut = _mutate_codon(mut, 143, "AGG")   # K143R
    result = R.call_substitutions(mut, reference=ref)
    from_string = {t for t in result["tokens"] if R.is_panel_token(t)}
    assert from_string == set(result["panel_hits"])
