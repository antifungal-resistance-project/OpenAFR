"""Offline tests for the temporal (prospective) holdout scaffold (issue #82, follow-up to #108).

The numeric run is release-gated — it needs a ChEMBL release newer than the 2026-08 tuning
corpus — so what has to be right NOW, and is provable offline, is:

  * the temporal split predicate: a molecule is in the holdout iff its earliest antifungal
    document year is STRICTLY after the frozen cutoff (an unknown year is not temporal);
  * the release gate itself: with no temporal set on disk, the grader must report
    `not-yet-run (release-gated)` and exit incomplete (2) — never a spurious pass or a crash;
  * the pre-registration is frozen (its pinned sha matches the committed file), so the rule,
    cutoff and interpretation cannot be edited after the fact;
  * the disjointness accounting includes the #108 external sets, so nothing already used can
    leak back into the temporal holdout.
"""
import hashlib
import pathlib

from scripts import make_temporal_holdout as mth
from scripts.validate_gate_external import classify_outcome  # reused verbatim by the temporal gate
from scripts.validate_gate_temporal import PREREG, PREREG_SHA

MIN_AUC = 0.70


def test_temporal_predicate_requires_strictly_after_cutoff():
    assert mth.is_temporal(mth.CUTOFF_YEAR + 1) is True
    assert mth.is_temporal(mth.CUTOFF_YEAR) is False        # same year overlaps the corpus
    assert mth.is_temporal(mth.CUTOFF_YEAR - 5) is False


def test_unknown_year_is_not_temporal():
    """Absence of a document year is not proof the molecule post-dates the corpus."""
    assert mth.is_temporal(None) is False


def test_earliest_document_year_takes_the_minimum():
    recs = [{"document_year": 2028}, {"document_year": 2027}, {"document_year": 2030}]
    assert mth.earliest_document_year(recs) == 2027


def test_earliest_document_year_ignores_missing_and_nonint_years():
    recs = [{"document_year": None}, {"foo": 1}, {"document_year": "2027"}, {"document_year": 2029}]
    assert mth.earliest_document_year(recs) == 2029
    assert mth.earliest_document_year([{"document_year": None}]) is None
    assert mth.earliest_document_year([]) is None


def test_a_recent_but_pre_cutoff_molecule_is_excluded():
    """End-to-end predicate: a 2026 deposition (the corpus's latest year) is NOT temporal."""
    assert mth.is_temporal(mth.earliest_document_year([{"document_year": 2026}])) is False
    assert mth.is_temporal(mth.earliest_document_year([{"document_year": 2027}])) is True


def test_cutoff_is_the_frozen_value():
    assert mth.CUTOFF_YEAR == 2026


def test_used_ids_includes_the_external_holdout_sets():
    """The temporal set must be disjoint from #108's external sets, so those ids must be counted
    among the 'already used' ids the builder excludes."""
    used = mth.used_ids()
    ext_actives = {mth._id_of(n) for _s, n in mth._read_smi("data/ligands/external_actives.smi")}
    ext_decoys = {mth._id_of(n) for _s, n in mth._read_smi("data/ligands/external_decoys.smi")}
    assert ext_actives and ext_actives <= used
    assert ext_decoys and ext_decoys <= used


def test_grader_is_release_gated_when_no_set_exists(monkeypatch, capsys):
    """No temporal set on disk today -> the grader reports not-yet-run and exits incomplete (2),
    the same honesty as the multi-baseline gate, rather than crashing or faking a pass."""
    assert not pathlib.Path(mth.ACTIVES_OUT).exists()      # nothing built yet (release-gated)
    monkeypatch.setattr("sys.argv", ["validate_gate_temporal.py"])
    from scripts.validate_gate_temporal import main
    assert main() == 2
    out = capsys.readouterr().out
    assert "not-yet-run (release-gated)" in out


def test_prereg_hash_matches_the_frozen_document():
    got = hashlib.sha256(pathlib.Path(PREREG).read_bytes()).hexdigest()
    assert got == PREREG_SHA


def test_three_outcome_rule_is_the_same_one_108_committed():
    """The temporal gate reuses #108's classify_outcome verbatim — pin that it still behaves."""
    assert classify_outcome(0.80, 0.71, MIN_AUC) == "PASS"
    assert classify_outcome(0.715, 0.630, MIN_AUC) == "SUPPORTED-BUT-WIDE"
    assert classify_outcome(0.69, 0.55, MIN_AUC) == "FAIL"
