"""Known-answer tests for the run-2 validation gate.

These lock the two functions the entire OpenAFR result rests on:

    metrics(ordered_flags) -> (AUC, [EF@1%, EF@5%, EF@10%], BEDROC)
        Pure scoring on a list of activity flags ALREADY sorted best-first.

    report(label, ranked)  -> (AUC, [EFs], BEDROC)
        Takes (name, key, is_active) rows, lower key = better, key None =
        unrankable. Sorts, appends unrankable LAST (never drops them), then
        calls metrics().

Every expected number below was hand-derived and cross-checked against RDKit
directly, so a future refactor that silently changes the math fails loudly.

    ranking                         AUC    EF@1/5/10%   note
    ------------------------------  -----  -----------  ---------------------
    10 actives, then 90 decoys      1.0    10/10/10     perfect screen
    90 decoys, then 10 actives      0.0    0/0/0        worst possible
    [active, decoy, active, decoy]  0.75   -            3 of 4 pairs ordered
    [1,0,0,1] (unrankable last)     0.5    -            active dragged to end
    [1,0,0]   (same, but DROPPED)   1.0    -            the inflation we forbid
"""
import pytest

import validate_gate2 as g


# --------------------------------------------------------------------------
# metrics(): pure scoring on a pre-ordered flag list
# --------------------------------------------------------------------------

def test_perfect_ranking_scores_max():
    # 10 actives ahead of 90 decoys. Textbook perfect screen.
    flags = [1] * 10 + [0] * 90
    auc, efs, bedroc = g.metrics(flags)
    assert auc == 1.0
    assert efs == [10.0, 10.0, 10.0]      # EF maxes at 1/prevalence = 1/0.1
    assert bedroc == pytest.approx(1.0)


def test_reversed_ranking_scores_zero():
    # Every active ranked behind every decoy. Worst possible.
    flags = [0] * 90 + [1] * 10
    auc, efs, bedroc = g.metrics(flags)
    assert auc == 0.0
    assert efs == [0.0, 0.0, 0.0]
    assert bedroc == pytest.approx(0.0)


def test_interleaved_ranking_known_auc():
    # [active, decoy, active, decoy]: of the 4 (active,decoy) pairs, 3 are
    # correctly ordered (only the 2nd active loses to the 1st decoy) -> 0.75.
    auc, _, _ = g.metrics([1, 0, 1, 0])
    assert auc == 0.75


def test_auc_is_symmetric_about_half():
    # Reversing a ranking must reflect AUC across 0.5. Guards the orientation
    # convention (best-first) from silently flipping in a refactor.
    flags = [1, 1, 0, 1, 0, 0, 0, 1, 0, 0]
    auc_fwd, _, _ = g.metrics(flags)
    auc_rev, _, _ = g.metrics(list(reversed(flags)))
    assert auc_fwd + auc_rev == pytest.approx(1.0)


# --------------------------------------------------------------------------
# report(): sorting + the "unrankable ranked LAST, never dropped" rule
# --------------------------------------------------------------------------

def test_report_sorts_by_key_ascending():
    # Deliberately unsorted input; best (lowest key) should win.
    ranked = [
        ("d1", 3.0, False),
        ("a1", 1.0, True),
        ("d2", 2.0, False),
    ]
    auc, _, _ = g.report("t", ranked)
    # ordered -> a1(1.0), d2(2.0), d1(3.0) -> flags [1,0,0] -> perfect
    assert auc == 1.0


def test_unrankable_active_ranked_last_not_dropped():
    # a2 has no usable pose (key None). The pre-registration rule says rank it
    # LAST, never drop it. That penalty is the whole point of the gate.
    ranked = [
        ("a1", 1.0, True),
        ("d1", 2.0, False),
        ("d2", 3.0, False),
        ("a2", None, True),   # unrankable active
    ]
    auc, _, _ = g.report("t", ranked)
    # ordered -> a1, d1, d2, a2 -> flags [1,0,0,1] -> AUC 0.5
    assert auc == 0.5
    # The contrast that proves the rule matters: had a2 been DROPPED, the
    # remaining [1,0,0] would score a perfect 1.0. Dropping unrankable actives
    # is the single easiest way to fake a passing gate.
    assert g.metrics([1, 0, 0])[0] == 1.0


def test_unrankable_decoy_also_ranked_last():
    # Symmetry: an unrankable decoy is also placed last (a[1] is None -> bad).
    ranked = [
        ("a1", 1.0, True),
        ("a2", 2.0, True),
        ("d1", None, False),
    ]
    auc, _, _ = g.report("t", ranked)
    # ordered -> a1, a2, d1 -> flags [1,1,0] -> perfect
    assert auc == 1.0


def test_report_is_input_order_independent():
    # Shuffling the input rows must not change the result: report() sorts
    # internally, so the caller's ordering is irrelevant (except exact ties,
    # see below).
    rows = [
        ("a1", 1.0, True),
        ("d1", 2.0, False),
        ("d2", 3.0, False),
        ("a2", 4.0, True),
    ]
    baseline = g.report("t", rows)
    shuffled = g.report("t", list(reversed(rows)))
    assert baseline == shuffled


def test_exact_tie_is_broken_by_input_order():
    # KNOWN LIMITATION (eng-review finding Q3): equal keys are broken by input
    # order because Python's sort is stable. With real float distances this
    # essentially never happens, but it means the metric is not fully
    # deterministic across filesystems (os.listdir order). A deterministic
    # secondary sort key (e.g. name) would remove this. This test pins the
    # CURRENT behavior so any future change is a conscious one.
    active_first = [("a", 1.0, True), ("d", 1.0, False)]
    decoy_first = [("d", 1.0, False), ("a", 1.0, True)]
    assert g.report("t", active_first)[0] == 1.0   # [1,0]
    assert g.report("t", decoy_first)[0] == 0.0    # [0,1]


def test_report_return_matches_metrics_on_reconstructed_order():
    # report() must return exactly metrics() applied to the order it built.
    # Internal-consistency guard: sorting and scoring can't drift apart.
    ranked = [
        ("d1", 2.0, False),
        ("a1", 1.0, True),
        ("a2", None, True),
        ("d2", 3.0, False),
    ]
    got = g.report("t", ranked)
    # hand-reconstruct: ok sorted [a1(1), d1(2), d2(3)] + bad [a2] -> [1,0,0,1]
    expected = g.metrics([1, 0, 0, 1])
    assert got == expected
