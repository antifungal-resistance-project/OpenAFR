"""Offline tests for the ERG11 re-caller's pre-run guards (scripts/recall_erg11.py).

Both guards exist to be trustworthy BEFORE the expensive part of a run:

  * the samtools version gate refuses an install whose `samtools consensus` (the pipeline's
    final step, >=1.13) is missing -- so that failure is instant, not discovered after a
    multi-GB read download, and
  * the `plan` subcommand reports exactly what a `fill` would fetch/re-call, using only the
    snapshot on disk -- so it runs anywhere (incl. hosts without the bioinformatics tools).

Neither guard touches the network or an external tool, so both are pinned here with the
tool probes faked -- the same offline discipline as test_recall_orchestration.py.
"""
import argparse
import importlib.util
import pathlib

import pytest

from openafr import earlywarning as ew

_SCRIPT = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "recall_erg11.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("recall_erg11", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def script():
    return _load_script()


# --- samtools version gate ---------------------------------------------------

def _fake_probe(version_text, help_text):
    """A _run_text stand-in that answers --version / --help and nothing else."""
    def _run_text(cmd):
        if cmd[:2] == ["samtools", "--version"]:
            return version_text
        if cmd[:2] == ["samtools", "--help"]:
            return help_text
        return ""
    return _run_text


# `samtools --help` lists (or omits) the consensus subcommand depending on version.
_HELP_WITH = "Commands:\n  sort\n  index\n  consensus   produce a consensus\n"
_HELP_WITHOUT = "Commands:\n  sort\n  index\n  mpileup   pileup\n"


def test_version_gate_passes_on_modern_samtools(script, monkeypatch):
    monkeypatch.setattr(script, "_run_text",
                        _fake_probe("samtools 1.17\nUsing htslib 1.17\n", _HELP_WITH))
    script._require_samtools_consensus()   # must not raise / exit


def test_version_gate_rejects_old_samtools(script, monkeypatch):
    # 1.10 predates `samtools consensus`: help omits it AND the version is below the floor.
    monkeypatch.setattr(script, "_run_text",
                        _fake_probe("samtools 1.10\nUsing htslib 1.10\n", _HELP_WITHOUT))
    with pytest.raises(SystemExit) as ei:
        script._require_samtools_consensus()
    msg = str(ei.value)
    assert "1.10" in msg and "consensus" in msg and "1.13" in msg


def test_version_gate_rejects_when_consensus_absent_even_if_version_ok(script, monkeypatch):
    # Belt-and-suspenders: a build that parses as new but does not list consensus is refused
    # on the capability, not the number -- the help listing is the ground truth.
    monkeypatch.setattr(script, "_run_text",
                        _fake_probe("samtools 1.15\n", _HELP_WITHOUT))
    with pytest.raises(SystemExit):
        script._require_samtools_consensus()


def test_version_gate_trusts_help_when_version_unparseable(script, monkeypatch):
    # If --version can't be parsed but consensus is clearly listed, don't false-block a
    # working install.
    monkeypatch.setattr(script, "_run_text",
                        _fake_probe("samtools (unknown build)\n", _HELP_WITH))
    script._require_samtools_consensus()   # must not raise


def test_require_tools_chains_the_version_gate(script, monkeypatch):
    # Presence check passes (everything "on PATH"), so _require_tools must then reach the
    # version gate and fail there -- proving the two checks are wired together.
    monkeypatch.setattr(script.shutil, "which", lambda t: "/usr/bin/" + t)
    monkeypatch.setattr(script, "_run_text",
                        _fake_probe("samtools 1.9\n", _HELP_WITHOUT))
    with pytest.raises(SystemExit) as ei:
        script._require_tools()
    assert "consensus" in str(ei.value)


# --- plan preflight ----------------------------------------------------------

def _snapshot(path, rows):
    records = []
    for r in rows:
        rec = {c: "" for c in ew.SNAPSHOT_COLUMNS}
        rec.update(r)
        records.append(rec)
    ew.write_snapshot(str(path), records)


def _plan_args(snapshot, limit=0, list_runs=False, random=False, seed=0):
    return argparse.Namespace(snapshot=str(snapshot), limit=limit, list_runs=list_runs,
                              random=random, seed=seed)


_ROWS = [
    {"isolate_key": "K1", "run_acc": "SRR1", "resistance_source": "pending:sra-recaller"},
    # first linked run is what fill picks; plan must count it the same way
    {"isolate_key": "K2", "run_acc": "SRR2,SRR2b", "resistance_source": "pending:sra-recaller"},
    # pending but no run_acc -> not attemptable
    {"isolate_key": "K3", "run_acc": "", "resistance_source": "pending:sra-recaller"},
    # already resolved -> not pending
    {"isolate_key": "K4", "run_acc": "SRR4", "erg11_call": "Y132F",
     "resistance_source": "sra-recaller:called"},
]


def test_plan_counts_match_fill_selection(script, tmp_path, capsys):
    snap = tmp_path / "snap.tsv"
    _snapshot(snap, _ROWS)

    # plan must not need the tools -- calling it without faking the gate proves that.
    rc = script.cmd_plan(_plan_args(snap))
    assert rc == 0
    out = capsys.readouterr().out
    assert "total isolates:      4" in out
    assert "pending:sra-recaller 3" in out      # K1,K2,K3
    assert "1 with no run_acc" in out            # K3
    assert "re-call:  2 isolate(s)" in out       # K1,K2 (K3 has no run, K4 resolved)

    # and it must not mutate the snapshot
    keys = [r["isolate_key"] for r in ew.read_snapshot(str(snap))]
    assert keys == ["K1", "K2", "K3", "K4"]


def test_plan_honors_limit(script, tmp_path, capsys):
    snap = tmp_path / "snap.tsv"
    _snapshot(snap, _ROWS)
    script.cmd_plan(_plan_args(snap, limit=1))
    out = capsys.readouterr().out
    assert "re-call:  1 isolate(s)" in out and "--limit 1" in out


def test_plan_list_runs_prints_first_linked_accession(script, tmp_path, capsys):
    snap = tmp_path / "snap.tsv"
    _snapshot(snap, _ROWS)
    script.cmd_plan(_plan_args(snap, list_runs=True))
    out = capsys.readouterr().out
    # the first linked run of each attemptable isolate, not the secondary SRR2b
    assert "SRR1" in out and "SRR2" in out
    assert "SRR2b" not in out                    # secondary run is not what fill downloads
    assert "SRR4" not in out                     # resolved row is not planned


# --- randomized sampling (--random / --seed) ---------------------------------
#
# Snapshot order tracks NCBI submission order, which clusters by study, so the head of
# the pending list is a biased sample. `_fill_order(shuffle=True)` draws a seeded random
# --limit slice instead; plan and fill share it so the preflight matches the real run.

def _pending_records(n):
    """n attemptable pending rows, K00.. / SRR00.., in snapshot order."""
    return [{c: "" for c in ew.SNAPSHOT_COLUMNS} | {
        "isolate_key": f"K{i:02d}", "run_acc": f"SRR{i:02d}",
        "resistance_source": "pending:sra-recaller"} for i in range(n)]


def _keys(rows):
    return [r["isolate_key"] for r in rows]


def test_fill_order_no_shuffle_is_snapshot_head(script):
    recs = _pending_records(50)
    assert _keys(script._fill_order(recs, limit=5)) == ["K00", "K01", "K02", "K03", "K04"]


def test_fill_order_shuffle_is_deterministic_for_a_seed(script):
    recs = _pending_records(50)
    a = _keys(script._fill_order(recs, limit=10, shuffle=True, seed=0))
    b = _keys(script._fill_order(recs, limit=10, shuffle=True, seed=0))
    assert a == b                                   # same snapshot + seed -> same sample
    assert len(a) == 10 and len(set(a)) == 10       # a real 10-isolate sample, no dupes


def test_fill_order_shuffle_differs_from_head_and_across_seeds(script):
    recs = _pending_records(50)
    head = _keys(script._fill_order(recs, limit=10))
    s0 = _keys(script._fill_order(recs, limit=10, shuffle=True, seed=0))
    s1 = _keys(script._fill_order(recs, limit=10, shuffle=True, seed=1))
    assert s0 != head          # the whole point: not the study-clustered head
    assert s0 != s1            # a different seed draws a different sample
    # every sampled row is still a real pending isolate (a permutation, nothing invented)
    assert set(s0) <= set(_keys(recs))


def test_fill_order_shuffle_then_limit_slices_the_shuffled_pool(script):
    recs = _pending_records(50)
    full = _keys(script._fill_order(recs, shuffle=True, seed=7))
    sliced = _keys(script._fill_order(recs, limit=8, shuffle=True, seed=7))
    assert sliced == full[:8]              # limit is applied AFTER the shuffle, not before
    assert sorted(full) == sorted(_keys(recs))   # full shuffle is a permutation of the pool


def test_plan_random_previews_the_same_sample_fill_would_take(script, tmp_path, capsys):
    recs = _pending_records(50)
    snap = tmp_path / "snap.tsv"
    ew.write_snapshot(str(snap), recs)

    script.cmd_plan(_plan_args(snap, limit=6, random=True, seed=3))
    out = capsys.readouterr().out
    assert "random sample, seed=3" in out
    assert "re-call:  6 isolate(s)" in out

    # plan's listed accessions must be exactly what a fill --random --seed 3 --limit 6 hits
    script.cmd_plan(_plan_args(snap, limit=6, random=True, seed=3, list_runs=True))
    listed = {ln for ln in capsys.readouterr().out.splitlines() if ln.startswith("SRR")}
    expected = {r["run_acc"] for r in script._fill_order(recs, limit=6, shuffle=True, seed=3)}
    assert listed == expected
