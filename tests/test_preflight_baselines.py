"""Tests for the offline multi-baseline preflight (scripts/preflight_baselines.py, #83).

The preflight is the go/no-go a cloud operator runs before spending time on the metal-aware
rescores. These lock its two load-bearing promises: it verifies the frozen surface (the
pre-registration + protocol.yaml hashes), and it reports the honest run-state without ever
invoking a scorer. Stdlib-only, so it runs on the same osx-arm64 host that cannot run the
baselines -- which is the whole point.
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
import preflight_baselines as pf  # noqa: E402


def test_frozen_hashes_currently_verify():
    # Drift guard: the pinned pre-registration + protocol must match their frozen SHAs.
    ok, lines = pf._check_frozen_sha(pf.PREREG_SHA_FILE)
    assert ok, "\n".join(lines)
    proto = pf._sha256(pf.ROOT / "protocol.yaml")
    from openafr import protocol
    assert proto == protocol.PROTOCOL_SHA


def test_sha256_matches_the_pinned_value():
    line = (pf.ROOT / pf.PREREG_SHA_FILE).read_text().split()
    want = line[0]
    assert pf._sha256(pf.ROOT / pf.PREREG) == want


def test_pose_count_none_when_dir_absent_int_when_present():
    assert pf._pose_count("work/screen2") is None            # gitignored / not regenerated
    assert pf._pose_count("work") == 0                        # a real dir, no .pdbqt at top


def test_check_frozen_sha_flags_a_mismatch(tmp_path, monkeypatch):
    # Build a throwaway tree with a file and a sha line carrying the WRONG hash, point
    # ROOT at it, and confirm the real helper reports not-ok (a DRIFT).
    (tmp_path / "payload.txt").write_text("real bytes\n")
    (tmp_path / "bad.sha256").write_text("0" * 64 + "  payload.txt\n")
    monkeypatch.setattr(pf, "ROOT", tmp_path)
    ok, lines = pf._check_frozen_sha("bad.sha256")
    assert ok is False
    assert any(pf.BAD in ln for ln in lines)


def test_check_frozen_sha_flags_a_missing_pinned_file(tmp_path, monkeypatch):
    (tmp_path / "gone.sha256").write_text("0" * 64 + "  not_here.txt\n")
    monkeypatch.setattr(pf, "ROOT", tmp_path)
    ok, lines = pf._check_frozen_sha("gone.sha256")
    assert ok is False
    assert any(pf.MISSING in ln for ln in lines)


def test_main_reports_run_steps_remain_in_the_current_tree(monkeypatch, capsys):
    # No poses regenerated + no scorers present -> exit 2 (the normal pre-run state),
    # never 0 (a spurious ready) and never 1 (drift) on an intact tree.
    monkeypatch.setattr(sys, "argv", ["preflight_baselines.py"])
    rc = pf.main()
    out = capsys.readouterr().out
    assert rc == 2
    assert "RUN-STEPS-REMAIN" in out
    assert "not-yet-run" in out
