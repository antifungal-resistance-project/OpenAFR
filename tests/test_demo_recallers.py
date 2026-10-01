"""Tests for the one-command re-caller demo (scripts/demo_recallers.py, issue C1 #162).

The demo is the external on-ramp: a tiny committed sample dataset the deterministic `call`
path runs on anywhere, with a committed EXPECTED.txt oracle. Two things are pinned here:

  * on the committed samples, every run matches EXPECTED.txt (the on-ramp actually works),
  * a wrong expected line makes --check FAIL (a drift-detector that cannot detect drift is
    worse than none) -- mirroring the tamper cases in test_repro.py.

The tamper case uses a temp EXPECTED fixture so it never touches the committed oracle.
"""
import pathlib

from scripts import demo_recallers as demo


def test_every_sample_matches_expected():
    """The four committed samples reproduce their committed call/source/panel lines."""
    res = demo.results()
    expected = demo.parse_expected(demo.EXPECTED_PATH.read_text())
    assert set(res) == set(expected), "sample set drifted from EXPECTED.txt"
    for label in res:
        assert res[label] == expected[label], (
            f"{label}: {res[label]} != expected {expected[label]}"
        )


def test_check_passes_on_committed_tree(monkeypatch):
    """`demo_recallers.py --check` exits 0 on the committed samples + EXPECTED.txt."""
    monkeypatch.setattr("sys.argv", ["demo_recallers.py", "--check"])
    assert demo.main() == 0


def test_check_detects_a_tampered_expected(tmp_path, monkeypatch):
    """A wrong expected line is caught: --check returns non-zero."""
    text = demo.EXPECTED_PATH.read_text().replace("Y132F", "Y132H")
    bad = tmp_path / "EXPECTED.txt"
    bad.write_text(text)
    monkeypatch.setattr(demo, "EXPECTED_PATH", bad)
    monkeypatch.setattr("sys.argv", ["demo_recallers.py", "--check"])
    assert demo.main() == 1


def test_samples_derive_from_the_pinned_references():
    """Each wild-type sample is the pinned reference CDS verbatim (no silent re-pin)."""
    root = pathlib.Path(demo.ROOT)

    def seq(path):
        lines = path.read_text().splitlines()
        return "".join(l for l in lines[1:] if l.strip())

    erg_ref = seq(root / "data/earlywarning/erg11_reference/erg11_cds.fasta")
    fks_ref = seq(root / "data/earlywarning/fks1_reference/fks1_cds.fasta")
    assert seq(demo.DEMO_DIR / "erg11_wildtype.fasta") == erg_ref
    assert seq(demo.DEMO_DIR / "fks1_wildtype.fasta") == fks_ref

    # The mutants differ from the reference by exactly one nucleotide (the single codon edit).
    erg_mut = seq(demo.DEMO_DIR / "erg11_Y132F.fasta")
    fks_mut = seq(demo.DEMO_DIR / "fks1_S639F.fasta")
    assert len(erg_mut) == len(erg_ref)
    assert len(fks_mut) == len(fks_ref)
    assert sum(a != b for a, b in zip(erg_mut, erg_ref)) == 1
    assert sum(a != b for a, b in zip(fks_mut, fks_ref)) == 1
