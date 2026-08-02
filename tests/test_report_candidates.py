"""Known-answer test for the candidate report (scripts/report_candidates.py).

The report is a presentation layer — it must not change the ranking and it must
not quietly drop the caveats that keep the list honest. So this locks: molecules
come out in the validated order (smallest N-Fe first, no-pose last), every
screened molecule appears exactly once, and the two non-negotiable framings
(hypothesis-not-a-hit, valid-only-after-the-gate-passes) are present. collect()
is exercised against real synthetic pose files so the parsing path is covered too.
"""
import pytest

from scripts.report_candidates import collect, render

# One receptor with a heme iron at the origin-ish; three ligands at known N-Fe
# distances so the expected order is unambiguous.
_RECEPTOR = "HETATM 9999 FE   HEM A 500      10.000  10.000  10.000  1.00  0.00          FE\n"


def _pose(score, n_xyz):
    """One-model pdbqt with a single nitrogen at n_xyz (or none if n_xyz is None)."""
    lines = ["MODEL 1", f"REMARK VINA RESULT:   {score:.1f}      0.000      0.000"]
    if n_xyz is not None:
        x, y, z = n_xyz
        lines.append(f"ATOM      1  N   LIG A   1    {x:8.3f}{y:8.3f}{z:8.3f}"
                     "  1.00  0.00           N")
    lines += ["ATOM      2  C   LIG A   1      50.000  50.000  50.000  1.00  0.00           C",
              "ENDMDL", ""]
    return "\n".join(lines)


@pytest.fixture
def screen(tmp_path):
    (tmp_path / "far.pdbqt").write_text(_pose(-7.0, (10.0, 10.0, 15.0)))    # 5.0 A
    (tmp_path / "close.pdbqt").write_text(_pose(-6.0, (10.0, 10.0, 12.0)))  # 2.0 A -> best
    (tmp_path / "nopose.pdbqt").write_text(_pose(-9.0, None))               # no N -> last
    rec = tmp_path / "rec.pdb"
    rec.write_text(_RECEPTOR)
    return str(tmp_path), str(rec)


def test_collect_orders_by_validated_criterion(screen):
    screen_dir, rec = screen
    fe, rows = collect(screen_dir, rec)
    assert fe == (10.0, 10.0, 10.0)
    # closest N-Fe first; the no-nitrogen molecule is last, not dropped.
    assert [r["name"] for r in rows] == ["close", "far", "nopose"]
    assert rows[0]["best_fe"] == pytest.approx(2.0)
    assert rows[-1]["best_fe"] is None


def test_best_score_does_not_decide_order(screen):
    # nopose has the strongest Vina score (-9.0) yet ranks last: score never decides.
    screen_dir, rec = screen
    _fe, rows = collect(screen_dir, rec)
    assert rows[-1]["name"] == "nopose"


def test_render_lists_every_molecule_once(screen):
    screen_dir, rec = screen
    fe, rows = collect(screen_dir, rec)
    md = render(screen_dir, rec, fe, rows, top=10)
    for name in ("close", "far", "nopose"):
        assert md.count(f"| {name} |") == 1
    assert "no pose" in md                       # unrankable shown, not hidden


def test_render_keeps_the_honesty_caveats(screen):
    screen_dir, rec = screen
    fe, rows = collect(screen_dir, rec)
    md = render(screen_dir, rec, fe, rows, top=10)
    assert "hypotheses, not hits" in md
    assert "validate_gate2.py" in md
    assert "applicability domain" in md.lower()


def test_render_reports_protocol_integrity(screen):
    # The frozen protocol is unmodified in the repo, so the report must say so.
    screen_dir, rec = screen
    fe, rows = collect(screen_dir, rec)
    md = render(screen_dir, rec, fe, rows, top=10)
    assert "unmodified since freezing" in md
