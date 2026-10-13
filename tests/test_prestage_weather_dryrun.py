"""Guard the FKS1 fill -> weather-page render path at scale (scripts/prestage_weather_dryrun.py).

The harness manufactures a release-shaped synthetic world (a raw-NCBI --tsv of N pending
isolates + a matching committed call table) and pushes it through the REAL render_weather CLI.
These tests run it small-but-wide so CI keeps the schema/scale contract green:

  * the overlay matches every isolate by isolate_key (schema contract -- #171 step 1),
  * the page leaves day-0 with an honest coverage clause (render contract -- #171 step 5),
  * coverage counts RESOLVED isolates (wild-type included), not just mutation-carriers, and
    prevalence is over that honest denominator -- the regression the dry-run first caught.

No network, no cloud: everything lives in a tmp dir the harness cleans up itself.
"""
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
import prestage_weather_dryrun as dry  # noqa: E402

AS_OF = "2026-01-01"       # a fixed baseline cut so the recent window is deterministic
WINDOW = 180
N = 300                    # small-but-wide: enough to exercise overlay + coverage at scale


def _render(tmp_path, n=N, resistant_frac=0.06, release=710):
    """Drive the harness end-to-end into tmp_path; return (proc, page_text, n_resistant)."""
    calls_dir = tmp_path / "fks1_calls"
    calls_dir.mkdir()
    tsv = tmp_path / "pull.tsv"
    table = calls_dir / f"{dry.ew.ORG_PREFIX}.{release}.tsv"
    out = tmp_path / "index.html"

    raw_rows, call_records = dry.synthesize(n, resistant_frac, AS_OF, WINDOW)
    dry._write_raw_tsv(tsv, raw_rows)
    sha, n_called = dry.ew.write_fks1_call_table(str(table), call_records)
    dry.ew.update_fks1_calls_index(str(calls_dir / "INDEX.tsv"), {
        "release_tag": f"{dry.ew.ORG_PREFIX}.{release}", "as_of": AS_OF,
        "window_days": WINDOW, "n_called": n_called, "sha256": sha,
        "written_at_utc": "2026-01-02T00:00:00Z",
    })
    proc, _elapsed = dry._run_render(tsv, table, AS_OF, WINDOW, out)
    return proc, (out.read_text() if out.exists() else ""), round(n * resistant_frac)


def test_render_never_crashes_and_overlays_every_isolate(tmp_path):
    proc, page, _ = _render(tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert page
    m = re.search(r"overlaid (\d+) FKS1 call", proc.stderr)
    assert m and int(m.group(1)) == N  # every synthetic isolate matched by isolate_key


def test_page_leaves_day0_with_full_window_coverage(tmp_path):
    # The whole window is re-called, so the page must show 100% coverage -- never the day-0
    # watching copy, and never an understated % that counts only mutation-carriers.
    _, page, _ = _render(tmp_path)
    assert f"{N} of {N} window isolates re-called (100%)" in page
    assert "calls current through" in page and AS_OF in page


def test_all_wild_type_window_is_a_clean_all_clear_not_day0(tmp_path):
    # A fully re-called, mutation-free window is a real all-clear: full coverage, no false
    # "awaiting re-call" / day-0 state just because nothing was positive.
    _, page, _ = _render(tmp_path, resistant_frac=0.0)
    assert f"{N} of {N} window isolates re-called (100%)" in page
    assert "awaiting re-call" not in page


def test_prevalence_is_over_the_resolved_denominator_not_carriers(tmp_path):
    # The regression the dry-run caught: with ~6% of a fully re-called window carrying a
    # hotspot call, each variant's recent prevalence must read ~its true share of ALL
    # re-called isolates (low single digits), never ~its share of carriers-only (~25%).
    _, page, n_resistant = _render(tmp_path, resistant_frac=0.08)
    shares = [float(x) for x in re.findall(r"(\d+\.\d+)%", page)]
    assert shares, "expected at least one prevalence percentage on the page"
    # No single hotspot can be more prevalent than the whole resistant fraction of the window.
    ceiling = 100.0 * n_resistant / N + 0.5
    assert max(shares) <= ceiling, (max(shares), ceiling)
