"""CI regression guard for the headline diagnostic-accuracy claim (issue #137, v3 PASS).

The echinocandin verdict arm was graded against phenotype on the 98-isolate PMC12323592 fixture
and, under the v3 high-PPV re-tier, cleared its pre-committed CLSI bar:
**VME 0/42 = 0.0%, ME 1/46 = 2.2%, 10 abstained -> verdict PASS** (work/RESULTS_diagnostic_accuracy.md).
That headline is load-bearing (docs/DIAGNOSTIC_GO_NO_GO.md, QUICKSTART_LAB.md, the JCM manuscript)
but -- unlike the caller's concordance headline (tests/test_fks1_concordance_replay.py) -- had no
test that re-derives it. This closes that asymmetry: it re-scores the frozen fixture through the
production metric + bar on every PR, so a regression in openafr.accuracy / openafr.verdict, or a
slip in the frozen pass bar, fails CI instead of silently invalidating the published claim.

SCOPE -- mirrors the concordance replay's discipline:
  * GUARDS the accuracy SCORING path (openafr.accuracy.evaluate + apply_bar), the frozen pass bar
    (BAR in scripts/validate_diagnostic_accuracy.py), and the frozen per-isolate verdicts the
    headline is built from (data/earlywarning/recaller_sanity/fks1_accuracy_98_v3.tsv, whose SHA
    is pinned via work/PREREG_diagnostic_accuracy_v3.sha256 and re-checked by scripts/repro.py).
  * Does NOT re-run the aligner/caller or re-derive the verdicts from reads -- that is the
    tool/data-gated grader (scripts/validate_diagnostic_accuracy.py). The v2->v3 re-tier
    derivation itself is covered by tests/test_accuracy_harvest_v3.py.
No network, no tools.
"""
import importlib.util
import pathlib

from openafr import accuracy as ac

ROOT = pathlib.Path(__file__).resolve().parent.parent
FIXTURE = ROOT / "data" / "earlywarning" / "recaller_sanity" / "fks1_accuracy_98_v3.tsv"
GRADER = ROOT / "scripts" / "validate_diagnostic_accuracy.py"


def _bar():
    """The frozen pass bar, read from the grader so the test and the grader cannot disagree."""
    spec = importlib.util.spec_from_file_location("validate_diagnostic_accuracy", GRADER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.BAR


def _records():
    """The frozen fixture as accuracy records (same mapping the grader uses)."""
    rows, header = [], None
    for line in FIXTURE.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if header is None:
            header = fields
            continue
        rows.append(dict(zip(header, fields)))
    return [{"drug_class": r.get("drug_class", "").strip(),
             "verdict": r.get("called_verdict", "").strip(),
             "phenotype": r.get("susceptibility", "").strip()} for r in rows]


def test_fixture_is_the_full_98_isolate_echinocandin_set():
    records = _records()
    assert len(records) == 98
    assert {r["drug_class"] for r in records} == {"echinocandin"}


def test_v3_diagnostic_accuracy_headline_still_holds():
    summary = ac.evaluate(_records())
    assert summary["drug_classes"] == ["echinocandin"]
    echo = summary["by_drug_class"]["echinocandin"]

    verdict, reasons = ac.apply_bar(echo, _bar())
    assert verdict == "PASS", f"diagnostic accuracy regressed: {verdict}\n" + "\n".join(reasons)

    # Scored 88 (42 R / 46 S), 10 abstained, nothing unresolved/unphenotyped — all excluded.
    assert echo["n_scored"] == 88
    assert echo["n_abstained"] == 10
    assert echo["n_unresolved"] == 0
    assert echo["n_no_phenotype"] == 0

    # Confusion from the write-up: TP=42 FN=0 FP=1 TN=45.
    c = echo["confusion"]
    assert (c["tp"], c["fn"], c["fp"], c["tn"]) == (42, 0, 1, 45)

    # The published CLSI headline: VME 0/42 = 0.0%, ME 1/46 = 2.2%.
    vme, me = echo["very_major_error"], echo["major_error"]
    assert (vme["k"], vme["n"]) == (0, 42) and vme["point"] == 0.0
    assert (me["k"], me["n"]) == (1, 46)
    assert round(me["point"], 3) == 0.022

    # Abstention 10/98 (characterised isolates with a categorical verdict).
    abst = echo["abstention"]
    assert (abst["k"], abst["n"]) == (10, 98)
