"""CI regression guard for the certified FKS1 concordance headline (issue #164, C3).

The FKS1 re-caller was certified against the 98-isolate PMC12323592 benchmark at
**100% sensitivity / 100% specificity / 100% exact-token identity / 98-of-98 resolved ->
pre-registered verdict PASS** (work/RESULTS_fks1_concordance.md, PR #151). This test re-scores
that result on every PR so a refactor of the scoring code, a slip in the pre-committed pass bar,
or an edit to a frozen artifact cannot silently invalidate the published claim.

SCOPE -- what this guards and what it does NOT:
  * GUARDS the concordance SCORING path (openafr.concordance.evaluate + apply_bar), the frozen
    pass bar (FKS1_CONCORDANCE_BAR), and the frozen per-isolate calls the headline is built from.
  * Does NOT re-run the aligner/caller. The real concordance run needs SRA downloads + alignment
    tools and cannot run in a stock GH Action; it lives in the cloud grader
    (scripts/validate_fks1_concordance.py) + work/RUNBOOK_fks1_run.md. The caller's
    genotype->panel LOGIC is covered offline by tests/test_fks1_caller.py.

Source of the per-isolate calls: the certified run's own runlog
(data/earlywarning/runlog/concordance-fks1.jsonl) -- the full 98-isolate PASS record, whose
`items` carry the exact HS1-panel `called`/`expected`/`resolved`/`resistant` the caller emitted.
We re-score those through the production scorer and bar. As an anti-tamper tie, each item's
`expected`/`resistant` is cross-checked against the frozen benchmark truth set. No network, no reads.
"""
import json
import pathlib

from openafr import concordance as cc

ROOT = pathlib.Path(__file__).resolve().parent.parent
SANITY = ROOT / "data" / "earlywarning" / "recaller_sanity"
BENCHMARK = SANITY / "fks1_benchmark_100.tsv"
RUNLOG = ROOT / "data" / "earlywarning" / "runlog" / "concordance-fks1.jsonl"

N_ISOLATES = 98


def _read_tsv(path):
    rows, header = [], None
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if header is None:
            header = fields
            continue
        rows.append(dict(zip(header, fields)))
    return rows


def _certified_run():
    """The full 98-isolate benchmark record from the concordance runlog (not a smoke run)."""
    records = [json.loads(line) for line in RUNLOG.read_text().splitlines() if line.strip()]
    full = [r for r in records
            if r.get("fixture") == BENCHMARK.name and r.get("n_isolates") == N_ISOLATES]
    assert len(full) == 1, (
        f"expected exactly one full {N_ISOLATES}-isolate concordance run in the runlog, "
        f"found {len(full)}")
    return full[0]


def test_runlog_matches_the_frozen_benchmark_truth_set():
    # Anti-tamper: the certified run's per-isolate expected/phenotype must match the committed
    # benchmark, so the replay cannot be gamed by editing the runlog in isolation.
    run = _certified_run()
    bench = {r["run_acc"]: r for r in _read_tsv(BENCHMARK)}
    assert len(bench) == N_ISOLATES
    items = {it["run_acc"]: it for it in run["items"]}
    assert items.keys() == bench.keys(), "runlog isolates differ from the benchmark truth set"
    for acc, b in bench.items():
        it = items[acc]
        assert cc._norm(it["expected"]) == cc._norm(b["expected_panel"]), acc
        assert bool(it["resistant"]) == (b["susceptibility"].strip().upper() == "R"), acc


def test_certified_concordance_headline_still_holds():
    run = _certified_run()
    summary = cc.evaluate(run["items"])
    verdict, reasons = cc.apply_bar(summary)

    assert verdict == "PASS", f"concordance regressed: {verdict}\n" + "\n".join(reasons)

    # Resolution: 0 unresolved, all 98 scored (work/RESULTS_fks1_concordance.md).
    assert summary["n_total"] == N_ISOLATES
    assert summary["n_resolved"] == N_ISOLATES
    assert summary["n_unresolved"] == 0

    # Published headline -- 100% on every gated metric.
    assert summary["sensitivity"]["point"] == 1.0
    assert summary["specificity"]["point"] == 1.0
    assert summary["identity"]["point"] == 1.0

    # Confusion counts from the write-up (sens 37/37, spec 61/61, zero false calls).
    c = summary["confusion"]
    assert (c["tp"], c["tn"], c["fp"], c["fn"]) == (37, 61, 0, 0)
