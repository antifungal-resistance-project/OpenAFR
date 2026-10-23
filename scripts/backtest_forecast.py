"""Grade the FKS1 prevalence forecast against its pre-registration, and write RESULTS.

Pre-registration: work/PREREGISTRATION_fks1_forecast_v2.md (the precision-weighted model;
v1 was work/PREREGISTRATION_fks1_forecast.md), sha frozen in
work/PREREG_fks1_forecast_v2.sha256. This grader re-checks that hash and REFUSES to certify if
the prereg moved after it was frozen (the same integrity guard validate_fks1_concordance.py
uses) -- so the bar can never be edited to fit the result.

Fully local: reads the committed release-710 call table offline, builds the panel-gated
monthly prevalence series (openafr.forecast), runs the walk-forward skill test, and writes
work/RESULTS_fks1_forecast.md. No NCBI pull, no GCP, no network.

Usage
-----
  python scripts/backtest_forecast.py                 # default table + RESULTS path
  python scripts/backtest_forecast.py --table <tsv> --out <md>
  python scripts/backtest_forecast.py --check         # integrity gate only, no write
"""
import argparse
import datetime
import hashlib
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from openafr import earlywarning as ew   # noqa: E402
from openafr import forecast             # noqa: E402

PREREG = ROOT / "work" / "PREREGISTRATION_fks1_forecast_v2.md"
PINS = ROOT / "work" / "PREREG_fks1_forecast_v2.sha256"
DEFAULT_TABLE = ROOT / "data" / "earlywarning" / "fks1_calls" / "PDG000000067.710.tsv"
DEFAULT_OUT = ROOT / "work" / "RESULTS_fks1_forecast.md"


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _check_integrity():
    """Re-check the frozen prereg hash. Returns (ok, messages). An unpinned prereg fails
    closed with the exact freeze command, so a first run cannot certify against nothing."""
    import os
    pinned = {}
    if PINS.exists():
        for line in PINS.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            digest, _, name = line.partition("  ")
            pinned[name.strip()] = digest.strip()
    rel = os.path.relpath(PREREG, ROOT)
    want = pinned.get(rel)
    got = _sha256(PREREG) if PREREG.exists() else None
    if got is None:
        return False, [f"MISSING: {rel} does not exist"]
    if want is None:
        return False, [f"UNPINNED: {rel} has no frozen hash in {PINS.name} "
                       f"(freeze it: `shasum -a 256 {rel} > work/{PINS.name}`)"]
    if want != got:
        return False, [f"HASH MOVED: {rel}\n    frozen {want}\n    now    {got}"]
    return True, [f"ok: {rel} matches frozen hash"]


def _pct(x):
    return "—" if x is None else f"{100 * x:.2f}%"


def _pp(x):
    return "—" if x is None else f"{100 * x:.2f} pp"


def _render_results(table, table_sha, series, bt, fc, skill):
    now = datetime.date.today().isoformat()
    series_rows = "\n".join(
        f"| {e['month']} | {e['n_resistant']} / {e['n_resolved']} | {_pct(e['prevalence'])} |"
        for e in series)
    fold_rows = "\n".join(
        f"| {f['train_end']} → {f['target_month']} | {_pct(f['actual'])} | "
        f"{_pct(f['model'])} | {_pct(f['naive_last'])} | {_pct(f['naive_pooled'])} | "
        f"{'✓' if f['in_interval'] else '✗'} |"
        for f in bt["folds"])
    verdict = "SKILL" if skill else "NO SKILL — trajectory only"
    if skill:
        fc_line = (f"One-step projection ({fc['next_month']}): **{_pct(fc['point'])}** "
                   f"(95% PI {_pct(fc['lo'])}–{_pct(fc['hi'])}).")
    else:
        pooled = fc["naive_pooled"] if fc else None
        fc_line = ("No forecast skill over a flat baseline at this series depth; the page shows "
                   f"the monthly trajectory only. Pooled prevalence ≈ {_pct(pooled)}.")
    return f"""# Results — FKS1 (echinocandin) prevalence forecast

**Run {now}.** Pre-registered in [PREREGISTRATION_fks1_forecast_v2.md](PREREGISTRATION_fks1_forecast_v2.md)
(hash re-checked by the grader before this file was written). Fully local — committed call
table read offline, no GCP.

- **Source table:** `{table}`
  sha256 `{table_sha}`
- **Series depth:** {len(series)} deposit-month bins ({series[0]['month']} → {series[-1]['month']}).

## P1 — monthly panel-gated prevalence (deposit date)

| Month | Panel-positive / resolved | Prevalence |
|---|---|---|
{series_rows}

Panel-gated (known-resistance tokens only) over the resolved denominator. Deposit-month axis,
**not** epidemiological incidence (see prereg limitation 1).

## P2 — walk-forward skill test

Pre-committed bar: model MAE must be strictly below **both** naive baselines.

| Metric | Value |
|---|---|
| Held-out folds | {bt['n_folds']} (min train {bt['min_train']}, horizon {bt['horizon']}) |
| Model MAE | {_pp(bt['model_mae'])} |
| Naive-last MAE | {_pp(bt['naive_last_mae'])} |
| Naive-pooled MAE | {_pp(bt['naive_pooled_mae'])} |

| Fold (train → target) | Actual | Model | Naive-last | Naive-pooled | In 95% PI |
|---|---|---|---|---|---|
{fold_rows}

## P3 — interval coverage

Actual-in-95%-PI over held-out months: **{_pct(bt['interval_coverage'])}** (descriptive; 7 bins
underpower a coverage *claim*).

## Verdict: {verdict}

{fc_line}

A deeper multi-release series (which would give the trend model something to fit) remains
blocked on repeated billable fills — the single named remaining blocker.
"""


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", default=str(DEFAULT_TABLE), help="committed FKS1 call table")
    ap.add_argument("--out", default=str(DEFAULT_OUT), help="RESULTS markdown path")
    ap.add_argument("--check", action="store_true",
                    help="run the prereg integrity gate only; do not grade or write")
    args = ap.parse_args()

    ok, msgs = _check_integrity()
    for m in msgs:
        print(m, file=sys.stderr)
    if not ok:
        print("REFUSING to certify: the pre-registration is missing/unpinned/moved.",
              file=sys.stderr)
        return 2
    if args.check:
        print("integrity OK.", file=sys.stderr)
        return 0

    table = pathlib.Path(args.table)
    records = list(ew.read_fks1_call_table(str(table)).values())
    series = forecast.monthly_prevalence_series(records)
    if len(series) < 2:
        print(f"ERROR: only {len(series)} month bin(s) — too short to forecast.",
              file=sys.stderr)
        return 1
    bt = forecast.backtest(series, min_train=3, horizon=1)
    fc = forecast.forecast_next(series, horizon=1)
    skill = forecast.has_skill(bt)

    out = pathlib.Path(args.out)
    # Render the table path repo-relative when it lives under ROOT, so the committed
    # artifact is identical no matter where the grader was invoked from.
    try:
        table_label = table.resolve().relative_to(ROOT)
    except ValueError:
        table_label = table
    out.write_text(_render_results(table_label, _sha256(table), series, bt, fc, skill))
    print(f"wrote {out} — verdict: {'SKILL' if skill else 'NO SKILL (trajectory only)'} "
          f"(model MAE {_pp(bt['model_mae'])} vs last {_pp(bt['naive_last_mae'])} / "
          f"pooled {_pp(bt['naive_pooled_mae'])}).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
