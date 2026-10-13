"""Pre-stage the FKS1 fill -> weather-page render path at full scale (issue #171).

Why this exists
---------------
The public Resistance-Weather page shows a real echinocandin picture by overlaying a
*committed* FKS1 call table (`data/earlywarning/fks1_calls/PDG...<N>.tsv`) onto a live NCBI
pull. That overlay path has only ever been exercised against the first fill (release 709,
**300 isolates**). The next fill is the full trailing-180-day window -- release 710,
**~3375 isolates**. That run is **billable GCP compute**: if `render_weather.py` or
`earlywarning.overlay_fks1_calls` chokes on the real table's *scale* or its isolate-key
*match rate*, we want to know now, against a synthetic table, not after paying for the fill.

This harness manufactures a release-710-shaped synthetic world -- a raw-NCBI `--tsv` of N
pending isolates plus a matching committed call table -- pushes it through the **real**
`scripts/render_weather.py` CLI exactly as the GitHub Action does (as-of = today-180, window
180), and asserts the page leaves the day-0 "watching" state with an honest coverage clause.
No cloud, no billable compute, nothing committed: everything lives in a temp dir.

It reuses the producer's own serializers (`ew.write_fks1_call_table`,
`ew.update_fks1_calls_index`) so the synthetic table can never drift from the schema a real
`recall_fks1.py fill` emits -- this doubles as the "confirm the schema consumers expect"
step (#171 step 1).

Usage
-----
  python scripts/prestage_weather_dryrun.py                 # 3375 isolates, ~6% resistant
  python scripts/prestage_weather_dryrun.py --n 300         # mimic the 709 fill scale
  python scripts/prestage_weather_dryrun.py --resistant-frac 0   # all-wild-type fill
  python scripts/prestage_weather_dryrun.py --keep          # leave the temp dir for eyeballing

Exit status is 0 on PASS, 1 on FAIL -- so it can gate a release.

Stdlib + openafr.earlywarning only, matching the track's house style.
"""
import argparse
import csv
import datetime
import pathlib
import re
import subprocess
import sys
import tempfile
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from openafr import earlywarning as ew  # noqa: E402

# Raw-NCBI metadata columns render_weather's --tsv path feeds to ew.normalize(). A live pull
# carries NO call -- normalize seeds every row pending -- so these are metadata-only.
_RAW_COLUMNS = ("target_acc", "Run", "biosample_acc", "asm_acc",
                "collection_date", "target_creation_date", "geo_loc_name", "lat_lon")

# Hotspot calls + their re-caller source strings, mirroring the committed 709 table's shape
# (sra-fks1-recaller:HS1=...,HS2=...,HS3=...). Which hotspot a call sits in is cosmetic here;
# the render path only cares that the source is non-pending and the call is non-empty.
_HOTSPOTS = [
    ("D642Y", "HS1"),
    ("S639F", "HS1"),
    ("M690I", "HS3"),
    ("R1354S", "HS3"),
]
_WILD_SOURCE = "sra-fks1-recaller:HS1=wild-type,HS2=wild-type,HS3=wild-type"

_GEOS = ["USA:California", "USA:Texas", "India:Delhi", "South Africa:Gauteng",
         "United Kingdom:England", "Spain:Madrid", ""]


def _called_source(hs):
    """The re-caller source string for a resistant call in hotspot `hs`."""
    return ",".join(f"{h}={'called' if h == hs else 'wild-type'}"
                    for h in ("HS1", "HS2", "HS3"))


def _spread_dates(n, as_of, window_days):
    """n target_creation_dates spread across the recent window (as_of, as_of+window_days].

    Mirrors the detector's window semantics (openafr/emergence.split_by_creation): a record is
    'recent' iff created strictly after the baseline cut and within window_days of it. The
    workflow sets as_of = today-180, so this window is the trailing 180 days up to today."""
    base = datetime.date.fromisoformat(as_of)
    # Offsets 1..window_days inclusive, cycled, so every row lands inside the recent window.
    return [(base + datetime.timedelta(days=1 + (i % window_days))).isoformat()
            for i in range(n)]


def synthesize(n, resistant_frac, as_of, window_days):
    """Build a synthetic population: a list of raw-NCBI rows and the matching resolved call
    records. Resistant isolates get a hotspot call; the rest an all-wild-type (but resolved)
    source so the whole window reads as re-called. isolate_key is the version-stripped
    target_acc, so the call table overlays cleanly onto the normalized pull."""
    dates = _spread_dates(n, as_of, window_days)
    n_resistant = round(n * resistant_frac)
    raw_rows, call_records = [], []
    for i in range(n):
        acc = f"PDT{i:08d}"
        target_acc = f"{acc}.1"
        run_acc = f"SRR{40000000 + i}"
        created = dates[i]
        raw_rows.append({
            "target_acc": target_acc,
            "Run": run_acc,
            "biosample_acc": f"SAMN{30000000 + i}",
            "asm_acc": "",
            "collection_date": created,
            "target_creation_date": created,
            "geo_loc_name": _GEOS[i % len(_GEOS)],
            "lat_lon": "",
        })
        if i < n_resistant:
            call, hs = _HOTSPOTS[i % len(_HOTSPOTS)]
            source = f"sra-fks1-recaller:{_called_source(hs)}"
        else:
            call, source = "", f"sra-fks1-recaller:{_WILD_SOURCE.split(':', 1)[1]}"
        call_records.append({
            "isolate_key": acc,
            "fks1_call": call,
            "fks1_resistance_source": source,
            "run_acc": run_acc,
            "target_creation_date": created,
        })
    return raw_rows, call_records


def _write_raw_tsv(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(_RAW_COLUMNS), delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _run_render(tsv, calls_table, as_of, window_days, out):
    """Invoke the real render CLI exactly as the workflow does; return (proc, elapsed_s)."""
    cmd = [sys.executable, str(ROOT / "scripts" / "render_weather.py"),
           "--gene", "fks1", "--tsv", str(tsv), "--fks1-calls", str(calls_table),
           "--as-of", as_of, "--window", str(window_days), "--out", str(out)]
    t0 = time.monotonic()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc, time.monotonic() - t0


def _check(cond, label, failures):
    mark = "PASS" if cond else "FAIL"
    print(f"  [{mark}] {label}")
    if not cond:
        failures.append(label)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=3375,
                    help="synthetic isolate count (default: 3375, the release-710 full window)")
    ap.add_argument("--resistant-frac", type=float, default=0.06,
                    help="fraction carrying a hotspot call (default: 0.06, ~19/300 like 709)")
    default_as_of = (datetime.date.today() - datetime.timedelta(days=180)).isoformat()
    ap.add_argument("--as-of", default=default_as_of,
                    help="baseline cut date (default: today-180, matching the workflow)")
    ap.add_argument("--window", type=int, default=180, help="recent-window days (default: 180)")
    ap.add_argument("--release", type=int, default=710,
                    help="PDG release number for the synthetic fill (default: 710)")
    ap.add_argument("--keep", action="store_true",
                    help="keep the temp dir (print its path) instead of cleaning up")
    args = ap.parse_args()

    n_resistant = round(args.n * args.resistant_frac)
    print(f"Pre-staging render: n={args.n} isolates, {n_resistant} resistant, "
          f"release {args.release}, as-of {args.as_of} (+{args.window}d window).")

    tmp = pathlib.Path(tempfile.mkdtemp(prefix="prestage-weather-"))
    calls_dir = tmp / "fks1_calls"
    calls_dir.mkdir()
    tsv = tmp / "synthetic_pull.tsv"
    calls_table = calls_dir / f"{ew.ORG_PREFIX}.{args.release}.tsv"
    out = tmp / "index.html"

    raw_rows, call_records = synthesize(args.n, args.resistant_frac, args.as_of, args.window)
    _write_raw_tsv(tsv, raw_rows)
    # Build the call table + INDEX via the SAME serializers recall_fks1.py fill uses, so the
    # synthetic table cannot drift from the schema the consumers expect (#171 step 1).
    sha256, n_called = ew.write_fks1_call_table(str(calls_table), call_records)
    ew.update_fks1_calls_index(str(calls_dir / "INDEX.tsv"), {
        "release_tag": f"{ew.ORG_PREFIX}.{args.release}",
        "as_of": args.as_of, "window_days": args.window, "n_called": n_called,
        "sha256": sha256,
        "written_at_utc": datetime.datetime.now(datetime.timezone.utc)
                                  .strftime("%Y-%m-%dT%H:%M:%SZ"),
    })

    proc, elapsed = _run_render(tsv, calls_table, args.as_of, args.window, out)

    print(f"\nrender_weather exited {proc.returncode} in {elapsed:.1f}s.")
    if proc.stderr.strip():
        print("  stderr:", proc.stderr.strip().replace("\n", "\n          "))

    failures = []
    _check(proc.returncode == 0, "render exits cleanly (never 500s)", failures)

    overlaid = None
    m = re.search(r"overlaid (\d+) FKS1 call", proc.stderr)
    if m:
        overlaid = int(m.group(1))
    _check(overlaid == args.n,
           f"overlay matched every isolate ({overlaid} of {args.n} by isolate_key)", failures)

    page = out.read_text() if out.exists() else ""
    _check(bool(page), "page was written", failures)

    # "Leaves day-0" = an honest coverage clause naming the re-called window, never the calm
    # "watching since" day-0 copy. render_weather renders "N of M window isolates re-called".
    cov_clause = re.search(r"(\d[\d,]*) of (\d[\d,]*) window isolates re-called", page)
    _check(bool(cov_clause), "page shows a coverage clause (not day-0 watching)", failures)
    if cov_clause:
        called_shown = int(cov_clause.group(1).replace(",", ""))
        total_shown = int(cov_clause.group(2).replace(",", ""))
        _check(total_shown == args.n and called_shown == args.n,
               f"coverage reads {called_shown}/{total_shown} = full window re-called", failures)

    _check(f"calls current through" in page and args.as_of in page,
           f"freshness clause reports the synthetic fill date ({args.as_of})", failures)

    if n_resistant > 0:
        _check("resistance" in page.lower() and "card" in page.lower() or "alert" in page.lower(),
               "page surfaces the resistance signal", failures)

    # Scale sanity: a full-window render should be well under a minute; flag if pathological.
    _check(elapsed < 60, f"render wall-time is sane at {args.n}-row scale ({elapsed:.1f}s)",
           failures)

    print()
    if failures:
        print(f"FAIL — {len(failures)} check(s) failed: {'; '.join(failures)}")
        if args.keep:
            print(f"temp dir kept: {tmp}")
        return 1
    print(f"PASS — the {args.n}-isolate release-{args.release} fill renders a real "
          f"echinocandin picture; the overlay/render/freshness path is proven at full scale.")
    if args.keep:
        print(f"temp dir kept (eyeball {out}): {tmp}")
    else:
        import shutil
        shutil.rmtree(tmp, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
