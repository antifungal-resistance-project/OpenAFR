"""Derive the #137 v3 accuracy fixture from the frozen v2 fixture (no reads, no GCP).

The v2 re-measure (work/RESULTS_diagnostic_accuracy.md) fixed under-detection (VME 0%) but
FAILED on over-calling (major-error 13.5%) because two known FKS1 positions -- D642Y and M690I
-- are LOW-PPV: their carriers split R/S rather than tracking resistance. v3 re-tiers the panel
(openafr/fks1_caller.panel_tier) so only high-PPV CORE markers detect resistance and a
low-PPV-only call abstains (UNCHARACTERIZED_VARIANT) instead of over-calling. See
work/PREREGISTRATION_diagnostic_accuracy_v3.md.

This script does NOT re-read any SRA data. It re-derives each isolate's VERDICT from the caller
output ALREADY FROZEN in the v2 fixture -- the per-isolate panel hits recorded verbatim in the
`note` column (panel_hits=...) -- by running the PRODUCTION tier-aware verdict path
(openafr.verdict.echinocandin_verdict). It is therefore a deterministic, offline transform of
frozen data: the v3 called_verdict is orchestration-derived, never hand-authored, exactly as the
prereg requires. A panel-negative isolate's verdict cannot change under a panel re-tiering, so it
is carried verbatim from v2 (and the invariant is asserted: only low-PPV-only calls flip).

Usage:
  python scripts/harvest_diagnostic_accuracy_v3.py           # writes the v3 fixture + prints a diff
  python scripts/harvest_diagnostic_accuracy_v3.py --check   # re-derive and verify, write nothing
"""
import argparse
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from openafr import fks1_caller as F        # noqa: E402
from openafr import verdict as vd           # noqa: E402

V2 = ROOT / "data" / "earlywarning" / "recaller_sanity" / "fks1_accuracy_98_v2.tsv"
V3 = ROOT / "data" / "earlywarning" / "recaller_sanity" / "fks1_accuracy_98_v3.tsv"

_PANEL_HITS_RE = re.compile(r"panel_hits=([^;]*)")

HEADER = ["isolate", "species", "run_acc", "gene", "drug_class",
          "called_verdict", "susceptibility", "note", "citation"]

V3_COMMENT = [
    "# #137 v3 verdict-accuracy fixture -- DERIVED (not re-read) from fks1_accuracy_98_v2.tsv by",
    "#   scripts/harvest_diagnostic_accuracy_v3.py. called_verdict is re-emitted by the PRODUCTION",
    "#   tier-aware echinocandin_verdict() over the SAME frozen caller output as v2 (the per-isolate",
    "#   panel hits in the v2 note column); susceptibility is carried verbatim. Only low-PPV-only",
    "#   calls (D642Y/M690I) change disposition -- RESISTANCE_MARKER_DETECTED -> UNCHARACTERIZED_VARIANT",
    "#   -- per work/PREREGISTRATION_diagnostic_accuracy_v3.md. No SRA reads, no GCP; a deterministic",
    "#   offline transform of frozen data.",
]


def _panel_hits_of(note):
    m = _PANEL_HITS_RE.search(note or "")
    if not m:
        return []
    raw = m.group(1).strip()
    if raw in ("", "-"):
        return []
    return [t for t in raw.split(",") if t]


def _read_v2():
    rows, header = [], None
    for line in V2.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if header is None:
            header = fields
            continue
        rows.append(dict(zip(header, fields)))
    return rows


def _v3_verdict(row):
    """Return (v3_verdict, tier_label) for one v2 row, re-derived via the production path."""
    hits = _panel_hits_of(row.get("note", ""))
    if not hits:
        # No panel hit -> a panel re-tiering cannot change this verdict; carry v2 verbatim.
        return row["called_verdict"].strip(), "-"
    core = [h for h in hits if F.is_core_panel_token(h)]
    # Faithful reconstruction for a panel-positive isolate: every hit is also a token. Any extra
    # non-panel tokens it might carry are irrelevant -- a core hit already wins, and a low-PPV-only
    # isolate correctly yields UNCHARACTERIZED_VARIANT from its low-PPV token in `tokens`.
    call = {"tokens": hits, "panel_hits": hits, "core_panel_hits": core, "uncalled_panel": []}
    verdict = vd.echinocandin_verdict(call=call)["verdict"]
    return verdict, ("core" if core else "low_ppv_abstain")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="re-derive and report the diff without writing the fixture")
    args = ap.parse_args()

    rows = _read_v2()
    out_lines = list(V3_COMMENT)
    out_lines.append("\t".join(HEADER))
    flips = []
    for r in rows:
        v3v, tier = _v3_verdict(r)
        v2v = r["called_verdict"].strip()
        if v3v != v2v:
            flips.append((r["isolate"], _panel_hits_of(r.get("note", "")), v2v, v3v))
        note = r.get("note", "")
        note = f"{note};v3_tier={tier}" if note else f"v3_tier={tier}"
        out_lines.append("\t".join([
            r["isolate"], r["species"], r["run_acc"], r["gene"], r["drug_class"],
            v3v, r["susceptibility"], note, r["citation"],
        ]))

    print(f"rows: {len(rows)}   verdicts changed v2 -> v3: {len(flips)}", file=sys.stderr)
    for iso, hits, v2v, v3v in flips:
        print(f"  {iso}  {','.join(hits) or '-'}: {v2v} -> {v3v}", file=sys.stderr)
    # Invariant: only low-PPV-only calls may flip, and only from detected -> abstain.
    for iso, hits, v2v, v3v in flips:
        assert v2v == vd.RESISTANCE_MARKER_DETECTED and v3v == vd.UNCHARACTERIZED_VARIANT, \
            f"unexpected flip for {iso}: {v2v} -> {v3v}"
        assert hits and all(not F.is_core_panel_token(h) for h in hits), \
            f"{iso} flipped but carries a core marker {hits}"

    if args.check:
        print("check-only: no file written", file=sys.stderr)
        return 0
    V3.write_text("\n".join(out_lines) + "\n")
    print(f"wrote {V3.relative_to(ROOT)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
