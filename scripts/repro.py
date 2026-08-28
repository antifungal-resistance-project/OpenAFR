"""One-command reproducibility + drift harness (issue #89).

The project's edge is provable discipline: a hash-frozen protocol, hash-frozen
pre-registrations, and a confidence stack that derives every published number from
committed inputs. As the candidate-confidence batch (#65-#88) piled axis on axis, a
single silent edit — a tuned pass bar in protocol.yaml, a pre-registration quietly
reworded after the fact, a per-axis JSON regenerated under changed logic — could break
a claim without anyone noticing. This harness re-checks the whole frozen surface in one
command and FAILS LOUDLY (exit 1) on any drift, so CI can gate every PR on it.

Three checks, each independent:

  1. PROTOCOL FREEZE. protocol.yaml (the docking box + search effort + the gate's pass
     bar — the numbers easiest to tune into a pass) matches its frozen SHA-256.

  2. PRE-REGISTRATION FREEZE. Every work/PREREG_*.sha256 line verifies — the
     pre-registration markdown AND the frozen input datasets each one pins (held-out
     SMILES, decoy sets). A mismatch means the rules or the data changed after freezing,
     so any "pass" graded against them is no longer credible.

  3. DOCKING-FREE DERIVATION. The consolidated per-candidate scorecard
     (work/candidates/shortlist_scorecard.json — the capstone artifact the whole stack
     feeds) is re-derived from its committed per-axis inputs and byte-compared. If any
     committed per-axis JSON drifted, or the merge logic changed, the bytes differ.

DELIBERATELY OUT OF SCOPE: re-docking. The ~43 MB of Vina pose trees are gitignored and
regenerate from scripts/ (see each RESULTS_*.md "Reproduce" block); re-running them needs
the docking toolkit and hours of compute, not a per-PR CI job. The pose-DEPENDENT per-axis
JSONs are treated here as frozen committed inputs — check 3 proves the scorecard still
derives from them unchanged, check 2 proves the rules that graded them were not moved.
This harness therefore guarantees no *silent* drift in anything that does not require
docking; a genuine re-dock is a separate, documented, manual step.

Usage:
    conda run -n openafr python scripts/repro.py          # exits 1 on any drift
    conda run -n openafr python scripts/repro.py --quiet   # only the final verdict + failures
"""
import glob
import hashlib
import json
import os
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from openafr import protocol

PASS = "PASS"
FAIL = "FAIL"


def _sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def check_protocol():
    """protocol.yaml matches its frozen SHA-256 (openafr.protocol's own guarantee)."""
    lines = []
    ok = _sha256(protocol.PROTOCOL) == protocol.PROTOCOL_SHA
    if ok:
        lines.append(f"  protocol.yaml unmodified (sha256 {protocol.PROTOCOL_SHA[:16]}...)")
    else:
        lines.append("  !! protocol.yaml HAS DRIFTED from its frozen SHA-256 !!")
        lines.append(f"     expected {protocol.PROTOCOL_SHA[:16]}...  got {_sha256(protocol.PROTOCOL)[:16]}...")
    return ok, lines


def _parse_sha256_file(path):
    """Yield (expected_hash, target_path) from a `sha256sum`-format file."""
    for raw in pathlib.Path(path).read_text().splitlines():
        raw = raw.strip()
        if not raw:
            continue
        expected, _, target = raw.partition("  ")   # sha256sum uses two spaces
        yield expected.strip(), target.strip()


def check_preregistrations():
    """Every work/PREREG_*.sha256 verifies: the pre-registration doc + its frozen inputs."""
    lines = []
    sha_files = sorted(glob.glob(str(ROOT / "work" / "PREREG_*.sha256")))
    if not sha_files:
        return False, ["  no work/PREREG_*.sha256 files found — expected 18+"]
    n_entries = 0
    drifted = []
    for sf in sha_files:
        for expected, target in _parse_sha256_file(sf):
            n_entries += 1
            tp = ROOT / target
            if not tp.exists():
                drifted.append(f"{target} (MISSING; pinned by {os.path.basename(sf)})")
            elif _sha256(tp) != expected:
                drifted.append(f"{target} (edited since freeze; pinned by {os.path.basename(sf)})")
    ok = not drifted
    if ok:
        lines.append(f"  {len(sha_files)} pre-registrations / {n_entries} frozen files all unmodified")
    else:
        lines.append(f"  !! {len(drifted)} of {n_entries} frozen files DRIFTED !!")
        lines += [f"     - {d}" for d in drifted]
    return ok, lines


def check_scorecard_reproduces():
    """The consolidated scorecard re-derives byte-for-byte from committed per-axis inputs.

    Runs the builder as a subprocess to a scratch path (so a broken import surfaces the
    same way it would for a user) and byte-compares against the committed card.
    """
    lines = []
    committed = ROOT / "work" / "candidates" / "shortlist_scorecard.json"
    if not committed.exists():
        return False, ["  work/candidates/shortlist_scorecard.json missing"]
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tmp:
        out = tmp.name
    try:
        proc = subprocess.run(
            [sys.executable, str(ROOT / "work" / "candidates" / "build_scorecard.py"), "--out", out],
            capture_output=True, text=True)
        if proc.returncode != 0:
            return False, ["  build_scorecard.py failed to run:", *("     " + l for l in proc.stderr.strip().splitlines()[-5:])]
        regenerated = pathlib.Path(out).read_bytes()
    finally:
        os.unlink(out)
    if regenerated == committed.read_bytes():
        lines.append("  scorecard re-derives byte-for-byte from committed per-axis inputs")
        return True, lines
    lines.append("  !! shortlist_scorecard.json does NOT re-derive from its committed inputs !!")
    lines.append("     a per-axis JSON or the merge logic changed without rebuilding the card")
    return False, lines


CHECKS = [
    ("1. Frozen protocol integrity", check_protocol),
    ("2. Frozen pre-registration integrity", check_preregistrations),
    ("3. Docking-free derivation reproduces", check_scorecard_reproduces),
]


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    quiet = "--quiet" in argv

    results = []
    if not quiet:
        print("=" * 78)
        print("OPENAFR REPRODUCIBILITY / DRIFT HARNESS (#89)")
        print("=" * 78)
    for title, fn in CHECKS:
        ok, lines = fn()
        results.append((title, ok))
        if not quiet or not ok:
            print(f"[{PASS if ok else FAIL}] {title}")
            for l in lines:
                print(l)
    all_ok = all(ok for _, ok in results)
    print("-" * 78)
    if all_ok:
        print(f"{PASS}: all {len(results)} checks green — no silent drift in the frozen surface.")
    else:
        failed = [t for t, ok in results if not ok]
        print(f"{FAIL}: {len(failed)} check(s) drifted: {'; '.join(failed)}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
