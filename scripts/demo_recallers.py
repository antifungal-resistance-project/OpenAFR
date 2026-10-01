"""One-command reproducible demo of both re-callers (issue C1, #162).

The on-ramp for an external user: run each re-caller on a tiny committed sample and watch
it produce a call -- no conda env, no network, no bioinformatics tools, any platform
(incl. osx-arm64). It drives the deterministic `call` subcommand of the two CLIs
(`scripts/recall_erg11.py`, `scripts/recall_fks1.py`) exactly as a human would type it, on
the four sample consensus FASTAs in `data/earlywarning/demo/` (two per gene: a wild-type
susceptible control and a single known resistance mutation). See that directory's README
for how each sample is derived from the pinned reference CDS.

This covers ONLY the reproducible-anywhere `call` path. The reads->consensus orchestration
(`recall`/`fill`/sanity) needs the bioconda toolchain + multi-GB SRA downloads and is
documented in the RUNBOOKs, not here.

Usage:
  python scripts/demo_recallers.py           # run all four samples, print a table
  python scripts/demo_recallers.py --check    # also verify each against EXPECTED.txt (CI)

With --check the script compares the emitted call/source/panel lines against the committed
`data/earlywarning/demo/EXPECTED.txt` oracle and exits non-zero on any mismatch, so the
samples and the callers can never silently drift apart.
"""
import argparse
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
DEMO_DIR = ROOT / "data" / "earlywarning" / "demo"
EXPECTED_PATH = DEMO_DIR / "EXPECTED.txt"

# (label, caller CLI, sample FASTA) -- the four demo runs, in table order.
CASES = [
    ("erg11_wildtype", "recall_erg11.py", "erg11_wildtype.fasta"),
    ("erg11_Y132F", "recall_erg11.py", "erg11_Y132F.fasta"),
    ("fks1_wildtype", "recall_fks1.py", "fks1_wildtype.fasta"),
    ("fks1_S639F", "recall_fks1.py", "fks1_S639F.fasta"),
]

# The three stdout lines of the `call` subcommand we treat as the call's identity. The FKS1
# per-window detail lines are informational and not part of the checked oracle.
KEY_PREFIXES = ("call:", "source:", "panel:")


def run_case(script, fasta):
    """Invoke a re-caller's deterministic `call` on one sample; return its full stdout."""
    cmd = [
        sys.executable, str(ROOT / "scripts" / script),
        "call", "--consensus-fasta", str(DEMO_DIR / fasta),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise SystemExit(f"re-caller failed on {fasta} (exit {proc.returncode})")
    return proc.stdout


def key_lines(stdout):
    """The call/source/panel identity lines of a `call` run, normalised to 'prefix value'."""
    out = {}
    for line in stdout.splitlines():
        s = line.strip()
        for p in KEY_PREFIXES:
            if s.startswith(p):
                out[p] = " ".join(s.split())
    return out


def results():
    """Run every case; return {label: {'call:': ..., 'source:': ..., 'panel:': ...}}."""
    return {label: key_lines(run_case(script, fasta)) for label, script, fasta in CASES}


def format_expected(res):
    """Serialise results to the EXPECTED.txt format: a `## label` block per case."""
    blocks = []
    for label, _script, _fasta in CASES:
        lines = [f"## {label}"]
        for p in KEY_PREFIXES:
            lines.append(res[label][p])
        blocks.append("\n".join(lines))
    return "\n\n".join(blocks) + "\n"


def parse_expected(text):
    """Inverse of format_expected: EXPECTED.txt text -> {label: {prefix: 'prefix value'}}."""
    out, label = {}, None
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("## "):
            label = s[3:].strip()
            out[label] = {}
        elif label and any(s.startswith(p) for p in KEY_PREFIXES):
            prefix = next(p for p in KEY_PREFIXES if s.startswith(p))
            out[label][prefix] = " ".join(s.split())
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="verify each run against data/earlywarning/demo/EXPECTED.txt; "
                         "exit non-zero on any mismatch")
    args = ap.parse_args()

    res = results()

    print(f"{'sample':<16} {'call':<10} {'source'}")
    print("-" * 72)
    for label, _script, _fasta in CASES:
        r = res[label]
        call = r["call:"].split(":", 1)[1].strip()
        if call.startswith("(none"):
            call = "-"
        source = r["source:"].split(":", 1)[1].strip()
        print(f"{label:<16} {call:<10} {source}")

    if not args.check:
        print(f"\nSamples: {DEMO_DIR.relative_to(ROOT)}/  "
              f"(--check verifies these against EXPECTED.txt)")
        return 0

    expected = parse_expected(EXPECTED_PATH.read_text())
    failures = []
    print()
    for label, _script, _fasta in CASES:
        got, want = res[label], expected.get(label, {})
        if got == want:
            print(f"PASS  {label}")
        else:
            print(f"FAIL  {label}")
            for p in KEY_PREFIXES:
                if got.get(p) != want.get(p):
                    print(f"        expected {want.get(p)!r}")
                    print(f"        got      {got.get(p)!r}")
            failures.append(label)

    if failures:
        print(f"\n{len(failures)} sample(s) drifted from EXPECTED.txt: "
              f"{', '.join(failures)}")
        return 1
    print(f"\nAll {len(CASES)} samples match EXPECTED.txt.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
