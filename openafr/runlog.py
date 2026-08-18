"""Append-only provenance log for the ERG11 re-caller's un-reproducible runs.

The two tool+network paths -- the sanity harness (`scripts/recaller_sanity.py`) and
`scripts/recall_erg11.py {recall,fill}` -- shell out to fasterq-dump/minimap2/samtools
and either only print (sanity) or mutate a snapshot in place (fill). Without a log each
run's outcome and provenance vanish when the terminal scrolls, and a re-`fill` silently
overwrites the last one's calls. This records one JSON object per run to an append-only
`.jsonl`, capturing enough to answer months later: *what did that run do, on which inputs,
with which external tools, against which code and reference?*

Same spirit as the git-tracked `snapshots/INDEX.tsv` (openafr.earlywarning.update_index):
the big regenerable blobs stay out of git; the small durable memory of "what happened" is
committed. One `.jsonl` line per run is tiny, so the run logs are meant to be committed.

Design contract -- **logging must never fail or alter the science run**:
  * Every provenance probe (git commit, tool versions, file digest) is best-effort and
    returns None on any error rather than raising.
  * `record_run` is a context manager that writes the entry on exit *even if the run body
    raised* -- a crashed run is logged with whatever items accumulated plus the error, not
    lost. It never suppresses the exception (it re-raises).
  * A failure to append the log line itself is swallowed (reported to stderr), so a
    read-only filesystem degrades tracking without breaking the run.

Stdlib only (json + hashlib + subprocess + platform), matching the rest of openafr.
"""
import contextlib
import hashlib
import json
import os
import pathlib
import platform
import subprocess
import sys
import uuid

from openafr import earlywarning as _ew  # utcnow_iso(): the one timestamp format we use

_ROOT = pathlib.Path(__file__).resolve().parent.parent
DEFAULT_LOG_DIR = _ROOT / "data" / "earlywarning" / "runlog"

# Schema version of a run-log line. Bump when the object shape changes so a reader can
# tell old lines from new ones. Append-only files will contain a mix across time.
SCHEMA = 1


# --------------------------- best-effort provenance -------------------------

def _first_line(text):
    for line in (text or "").splitlines():
        line = line.strip()
        if line:
            return line
    return None


def git_commit(root=_ROOT):
    """(short_sha, dirty) for the working tree, or (None, None) if git is unavailable."""
    try:
        sha = subprocess.run(["git", "-C", str(root), "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=5, check=True).stdout.strip()
        porcelain = subprocess.run(["git", "-C", str(root), "status", "--porcelain"],
                                   capture_output=True, text=True, timeout=5, check=True).stdout
        return (sha or None), bool(porcelain.strip())
    except Exception:
        return None, None


def tool_version(tool):
    """First line of `<tool> --version`, or None if the tool is absent/errors."""
    try:
        out = subprocess.run([tool, "--version"], capture_output=True, text=True, timeout=5)
        # Some tools print the banner to stderr (or exit nonzero while still reporting).
        return _first_line(out.stdout) or _first_line(out.stderr)
    except Exception:
        return None


def tool_versions(tools):
    return {t: tool_version(t) for t in tools}


def file_digest(path):
    """sha256 of a file's bytes (e.g. the reference CDS), or None if unreadable."""
    try:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                h.update(chunk)
        return h.hexdigest()
    except Exception:
        return None


# ------------------------------- the record ---------------------------------

class RunRecord:
    """Accumulates one run's outcome. Populate `items` as the run proceeds and call
    `set(...)` for the run-level summary; the context manager writes it out on exit."""

    def __init__(self, entry):
        self.entry = entry

    def add_item(self, **fields):
        """Record one per-unit outcome (a strain, an isolate). Appended in run order."""
        self.entry.setdefault("items", []).append(fields)
        return self

    def set(self, **fields):
        """Set/overwrite run-level fields (status, summary counts, snapshot path, ...)."""
        self.entry.update(fields)
        return self


def _append(log_path, entry):
    """Append one JSON line. Best-effort: a write failure warns but never raises."""
    try:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(entry, sort_keys=True) + "\n")
    except Exception as e:  # read-only FS, permissions, ... -- degrade, don't break the run
        print(f"runlog: could not append run record to {log_path}: {e}", file=sys.stderr)


@contextlib.contextmanager
def record_run(kind, *, reference_path=None, tools=(), log_dir=None, extra=None):
    """Context manager that logs one run of `kind` to `<log_dir>/<kind>.jsonl`.

    Captures provenance up front (timestamp, code commit + dirty flag, tool versions,
    reference digest, host/python), yields a RunRecord to fill in, and appends the entry
    on exit -- even on exception (status='error', with the error recorded and re-raised).

    kind            "sanity" | "fill" | "recall" (also the log filename stem)
    reference_path  the ERG11 reference CDS FASTA to digest (proves what was called against)
    tools           external binaries whose --version to capture (the reproducibility layer)
    log_dir         override the default data/earlywarning/runlog (tests point elsewhere)
    extra           dict merged into the entry at start (e.g. cli args, snapshot path)
    """
    log_dir = pathlib.Path(log_dir) if log_dir is not None else DEFAULT_LOG_DIR
    sha, dirty = git_commit()
    entry = {
        "schema": SCHEMA,
        "kind": kind,
        "run_id": uuid.uuid4().hex[:12],
        "started_at": _ew.utcnow_iso(),
        "code_commit": sha,
        "code_dirty": dirty,
        "tools": tool_versions(tools),
        "host": platform.node() or None,
        "platform": platform.platform(),
        "python": platform.python_version(),
        "status": "ok",
        "items": [],
    }
    if reference_path is not None:
        entry["reference"] = {
            "name": os.path.basename(str(reference_path)),
            "sha256": file_digest(reference_path),
        }
    if extra:
        entry.update(extra)

    rec = RunRecord(entry)
    try:
        yield rec
    except BaseException as e:  # noqa: BLE001 -- log the crash, then let it propagate
        entry["status"] = "error"
        entry.setdefault("error", str(e).splitlines()[0][:200] if str(e) else type(e).__name__)
        raise
    finally:
        entry["finished_at"] = _ew.utcnow_iso()
        _append(log_dir / f"{kind}.jsonl", entry)


def read_runs(kind, log_dir=None):
    """Read back every run of `kind`, oldest first. Skips any malformed line rather than
    failing -- an append-only log is only as robust as its worst line."""
    log_dir = pathlib.Path(log_dir) if log_dir is not None else DEFAULT_LOG_DIR
    path = log_dir / f"{kind}.jsonl"
    runs = []
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    runs.append(json.loads(line))
                except json.JSONDecodeError:
                    continue
    except FileNotFoundError:
        pass
    return runs
