# ERG11 re-caller run log — append-only provenance of every tool+network run

The durable, git-diffable memory of *what actually ran*. The re-caller's two
un-reproducible paths (they shell out to fasterq-dump/minimap2/samtools and hit the
network) print to the terminal and, for `fill`, overwrite a snapshot in place — so
without this log each run's outcome and provenance vanish when the terminal scrolls, and
a re-`fill` silently erases the last one's calls. Written by `openafr/runlog.py`.

Same convention as `../snapshots/INDEX.tsv`: the big regenerable blobs stay out of git;
this small record of "what happened" is committed. One JSON line per run is tiny.

## Files

```
runlog/
  README.md        # this file
  sanity.jsonl     # one line per `scripts/recaller_sanity.py` run
  fill.jsonl       # one line per `scripts/recall_erg11.py fill` run
  recall.jsonl     # one line per `scripts/recall_erg11.py recall` run
```

Each `.jsonl` is **append-only**: never rewrite past lines. Commit new lines after a run.
The files may not exist until the first run of that kind (the tool+network paths only run
on a Linux/x86 host with the bioconda tools — see `../../../scripts/README_recaller_sanity.md`).

## One run's fields

| field | meaning |
|---|---|
| `schema` | run-log line schema version (currently `1`) |
| `kind` | `sanity` \| `fill` \| `recall` |
| `run_id` | random 12-hex id for this run |
| `started_at` / `finished_at` | UTC ISO-8601 (`utcnow_iso` format) |
| `status` | `ok` \| `fail` \| `error` (`error` = the run body crashed; see `error`) |
| `code_commit` / `code_dirty` | repo HEAD short-SHA and whether the tree had uncommitted changes |
| `reference` | `{name, sha256}` of the ERG11 reference CDS the call was made against |
| `tools` | `{fasterq-dump, minimap2, samtools}` → their `--version` first line (`null` if absent) |
| `host` / `platform` / `python` | where it ran |
| `summary` | run-level counts (kind-specific: sanity → pass/fail/inconclusive; fill → called/partial/failed/attempted) |
| `items` | per-unit outcomes in run order (per strain for sanity; per isolate for fill) |
| *(fill)* `snapshot`, `snapshot_sha256_before/after`, `wrote` | which snapshot, its content hash before and after, and whether the file was actually written (`false` on `--dry-run`) |

`status: "error"` with equal before/after fill digests and `wrote: false` is the honest
signature of a crashed mid-batch run that changed nothing on disk.

## Reading it back

```python
from openafr import runlog
for run in runlog.read_runs("sanity"):      # oldest first, malformed lines skipped
    print(run["started_at"], run["status"], run.get("summary"))
```

`work/RESULTS_recaller_sanity.md` is the human-readable companion: the curated PASS/FAIL
history and how to interpret an `INCONCLUSIVE`.
