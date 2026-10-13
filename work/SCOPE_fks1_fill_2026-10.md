# FKS1 fill — scope & cost card for the 2026-10 window (issue #171)

**Status: SCOPED, NOT RUN.** Everything below is the non-billable preflight. The actual
re-call (`scripts/recall_fks1.py fill`) is **billable GCP compute and needs explicit go-ahead**
(see the gate at the bottom). Produced 2026-10-11 on a local host with NCBI network access;
no cloud VM was started.

## Why now

The live weather page has gone dark. The only committed fill (`PDG000000067.709`,
`as_of 2026-04-11`, 300-isolate sample) has fully aged out of the page's trailing-180-day
window. Rendered at today's date the page honestly reads:

> *"no called isolates in the recent window — cannot detect mutation emergence yet"* · QUIET,
> 24 days with nothing new · as of 2026-10-11 · window 180 d

That is the freshness-honesty machinery (#176) doing its job — but the product's whole point is
a **current** echinocandin picture. Only a fresh fill re-lights it.

## The current window (measured, non-billable)

```
# 1. pull  (ran — non-billable)
python3 scripts/snapshot_ncbi_auris.py pull
  -> release PDG000000067.710: 30476 isolates, 95.9% SRA-callable

# 2. plan  (ran — non-billable, no tools/network)
python3 scripts/recall_fks1.py plan data/earlywarning/snapshots/PDG000000067.710.tsv
```

| Quantity | Value |
|---|---|
| Current NCBI release | **PDG000000067.710** (was 709) |
| Recent window | as-of **2026-04-14**, 180 d |
| In-window isolates | **3922** |
| Pending-with-reads (fill would re-call) | **3375** = 3375 SRA downloads |
| Prior committed sample (709) | 300 isolates — now fully aged out |

## Cost & runtime (extrapolated from the runlog)

Per-isolate wall time from `data/earlywarning/runlog/fill-fks1.jsonl` (two real fills):
450-iso run = 235 s/iso; 300-iso run = 195 s/iso → **~215 s/iso**, network-bound (serial SRA
download dominates, so a bigger VM barely helps). VM ref: `e2-standard-2` (~$0.067/hr on-demand,
~$0.020/hr spot, us-central1).

| Scope | `--limit` | Coverage of window | Serial wall time | On-demand $ | Spot $ |
|---|---|---|---|---|---|
| Match prior cadence | `--limit 300 --random` | 7.6% | ~18 h | ~$1.2 | ~$0.4 |
| Mid sample | `--limit 500 --random` | 13% | ~30 h | ~$2.0 | ~$0.6 |
| **Full window** | *(no limit)* | 100% | **~200 h ≈ 8.4 d** | ~$13.5 | ~$4.0 |

Full window single-VM is ~8 days unattended. To get the full picture in ~1 day, **shard** the
3375 run-accessions across ~8 VMs (each a `fill --limit` over a disjoint slice; the call table
is idempotent and merges by `isolate_key`). Spot pricing is tempting but an 8-day unattended run
risks preemption mid-download — prefer on-demand for the long single-VM path, or spot only for
short sharded slices.

### Recommendation

The product promises a real surveillance picture, not a prevalence sample, so the **full window**
is the right target — run it **sharded across ~8 spot VMs (~1 day, ~$4)**, or a single on-demand
VM (~8 days, ~$13.5) if simpler ops are worth more than wall-clock. If you only want to re-light
the page cheaply and defer full coverage, `--limit 500 --random` (~$2, 30 h) refreshes the
precedent at better coverage than the old 300.

## Ready-to-run command block (RUN REQUIRES GO-AHEAD)

On the GCP Linux/x86 VM from `scripts/cloud_recaller_setup.sh` (see `work/RUNBOOK_fks1_run.md` §0
for host setup), from the OpenAFR clone:

```bash
bash scripts/cloud_recaller_setup.sh && conda activate openafr
python -m pytest -q

python scripts/snapshot_ncbi_auris.py pull
SNAP=data/earlywarning/snapshots/PDG000000067.710.tsv   # confirm this is newest

# *** BILLABLE *** choose ONE:
python scripts/recall_fks1.py fill "$SNAP"                      # full window (~8 d / ~$13.5 on-demand)
# python scripts/recall_fks1.py fill "$SNAP" --limit 500 --random   # cheap refresh (~30 h / ~$2)

python scripts/recaller_sanity_fks1.py                         # caller sanity gate
python scripts/render_weather.py --as-of "$(date -d '180 days ago' +%F)" --window 180 \
    --out /tmp/index.html                                      # confirm page leaves day-0

git add data/earlywarning/fks1_calls/ data/earlywarning/snapshots/INDEX.tsv \
        data/earlywarning/runlog/fill-fks1.jsonl
git commit -m "#171: production FKS1 fill — dated call table (release 710)"
git push   # weekly weather.yml overlays it automatically
```

**Sharding variant** (full window, ~8 VMs, ~1 day): on VM *i* of *n*, run
`python scripts/recall_fks1.py plan "$SNAP" --list-runs` to get the 3375 accessions, split into
*n* disjoint accession lists, and give each VM its slice. Each VM commits its own call-table
rows; the table merges idempotently by `isolate_key`.

## What this card did NOT do

- Did **not** run `recall_fks1.py fill` — no SRA downloads, no VM, no spend.
- Only the release-710 `snapshots/INDEX.tsv` provenance row is committed (the multi-MB snapshot
  blob stays gitignored and is re-pulled on the VM). The `709` call table is untouched.
