# Runbook — the production FKS1 fill that lights up the weather page (issue #171)

This is the **surveillance-cadence** fill: the one step between a research repo and a live
product. The render path is FKS1-wired and `weather.yml` publishes weekly, but the page shows the
day-0 "watching since…" state because no FKS1 calls exist yet. This runbook produces a dated,
committed **FKS1 call table** (`data/earlywarning/fks1_calls/`) that the public page overlays onto
its live pull — same data, same source, now lit.

It reuses the host setup and tool gate from [`RUNBOOK_fks1_run.md`](RUNBOOK_fks1_run.md) (the
first-prevalence run); read its §0 once. The difference here: the fill is **scoped to the
trailing-180-day recent window the page displays** (not the full ~29k archive), so the billable
GCP cost tracks what is shown.

**Billable cloud compute — needs explicit go-ahead before step 3.** The reads→consensus
toolchain is bioconda linux-64 only and pulls multi-GB SRA downloads per isolate; it cannot run
in a stock GitHub Action. Run on the Google Cloud (or equivalent) Linux/x86 VM from
[`scripts/cloud_recaller_setup.sh`](../scripts/cloud_recaller_setup.sh).

## The known gotcha: disk, not compute

Each isolate's workdir is deleted after it is called, so disk stays flat at ~one isolate's
scratch (~2–3 GB) on a **native x86 host**. The old local-Docker-under-emulation path does *not*
reclaim VM raw-disk across runs and overflows even on the sanity set — that is why this runs on a
real VM, not a laptop. Budget bandwidth + unattended runtime for the SRA downloads, not disk.

## Sequence (on the VM, from the OpenAFR clone)

```bash
# 0. Setup (once) — see RUNBOOK_fks1_run.md §0; must be a linux-64 host with samtools >= 1.13.
bash scripts/cloud_recaller_setup.sh && conda activate openafr
python -m pytest -q                      # tested core green here

# 1. Pull + write the current NCBI release snapshot (non-billable download; big blob gitignored,
#    INDEX.tsv row committed).
python scripts/snapshot_ncbi_auris.py pull        # writes data/earlywarning/snapshots/PDG….tsv

SNAP=data/earlywarning/snapshots/PDG000000067.<N>.tsv   # the file step 1 just wrote

# 2. Preflight — the BOUNDED workload (no tools/network; runnable anywhere). Confirms the fill
#    re-calls only the recent-window pending-with-reads isolates, i.e. exactly what the page shows.
python scripts/recall_fks1.py plan "$SNAP"              # defaults: --window-days 180, as-of today-180
#   -> prints "scope: recent window as-of … 180d (K isolate(s))" and "fill would re-call: M"

# 3. *** BILLABLE *** The windowed GCP fill. Writes the snapshot AND the committed call table.
python scripts/recall_fks1.py fill "$SNAP"              # same window default as plan
#   -> wrote data/earlywarning/fks1_calls/PDG000000067.<N>.tsv (…, K resolved call(s))
#      and appends a provenance row to data/earlywarning/fks1_calls/INDEX.tsv

# 4. Sanity gate — the caller still behaves on the fixtures (detection-of-record unchanged).
python scripts/recaller_sanity_fks1.py

# 5. Verify the page LEAVES day-0, locally, against the committed call table overlaid on a live
#    pull (same path CI runs; auto-discovers the newest table).
python scripts/render_weather.py --as-of "$(date -d '180 days ago' +%F)" --window 180 \
    --out /tmp/index.html
#   -> stderr: "overlaid N FKS1 call(s) from …"; open /tmp/index.html — it shows a real
#      echinocandin picture, not "watching since …".

# 6. Commit ONLY the small durable artifacts (the big snapshot blob stays gitignored):
git add data/earlywarning/fks1_calls/ data/earlywarning/snapshots/INDEX.tsv \
        data/earlywarning/runlog/fill-fks1.jsonl
git commit -m "#171: production FKS1 fill — commit the dated call table (release <N>)"
git push     # the weekly weather.yml render overlays it automatically; no new secrets
```

`--as-of` / `--window-days` default so the fill scope equals the page window; pass `--no-window`
only to re-call the whole pending pool (e.g. for a concordance/prevalence study, not the page).

## Cadence

**Monthly is enough.** NCBI's median arrival lag for *C. auris* is ~97 days, so new in-window
isolates accrue slowly; a weekly fill would mostly re-pay for an unchanged window. The weekly
`weather.yml` keeps *publishing* the latest committed table regardless — only the (billable) fill
is monthly. Each fill is idempotent and re-emits the union of calls-so-far, so a missed month
simply catches up on the next run.

Between fills the page is honest about how current it is: `render_weather.py` passes the newest
`INDEX.tsv` row to `weather.render_page`, which surfaces a coverage/freshness clause —
`echinocandin calls: N of M window isolates re-called (X%) · calls current through <as_of>`, or a
muted *awaiting re-call* when the fill has aged out of the displayed window — and the footer states
this monthly contract. So a lapsed fill reads as *stale*, never as a confident current picture;
that is the signal that it is time to re-run this runbook.
