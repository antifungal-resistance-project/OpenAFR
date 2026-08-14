# Early-warning delivery (issue #26)

The delivery tip of the C. auris azole early-warning track. The scheduled pipeline
(`.github/workflows/earlywarning.yml`) pulls the latest NCBI release, runs
detect → map → estimate → compose, and **delivers only on genuine new news**. This
directory is the git-tracked memory of that delivery.

## Layout

```
data/earlywarning/digest/
  README.md      # this file (committed)
  STATE.json     # committed: the most-urgent tier ever delivered, per mutation
  DIGEST.md      # committed on delivery: the latest full alert brief + what's new
```

## What "new news" means

A run delivers only when it has something a human hasn't already been told:

- a mutation **never delivered** an alert for before, or
- a mutation whose so-what tier has **escalated** above the most urgent tier ever
  delivered for it (context → watch, watch → act-now).

A de-escalation (act-now → watch) is deliberately *not* news — we already raised the
louder alarm. Everything else is a **quiet no-op**: no commit, no issue, no noise. The
per-mutation memory that makes this decision lives in `STATE.json`, so its git history
*is* the delivery log — every escalation is a one-line diff on the run that made it.

## Today's honest state

On real NCBI snapshots every `erg11_call` is still empty (calls are pending the ERG11
re-caller, #23), so the composer returns "no calls to assess" and the schedule
correctly stays silent. It lights up automatically the moment #23 lands a real call.
The machinery is proven now against the synthetic #22 fixture:

```bash
python scripts/deliver_alert.py demo --markdown   # deliver -> re-run no-op -> escalation
```

## Commands

```bash
# Pull the latest release and deliver if there's new news:
python scripts/deliver_alert.py run --as-of 2025-01-01

python scripts/deliver_alert.py run --as-of 2025-01-01 --dry-run   # decide, write nothing
python scripts/deliver_alert.py run --snapshot <snap.tsv> --as-of 2025-01-01
```

Library: `openafr/delivery.py`. Tests: `tests/test_delivery.py`.
