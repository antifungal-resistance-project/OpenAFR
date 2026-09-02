# Early-warning delivery (issue #26)

The delivery tip of the C. auris early-warning track. The scheduled pipeline
(`.github/workflows/earlywarning.yml`) pulls the latest NCBI release, runs
detect → map → estimate → compose, and **delivers only on genuine new news**. This
directory is the git-tracked memory of that delivery.

## Layout

```
data/earlywarning/digest/
  README.md          # this file (committed)
  STATE.json         # committed: most-urgent tier ever delivered, per mutation (ERG11)
  DIGEST.md          # committed on delivery: latest ERG11 brief + what's new
  STATE_FKS1.json    # committed: same memory for the FKS1/echinocandin track (v2)
  DIGEST_FKS1.md     # committed on delivery: latest FKS1 (detection-only) brief
```

The two genes deliver **independently** and keep **separate** state, so an FKS1 alert
never overwrites the ERG11 memory (or vice-versa). FKS1 is detection-only: its brief
carries no structural fit verdict and never reaches the ACT-NOW tier (echinocandins
coordinate no metal, so the CYP51/heme-iron geometry does not transfer).

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

On real NCBI snapshots every `erg11_call` / `fks1_call` is still empty (calls are
pending the respective re-caller), so the composer returns "no calls to assess" and the
schedule correctly stays silent. It lights up automatically the moment a fill lands a
real call. The machinery is proven now against synthetic fixtures:

```bash
python scripts/deliver_alert.py demo --markdown             # ERG11: deliver -> no-op -> escalation
python scripts/deliver_alert.py demo --gene fks1 --markdown # FKS1 (detection-only)
```

## Commands

```bash
# Pull the latest release and deliver if there's new news:
python scripts/deliver_alert.py run --as-of 2025-01-01              # ERG11 (default)
python scripts/deliver_alert.py run --gene fks1 --as-of 2025-01-01  # FKS1 / echinocandin

python scripts/deliver_alert.py run --as-of 2025-01-01 --dry-run   # decide, write nothing
python scripts/deliver_alert.py run --snapshot <snap.tsv> --as-of 2025-01-01
```

Library: `openafr/delivery.py`. Tests: `tests/test_delivery.py`, `tests/test_fks1_delivery.py`.
