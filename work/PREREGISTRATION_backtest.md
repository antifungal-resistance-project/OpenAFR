# Pre-registration — early-warning backtest / honest validation (issue #27)

**Written 2026-08-14, BEFORE computing a single backtest number.** Everything below
is fixed in advance. Deviations must be recorded as deviations, not edited into this
document. This mirrors the discipline of [PREREGISTRATION_run2.md](PREREGISTRATION_run2.md)
and [PREREGISTRATION_auris_Y132F.md](PREREGISTRATION_auris_Y132F.md): the criterion and
the interpretation of *both* outcomes are committed now, so no result can be spun after
the numbers exist.

This is the credibility test for the whole resistance early-warning track (issues
#20–26). Its charter (issue #27): *prove the early-warning claim, or honestly fail it.*

## The claim under test

The track promises an **early warning**: that replaying historical NCBI snapshots, the
pipeline would have flagged an emerging *C. auris* ERG11 azole-resistance mutation
(the canonical example is **Y132F**) **before** clinical / published surveillance
characterised it — with a positive **lead time**.

## The load-bearing constraint, stated before we start

The spike (#20, [RESULTS_earlywarning_spike.md](RESULTS_earlywarning_spike.md)) and the
detector write-up (#22, [RESULTS_emergence.md](RESULTS_emergence.md)) already establish,
on the record, that:

- NCBI runs **no** AMR pipeline on *C. auris*: `erg11_call` is empty for **every** real
  isolate. The resistance call was always ours to manufacture from linked SRA reads via
  an **ERG11 re-caller**.
- **That re-caller was never built.** It appears nowhere in the shipped milestone
  (items 1–8 are spike, snapshot, emergence, mapping, structural, alert, delivery,
  backtest); issue #23 became the pocket *mapping*, not the re-caller. So the
  resistance-call **input** to the entire pipeline does not exist on real data.

Therefore the emergence-claim backtest **on real resistance calls is structurally
impossible today** — not because the method is broken, but because its input is absent.
Pretending otherwise (e.g. defaulting an uncalled isolate to wild-type, or hand-labelling
history) would manufacture the very result we are supposed to be honestly testing. We
refuse to. Instead we pre-register three things we *can* test honestly, and commit to
reporting the blocked half as the blocked half.

## What we test (three parts, each with its outcome pre-committed)

### H1-arrival — the lead-time *budget* (REAL data, the measurable half)

The spike (#20) explicitly deferred to #27 the one arrival-axis number nobody has yet
measured: **the deposit lag** `L_deposit = collection_date → target_creation_date`
(sample taken → isolate visible in NCBI). This bounds *any* achievable lead time: even a
perfect re-caller cannot warn earlier than the data becomes visible. We measure it on the
real latest snapshot.

**Metric.** Over isolates with *both* a parseable `collection_date` and
`target_creation_date`, `L_deposit` in days. Report median, p25/p75/p90, and the fraction
visible within 365 days of collection. Isolates missing either date are **excluded from
the ratio and reported as a count** (never imputed) — the #22 undated discipline.

**Pre-committed criterion (fixed now):**

    PASS-arrival  if  median L_deposit <= 365 days
                  AND fraction visible within 365 days >= 0.60
    FAIL-arrival  otherwise

**Interpretation, committed now (both directions):**

- **PASS-arrival.** NCBI visibility typically trails collection by under a year. A
  lead-time budget of that order is *compatible* with catching an emerging lineage
  ahead of the published record — under the stated assumption below. The arrival axis
  does **not** by itself defeat the early-warning claim.
- **FAIL-arrival.** NCBI visibility itself lags too much (median > 1 year). Then the
  "early" in early-warning is **not defensible on the arrival axis**, regardless of how
  good the re-caller is — a genuine, publishable negative for the track's premise.

**Stated assumption (declared, not measured here).** We compare `L_deposit` against the
*published-surveillance* lag for a novel *C. auris* resistance lineage, which the
literature puts on the order of **1–3 years** (clade characterisation, retrospective
genotype panels, and mutation-panel papers typically trail sample collection by that
much). This is an **assumption**, not a quantity this backtest measures; it is stated so
the reader can substitute their own number. The verdict above is gated only on the
measurable `L_deposit`.

### H2-null — the walk-forward on REAL data (the honest null)

Replay the real snapshot as a walk-forward: step an `as-of` date monthly across the
snapshot's real `target_creation_date` range, running `detect_emergence` at each step
with the frozen defaults (`window_days=180, min_count=3, min_delta=0.05,
backlog_frac=0.5`).

**Pre-committed expectation.** Because `erg11_call` is empty at every timepoint, **every
step must return the honest `note` ("no called isolates … cannot detect yet") and zero
signals.** We commit *in advance* that this is the correct output and NOT a pass of the
early-warning claim: a stable "cannot detect yet" across the whole timeline is the
pipeline correctly refusing to emit a false all-clear. If any step instead produced a
signal from empty calls, that would be a **bug** (a fabricated warning) and a FAIL of the
method's honesty contract.

### H3-mechanism — lead time on a LABELED replay (clearly synthetic)

To locate the blocker — is the negative "method can't produce lead time" or "method is
starved of input"? — we run the *same* walk-forward on a **labeled reconstruction**: a
snapshot whose `erg11_call` is populated on a small, explicitly-synthetic cohort in which
a Y132F-shaped mutation first appears at a known collection date and becomes NCBI-visible
after a realistic deposit lag. A configurable `reference_date` stands in for
"clinical/published recognition."

**Pre-committed criterion (fixed now):**

    PASS-mechanism  if the walk-forward first flags the target mutation at an as-of
                    date STRICTLY BEFORE reference_date  (lead_time_days > 0),
                    using the frozen detector defaults and no per-run tuning.
    FAIL-mechanism  otherwise.

**Interpretation, committed now.** This is a **mechanism demonstration, not evidence for
the real claim** — its dates are illustrative, exactly as `detect_emergence demo` is
synthetic. PASS-mechanism means the pipeline *does* convert a first-appearance into a
dated, positive-lead-time flag once calls exist, so H2-null is attributable to the
**missing re-caller**, not to a broken detector. It grants the track **no** real-data
credibility; only a re-caller applied to historical SRA reads could.

## Pre-committed overall verdict

The track's early-warning claim is reported as **NOT-YET-VALIDATED (blocked on the ERG11
re-caller)** unless *real resistance calls* exist and a real walk-forward shows a positive
lead time. That real validation is out of scope here because the input does not exist.
What #27 delivers is: (a) the real arrival-lead-time budget (H1, PASS or FAIL on its own
terms), (b) the honest real-data null (H2), and (c) the isolated proof that the mechanism
produces lead time when fed calls (H3). An honest "not yet — here is exactly what is
missing and here is the arrival budget it would have to work within" is, per the project
ethos, a valid and publishable outcome.

## Method / reproducibility (fixed)

- Library `openafr/backtest.py`, CLI `scripts/backtest_earlywarning.py`, tests
  `tests/test_backtest.py`. Stdlib only, matching the rest of the track.
- Real snapshot: the latest NCBI *C. auris* release pulled by the store (#21),
  `scripts/snapshot_ncbi_auris.py pull`. The release tag + content SHA-256 of the exact
  snapshot used are recorded in [RESULTS_backtest.md](RESULTS_backtest.md) and in
  `data/earlywarning/snapshots/INDEX.tsv`, so the arrival numbers are rebuildable.
- Detector thresholds are the frozen track defaults above; the backtest introduces **no**
  new tunable that could be reverse-fit to a desired lead time.

## Declared limitations (fixed in advance, true of every outcome)

1. **No real resistance calls exist**, so H2 is a null by construction and H3 is
   synthetic. Neither can, or claims to, validate the early-warning claim on real data.
2. **Deposit lag is a *necessary, not sufficient* condition** for lead time. Even a short
   `L_deposit` leaves the re-caller latency and the published-surveillance-lag assumption
   between us and a real warning; H1 bounds the budget, it does not deliver the warning.
3. **US-centric feed (#20).** Any lead-time statement inherits the feed's ~88% USA skew;
   it is "US early warning + thin international," not global.
4. **The published-surveillance-lag figure (1–3 y) is an assumption, not a measurement**
   in this repo; the H1 verdict is gated only on the measured `L_deposit`.
