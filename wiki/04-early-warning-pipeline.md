# The Early-Warning Pipeline (Track 2)

*Audience: you understand the [biology](01-antifungal-resistance-primer.md) and the
[overview](02-project-overview.md). This page walks track 2 down to the level of individual
modules and functions.*

---

## What track 2 is trying to do

Track 1 makes new drugs. Track 2 is the **radar**: it watches public genome data from
hospitals worldwide, tries to catch *C. auris* azole-resistance mutations **as they emerge** —
before they show up characterized in the published surveillance record — and for each one
attaches a **structural so-what**: *if this mutation spreads, does fluconazole still fit the
pocket, and how much worse?*

That structural verdict is what keeps track 2 from being a bare "here's a new mutation" novelty
feed. It **reuses track 1's model of the CYP51 pocket** to say something actionable about the
drugs we already have. Same enemy, same enzyme, shared structural moat.

---

## The honest status you must know before reading the code

The entire chain is **built, tested, and pulls real NCBI data.** But it **cannot yet produce a
validated warning**, for one specific, well-understood reason discovered in the opening spike:

> **NCBI runs no resistance-calling pipeline on *C. auris*.** Its Pathogen Detection feed is a
> *metadata + genome-pointer* feed — it tells you an isolate exists, where and when it was
> collected, and points to its raw sequencing reads. It does **not** tell you which resistance
> mutations that isolate carries. The `AMR_genotypes` field is empty for all ~29,000 isolates.

So the resistance signal — the `erg11_call` for each isolate — has to be **manufactured by this
project** from the raw reads, via an **ERG11 re-caller**. The re-caller's *deterministic core*
is built and tested; its *reads-to-consensus orchestration* exists but has not yet been run at
scale on real data (it needs a Linux/x86 host, bioinformatics tools, and large downloads).

Because the calls are absent on every real isolate today, the pre-registered backtest returns
**NOT-YET-VALIDATED** — and this is *correct behavior*, not a bug: fed absent calls, the
detector refuses to invent a warning; fed synthetic calls, it fires **397 days before** the
reference date. **The machinery works; the gap is the input.** Keep this framing in mind as you
read — most of the "honesty constraints" in the code exist to make sure an *absent* call is
never mistaken for a *reassuring* one.

---

## The pipeline stages

Each stage is a library module in `openafr/` plus a CLI in `scripts/`. Data flows left to
right; each module has a `work/RESULTS_*.md` write-up.

```
   NCBI Pathogen        ┌──────────────┐   ┌──────────────┐   ┌────────────┐   ┌────────────┐   ┌───────────┐   ┌────────────┐
   Detection feed  ───▶ │  snapshot    │─▶ │  emergence   │─▶ │  mapping   │─▶ │ structural │─▶ │  alert    │─▶ │  delivery  │
   (real data)          │ earlywarning │   │  emergence   │   │  mapping   │   │ structural │   │  alert    │   │  delivery  │
                        └──────────────┘   └──────────────┘   └────────────┘   └────────────┘   └───────────┘   └────────────┘
                              │                    ▲                                                                    │
                              │                    │                                                                    ▼
                        erg11_call is EMPTY   ┌────┴─────┐                                                     commits a digest +
                        (NCBI gives none)     │ recaller │  ◀── THE MISSING PIECE: reads → erg11_call          files a GitHub issue
                                              │ recaller │      (core built & tested; real run pending)        only on genuine new news
                                              └──────────┘
                                                    ▲
                              ┌─────────────────────┴──────────────────────┐
                              │  backtest — the credibility test for the    │
                              │  whole track (H1 arrival / H2 null / H3 mechanism) │
                              └─────────────────────────────────────────────┘
```

| Stage | Module | CLI | One-line job |
|---|---|---|---|
| snapshot store | `openafr/earlywarning.py` | `scripts/snapshot_ncbi_auris.py` | pull NCBI → normalized, hash-versioned per-release snapshot + baseline/diff |
| emergence detection | `openafr/emergence.py` | `scripts/detect_emergence.py` | flag first-appearing / rising ERG11 substitutions (the signal NCBI can't give) |
| mutation → pocket mapping | `openafr/mapping.py` | `scripts/map_mutation.py` | where the residue sits in the modeled azole pocket; can it be built? |
| structural so-what | `openafr/structural.py` | `scripts/structural_sowhat.py` | does fluconazole still fit if this spreads — measured, estimated, or honest no-call |
| alert composition | `openafr/alert.py` | `scripts/compose_alert.py` | the one-screen brief a clinician acts on |
| delivery + schedule | `openafr/delivery.py` | `scripts/deliver_alert.py` | deliver only on genuine new news; a GitHub Action runs it |
| backtest / validation | `openafr/backtest.py` | `scripts/backtest_earlywarning.py` | the pre-registered credibility test for the whole track |
| **ERG11 re-caller** | `openafr/recaller.py` | `scripts/recall_erg11.py` | **the missing input:** reads → resistance mutation call |
| run provenance log | `openafr/runlog.py` | (used by the above) | append-only record of the un-reproducible tool/network runs |

---

## Stage 1 — Snapshot store (`openafr/earlywarning.py`)

**Why it exists:** to detect that a mutation is *new* or *rising*, you need memory — you must
know what you'd seen before. This module turns each NCBI pull into a **normalized, versioned,
byte-stable snapshot** and gives the two operations everything downstream relies on.

Key design decisions (all consequences of the honesty philosophy):

- **Stable isolate identity.** NCBI's `target_acc` carries a version suffix (`PDT000585001.1`)
  that bumps when NCBI re-processes the same isolate. The stable key strips it
  (`PDT000585001`) so the same isolate is recognized across releases.
- **Byte-stable per release.** Rows are sorted by isolate key and the pull *timestamp is kept
  out of the data file*. Re-pulling the same NCBI release yields an *identical* file, so a
  `git diff` between releases shows only real change and re-running is idempotent. Provenance
  (timestamp, tool version, source URL, content hash) lives in a separate `INDEX.tsv` manifest.
- **The resistance field is explicitly *pending*, never faked.** Because of the spike finding,
  every isolate's `erg11_call` starts **empty**, and a `resistance_source` column records
  *why* — `pending:sra-recaller` (reads exist, awaiting the re-caller) vs. `pending:no-reads`.
  An uncalled isolate is never defaulted to wild-type.
- **Stdlib only** (`urllib` + `csv` + `json` + `hashlib`). A persistence layer shouldn't drag
  pandas into a repo whose whole dependency surface is rdkit + the docking toolkit.

Functions to know: `fetch_release_rows` / `normalize` (raw NCBI → clean records, dropping
un-keyable rows and *reporting* the count rather than inventing them), `write_snapshot` /
`read_snapshot` / `snapshot_sha256` (the byte-stable store), `baseline_as_of` (isolates NCBI
had made visible on or before a date), `diff_snapshots` (added / removed / revised between two
pulls), and `update_index` (the small git-tracked provenance manifest).

> **Two dates, kept rigorously separate — this matters everywhere downstream.**
> `collection_date` = when the *sample was taken* (biological reality). `target_creation_date`
> = when NCBI *made it visible* (what you could have seen in a prior pull). Emergence is
> measured against visibility, not sample age, and the gap between the two dates is itself a
> signal (see the "backlog" logic below).

## Stage 2 — Emergence detection (`openafr/emergence.py`)

**Why it exists:** this is the core method NCBI does not provide. Given the versioned
snapshots, it answers: *compared with what we'd seen before, which ERG11 substitutions are new,
and which are climbing?*

Two signals, defined honestly:

- **FIRST-APPEARANCE** — a substitution present among recently-visible isolates that was never
  seen in the baseline window.
- **RISING** — a substitution whose share of *called* isolates in the recent window exceeds its
  baseline share by at least a threshold (`min_delta`).

Both require a minimum absolute count (`min_count`, default 3) so a single sequenced isolate
can't raise an alarm.

The honesty machinery inside it:

- **The denominator is *called* isolates only.** An isolate with an empty `erg11_call` is
  **not** counted as susceptible — proportions are taken over called isolates, and the calling
  coverage is reported alongside every verdict. So "6% carry Y132F" can never be silently read
  off 2%-called data. On today's real snapshots coverage is 0%, so `detect_emergence` returns
  no signals and *says so* — an honest "cannot detect yet," not a false "all clear."
- **The backlog guard (`assess_backlog`).** A lab depositing a years-old backlog makes *old*
  isolates appear *now* — a spike in the recent window that's really archival, not real-time
  emergence. Using the two-dates discipline, each flagged mutation is checked: if most of its
  recent carriers were *collected* long before they became *visible*, the signal is annotated
  `possible_backlog` with the reason — rather than dropped or trusted blindly.
- `parse_call("Y132F,K143R") → {"Y132F", "K143R"}` — tolerant parsing that never invents a
  call and never silently drops a nonempty one.

Main entry point: `detect_emergence(records, as_of, …)` returns a fully self-describing dict —
the thresholds used, the coverage numbers, and the ranked list of flagged signals — so a
verdict explains itself.

## Stage 3 — Mutation → pocket mapping (`openafr/mapping.py`)

**Why it exists:** a flagged token like "Y132F is rising" is not yet a *structural* statement.
This module answers: *given a flagged substitution, where is that residue in the modeled azole
pocket, does it touch the drug, and can we build the mutant and dock it?* **This is where the
moat starts — no other genomic surveillance feed does this.**

It reuses the exact receptor track 1 docks against (`work/receptor_A.pdb`, the 5TZ1 chain-A
CYP51 scaffold) — same heme, same iron, same coordinate frame — so a mapping verdict lives in
the same structural universe as every gate run.

How the pocket is defined — **measured, not asserted:** a residue "lines the pocket" if any of
its heavy atoms comes within a contact cutoff (default 4.5 Å, a standard van-der-Waals contact
distance) of any atom of the co-crystallized azole (`work/ref_VT1.pdb`). So "Y132 contacts the
drug" is *computed from the structure*, not claimed.

Three honesty constraints (a recurring shape across track 2):

1. **Assert the wild-type, don't assume it.** The token's wild-type letter must match the
   residue actually present at that position in the model. If residue 132 isn't Tyrosine, the
   mapping is *refused* with the observed residue reported — a numbering/structure mismatch is
   surfaced, never papered over.
2. **Non-residue tokens are refused, not forced.** Something like `TR34` (a promoter/tandem-
   repeat call, not a point substitution) cannot map to a single residue, and the module says
   so rather than picking a residue for it.
3. **Buildability is honest about the environment.** The pinned environment has no side-chain
   *repacker*, so only pure atom-**deletion** mutations can be built deterministically — Y132F
   (Phe = Tyr minus a hydroxyl) and V125A (Ala = Val minus two methyls). A mutation that *adds*
   atoms (K143R, F126L) is **mapped** (we know where it is and whether it touches the drug) but
   flagged as needing a repacker to *build* — mapping a residue and building its mutant side
   chain are two different claims, kept separate.

Main entry point: `map_mutation(token, residues, ligand_atoms, fe, …)` → a self-describing
verdict (mappable? distance to drug, distance to iron, pocket contact?, buildable?).

## Stage 4 — Structural so-what (`openafr/structural.py`)

**Why it exists:** the mapper says *where*; the clinician asks *so what — does fluconazole
still fit, and how much worse?* This module turns a mapping verdict into a fit-consequence
verdict. It is deliberately narrow about what it will claim, because **a wrong "the drug still
works" is the single most dangerous output the whole track could emit.**

Three evidence tiers, strongest to weakest — and every branch resolves to an explicit
confidence:

1. **EMPIRICAL-DOCK (`confidence: measured`).** For a mutation actually run through track 1's
   *C. auris* port pipeline (build the mutant, dock the held-out azoles, grade the geometric
   gate), it reports the *measured, pre-registered* outcome — not a re-run. Right now that is
   **Y132F alone** (`EMPIRICAL_DOCKS`). Even here it ships every caveat the port declared: the
   measured finding is that broad iron-reach *survives* (azoles as a class still coordinate the
   iron) while the top-rank *advantage* collapses — and it flags that the port docked an azole
   *class*, not fluconazole specifically.
2. **GEOMETRIC-ESTIMATE (`confidence: estimate`).** For a deletion-buildable mutation not yet
   docked, it gives a *direction, never a number*: `lost_pocket_contacts` deterministically
   measures which removed atoms were in contact with the co-crystal drug — remove contacting
   atoms → predicted *weaker* fit; remove only second-shell atoms → minimal *direct* effect.
   The magnitude waits on a dock, and it prints the exact command to get it.
3. **NO-CALL (`confidence: none`).** For a mutation the environment can't even build (an
   add-atom substitution needing a repacker, like K143R/F126L) or a token that doesn't map at
   all (TR34), it **refuses to call the fit.** A pocket-contact-but-unbuildable mutation is
   flagged "plausible but unresolved"; a second-shell one, "low prior." It never invents a
   ΔΔG it cannot compute.

The honesty bar in one line: *"the drug still fits" is the claim it is most reluctant to make,
and "no confident call" is a first-class, expected outcome.*

## Stage 5 — Alert composition (`openafr/alert.py`)

**Why it exists:** the three stages before it each answer one question and stop. None of them
is an *alert* — the assembled, human-readable statement a clinician or public-health analyst
can act on without opening a single PDB file:

> *"F126L is newly appearing on 4 isolates across 2 countries, first collected 2024; it lines
> the drug-contact surface and predicts weaker azole fit (estimate)."*

`compose_alerts` runs the whole detect → map → estimate loop, joins each emergence signal to
its structural verdict, and enriches it with the epidemiology the detector dropped: which
regions, first-seen date (computed over *all* carriers, not just the recent window), total
carriers. It reports first-*collected* and first-*visible* as **separate** fields — the
two-dates discipline again — so the reader can see the backlog gap rather than have it averaged
away.

It assigns a **so-what priority**, derived from the two upstream verdicts, never invented:

- **ACT-NOW** — a genuine (non-backlog) signal whose structural verdict predicts a *weaker*
  azole fit (measured or estimated). The dangerous combination: real spread **and** a
  defensible drug-fit consequence.
- **WATCH** — a genuine signal whose fit consequence is minimal, uncertain, or a no-call. Real
  spread, but no defensible fit story yet.
- **CONTEXT** — a possible-backlog signal. Reported in full (never hidden) but de-prioritized,
  because the backlog guard says the spike is likely archival.

Priority is a *triage aid, not a clinical claim* — every alert still carries the raw emergence
numbers and the structural caveats so the reader can overrule it. Renders two ways:
`render_markdown` (for a human) and `render_json` (for a machine).

## Stage 6 — Delivery + schedule (`openafr/delivery.py`)

**Why it exists:** an alert nobody sees is not an early warning. This is the layer that makes
the alert *find a human, on a cadence, minimally.* The whole point is **the quiet no-op**:

> A watch that cries every week is a watch nobody reads. The vast majority of runs surface the
> same endemic mutations at the same tier as last time — and those runs must produce
> **nothing**: no commit, no issue, no noise.

It delivers only on **genuine news**, defined narrowly: a mutation never alerted before, *or* a
mutation whose so-what tier has **escalated** above the most-urgent tier ever delivered for it
(watch → act-now, etc.). A *de-escalation* is deliberately **not** news — you already raised
the louder alarm; walking it back quietly is correct. This "most-urgent-ever" memory lives in a
tiny committed `STATE.json`, so the decision is reproducible, git-diffable, and a fresh checkout
picks up where the last delivery left off — *the git history of that file is the delivery log.*
When there *is* news, it renders a digest (committed to the repo) and a GitHub-issue body.
`.github/workflows/earlywarning.yml` runs the whole thing on a schedule.

## Stage 7 — Backtest / validation (`openafr/backtest.py`)

**Why it exists:** this is the credibility test for the *whole* track. The promise is an *early
warning* — replayed over history, it would flag an emerging mutation (canonically Y132F) with a
positive lead time. Testing that on real data needs the one thing the track hasn't run yet: the
ERG11 re-caller. So this module **refuses to fake the missing input** and instead measures the
three things that *are* honest (all pre-registered in `PREREGISTRATION_backtest.md`):

- **H1-arrival** — the lead-time *budget*: the deposit lag from `collection_date` →
  `target_creation_date` over **real** isolates. Even a perfect re-caller can't warn before the
  data is visible, so this bounds any achievable lead time. Pre-registered pass bar: median lag
  ≤ 365 days **and** ≥ 60% visible within 365 days. This one **passes on real data** (median
  ~97-day lag) — which is what justifies building the re-caller at all.
- **H2-null** — a walk-forward over the **real** snapshots: every step returns the honest empty
  state, because calls are absent. **This is the correct output, not a failure** — the pipeline
  refusing to emit a false all-clear.
- **H3-mechanism** — the *same* walk-forward over a labelled **synthetic** replay, proving the
  detector converts a first-appearance into a dated, positive-lead-time flag *once calls
  exist*. This isolates H2's emptiness to the missing re-caller, not a broken method.

It also computes the track's eventual deliverable number, the **azole event frequency**
(`panel_prevalence`): once the re-caller has populated calls, what fraction of the isolates it
could actually resolve carry a known panel substitution — reported with a **Wilson score
confidence interval** (chosen over the textbook Wald interval precisely because the true
proportion is near an extreme and the first samples are small, where Wald misbehaves). Until
`fill` populates calls it is honestly "NOT MEASURED YET." `wilson_halfwidth_worst_case` even
lets you size the run *before* peeking at the data — the same decide-n-in-advance discipline as
the pre-registrations.

---

## The missing piece — the ERG11 re-caller (`openafr/recaller.py`)

This is **the single blocker for the whole track** and the top item in [TODOS.md](../TODOS.md).
It's the thing that would turn every empty `erg11_call` into a real resistance call.

**The boundary it draws.** Turning raw reads into a per-isolate ERG11 *consensus coding
sequence* is the messy, tool-and-network half (`fasterq-dump` / `minimap2` / `samtools`, in
`scripts/recall_erg11.py`, gated on those binaries exactly like the docking scripts gate on
`vina`). `openafr/recaller.py` starts *after* that: **a consensus CDS in, substitution tokens
out.** Everything on this side is **stdlib-only and deterministic**, so it's unit-tested against
the pinned reference; the un-reproducible tool-orchestration side is kept out of the tested core
and honestly labelled.

Its four honesty constraints (the same shape as the other modules):

1. **Assert the reference, don't assume it.** Every call is relative to the pinned B8441 ERG11
   reference (accession XM_085597798.1, SHA-256 verified at load). `load_reference` also
   **refuses** a reference whose hotspot residues don't translate to the panel's wild-type
   letters (V125/F126/Y132/K143) — a numbering or sequence swap is surfaced, never silently
   mis-called. *This is what makes an emitted "Y132F" trustworthy: the code has proven residue
   132 really is a Tyrosine in its reference.*
2. **An uncalled codon is uncalled, never "wild type."** A codon containing an ambiguous or gap
   base can't be translated to a confident residue, so that position is reported *uncalled* — it
   never contributes a reassuring "no change." An isolate that couldn't be sequenced across a
   hotspot says so.
3. **Frame is asserted, not guessed.** The consensus must be an in-frame CDS the same length as
   the reference. A length that isn't a multiple of 3, or that disagrees with the reference
   length (a frameshift/indel), is *refused* with the reason — never forced into a frame.
4. **Tag the panel, don't invent novelty.** Every non-synonymous difference from the reference
   is emitted as a token. Known panel positions are *tagged* so the alert can say "this is a
   named resistance mutation," but non-panel substitutions (clade-defining polymorphisms,
   variants of unknown significance) are emitted verbatim — never suppressed, never asserted to
   be resistance. Deciding which token is *novel/rising* is the emergence detector's job, not
   the re-caller's.

Key functions: `translate` (standard genetic code; an ambiguous codon → `X`, i.e. uncalled),
`load_reference` (hash + numbering check), `call_substitutions` (consensus → tokens + panel
hits + uncalled positions), and `format_call` (renders the result into the `(erg11_call,
resistance_source)` pair that drops straight into a snapshot row — where an empty call is never
ambiguous between "sequenced, wild type" and "could not sequence").

## Provenance logging (`openafr/runlog.py`)

The re-caller's tool+network paths (`recaller_sanity.py`, `recall_erg11.py {recall,fill}`)
either only print or mutate a snapshot **in place**. Without a log, each run's outcome vanishes
when the terminal scrolls, and a re-`fill` silently overwrites the last one's calls. This module
records one JSON object per run to an append-only `.jsonl` — capturing code commit + dirty flag,
external tool versions, the reference digest, host/python — so months later you can answer *what
did that run do, on which inputs, with which tools, against which code?* Its design contract:
**logging must never fail or alter the science run.** Every provenance probe is best-effort
(returns `None` on any error rather than raising); `record_run` is a context manager that writes
the entry on exit *even if the run body crashed*, and re-raises. Same spirit as the snapshot
`INDEX.tsv`: regenerable blobs stay out of git, the small durable "what happened" is committed.

---

## Try it without NCBI (synthetic demo mode)

Every stage has a `demo` mode so a fresh clone can watch the chain fire end to end without
pulling reads or waiting on the missing re-caller:

```bash
python scripts/detect_emergence.py demo
python scripts/structural_sowhat.py estimate Y132F K143R F126L V125A TR34
python scripts/deliver_alert.py demo --markdown        # deliver → re-run no-op → escalation
python scripts/backtest_earlywarning.py mechanism      # detector produces lead time once fed calls
```

For the real run order (pull a fresh snapshot → offline preflight → sanity control → `fill` in
small batches), see [TODOS.md](../TODOS.md) and `work/RUNBOOK_recaller_run.md`.

---

## Two caveats the track says out loud

- **Coverage is US-centric** (~88% USA). "Global early warning" is honestly "US early warning +
  thin international." The alert copy reflects that.
- **The whole track is NOT-YET-VALIDATED** until the ERG11 re-caller runs on real data. The
  scaffold is real, the gap is named precisely, and nothing claims validation it hasn't earned.

---

Next: the **[Codebase Reference →](05-codebase-reference.md)** for the full file map, or the
**[Glossary →](06-glossary.md)** if a term tripped you up.
