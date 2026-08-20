# Project Overview

*Audience: you've read (or already know) the [biology primer](01-antifungal-resistance-primer.md)
and want the map of what OpenAFR actually is, before diving into code.*

---

## What OpenAFR is, in one paragraph

OpenAFR is an **open, self-validating computational pipeline** aimed at drug-resistant
*Candida auris* and its azole-target enzyme CYP51/ERG11. It has two tracks against that same
enemy: **(1) drug discovery** — a triage tool that ranks candidate molecules by a
mechanistic geometric criterion instead of trusting a docking program's score, and proves
that ranking is trustworthy before using it; and **(2) genomic early-warning surveillance** —
a pipeline that watches public genome databases for emerging resistance mutations and attaches
a structural "does our drug still work?" verdict to each one. The unifying thread is a strict
discipline of **honesty**: nothing here claims to be validated until it has earned it, and the
code is built to make cheating hard and self-deception visible.

---

## The founding rule

Everything unusual about this repository descends from one rule, stated in the main README:

> **The pipeline must prove it can re-discover the drugs we already know work, on molecules
> it has never seen, before it is allowed to rank anything new. If it cannot, the validation
> gate exits nonzero and blocks the screen.**

Why is this rule necessary? Because computational drug screening is *notorious* for
producing confident-looking rankings that mean nothing. Docking scores are noisy; published
hit rates are low. Anyone can point a docking program at a library, sort by score, and
publish a ranked list — with zero evidence the ranking is better than random. This project
treats *that evidence* as the actual deliverable. The molecules are secondary; the proof that
the method works is the product.

This is why you'll see, throughout the code:

- **A "gate"** — a program (`validate_gate2.py`) that can deliberately **fail** (`exit 1`),
  blocking the pipeline, when the method can't re-find known drugs.
- **Pre-registration** — the rules and pass/fail thresholds are written down and
  cryptographically frozen (SHA-256 hashed) *before* any data is generated, so they can't be
  quietly moved to make a result pass.
- **Honesty constraints** repeated in module after module — an uncalled thing is never
  assumed to be the convenient default; a limitation is stated in the same breath as a result.

If those design choices seem excessive at first, that's the sign you've correctly understood
how easy it is to fool yourself in this field.

---

## The two tracks

```
                        C. auris  +  CYP51 / ERG11
                        (same enemy, same enzyme)
                          /                    \
        ┌────────────────┘                      └────────────────┐
        │                                                        │
  TRACK 1: DRUG DISCOVERY                    TRACK 2: EARLY-WARNING SURVEILLANCE
  "Find new weapons."                        "Watch the enemy evolve."
                                             
  Rank candidate molecules by how            Watch NCBI genome data for emerging
  closely they bring a nitrogen to           azole-resistance mutations, and for
  the heme iron — a mechanistic,             each one answer: does fluconazole still
  geometric criterion — not by the           fit the pocket if this spreads?
  docking score. Prove the ranking
  works on held-out known drugs              Reuses Track 1's structural model of
  before trusting it.                        the pocket to give each mutation a
                                             "structural so-what."
  STATUS: validated core, works.             STATUS: full scaffold built & tested,
  Needs a wet-lab partner (the one           blocked on one missing piece (the
  critical open dependency).                 ERG11 re-caller). Honestly NOT-YET-
                                             VALIDATED until that exists.
```

### Track 1 — Drug discovery (the validated one)

The distinctive move: **it does not trust the docking program's affinity score.** Instead it
ranks molecules by a *geometric, mechanistic* criterion — the **N–Fe distance**, how closely
a ligand's nitrogen reaches the heme iron, which is the literal mechanism by which azoles
kill (see the [primer](01-antifungal-resistance-primer.md#3-how-the-main-antifungal-drugs-work-the-azoles)).

The headline result, from a pre-registered blind test:

| Ranking method | AUC | Verdict |
|---|---|---|
| **N–Fe distance** (the geometric criterion) | **0.794** | **PASS** — meaningfully better than random |
| AutoDock Vina's own affinity score | 0.471 | worse than random, on the *identical* poses |

(**AUC** = a standard 0–1 measure of ranking quality; 0.5 is a coin flip, 1.0 is perfect.
See the [glossary](06-glossary.md).) The finding in one line: **geometry beats the scoring
function** on this system. A permutation test put the odds of that happening by chance at
p ≈ 0.0028.

Track 1 carries one honest, self-inflicted caveat: a *later* pre-registered reliability check
found that one flashy sub-metric (EF@1%, the razor-top-1% enrichment) was inflated by docking
each molecule only once, and collapses under repeated sampling. The **broad** ranking power
(the AUC) survives; the razor-top ranking does not. So the honest headline is *"a trustworthy
top-~5% shortlist, not a razor-top ranking."* The fact that the project **found and published
its own inflated number** is exactly the discipline the whole thing is about.

→ Full code-level walkthrough: **[The Drug-Discovery Pipeline](03-drug-discovery-pipeline.md)**

### Track 2 — Early-warning surveillance (built, honestly blocked)

A second pipeline watches **NCBI Pathogen Detection** (a public feed of pathogen genome
metadata) for *C. auris* isolates, and tries to flag *emerging* azole-resistance mutations
early — then attaches a structural verdict: *if this mutation spreads, does fluconazole still
fit the pocket?* Track 1 finds new drugs; track 2 watches resistance evolve against the drugs
we already have.

The whole chain — snapshot → emergence detection → mutation-to-pocket mapping → structural
so-what → alert → scheduled delivery → backtest — is **built, tested, and pulls real NCBI
data.** But it cannot yet produce a *validated* warning, and the docs say so plainly. The
load-bearing discovery:

> **NCBI runs no resistance-calling pipeline on *C. auris*.** Its feed tells you an isolate
> *exists* and points to its raw sequencing reads, but it carries **no resistance call** — the
> `AMR_genotypes` field is empty for all ~29,000 isolates.

So the resistance signal has to be **manufactured by this project**, from the raw reads, via
an **ERG11 re-caller** (reads → resistance mutation call). The *deterministic core* of that
re-caller is built and tested; the messy reads-to-consensus orchestration exists but has not
yet been run at scale on real data. That single missing run is the one blocker for the whole
track, and the pipeline is honest about it: fed absent calls, the detector correctly refuses
to invent a warning; fed synthetic calls, it fires hundreds of days early — proving the
machinery works and **the gap is the input, not the method.**

→ Full code-level walkthrough: **[The Early-Warning Pipeline](04-early-warning-pipeline.md)**

---

## What's genuinely novel (told honestly)

The [VISION](../VISION.md) is blunt about this, and so is this wiki:

- **The ingredients are not new.** Rescoring docking poses by a mechanistic/geometric
  constraint (metal-coordination, pharmacophore constraints) is a known, decades-old idea.
  This project is not the first to distrust a docking score.
- **The discipline and the specific proven result are the contribution.** A pre-registered,
  hash-frozen, blind-holdout, permutation-tested screen — with a gate that documents a run
  that *failed first* — is rare, especially in open/amateur docking repositories that publish
  ranked lists with no evidence at all.

Frame it as *"a known idea, actually proven, honestly, on a target that matters"* — not as a
novel algorithm. **The honesty is the edge.**

Two parts of the method are broadly reusable beyond antifungals:

1. **Mechanism-anchored geometric rescoring** — works for any target with a knowable
   structural mechanism (metalloenzymes, clear pharmacophores).
2. **The self-validating gate** — fully domain-independent: *"don't trust a screening ranking
   until it re-discovers what you already know."* This applies to essentially any virtual
   screen and is arguably the most portable idea in the repo.

---

## Repository layout at a glance

```
README.md            what the pipeline does + how to run it (start here for running)
VISION.md            strategy, novelty, where it goes next
TODOS.md             deferred work, each item written to be picked up cold
protocol.yaml        THE frozen run protocol (docking box, search effort, pass bar) — hash-pinned
environment.yml      version-pinned conda toolchain (rdkit, openbabel, vina, pytest)

openafr/             the shared Python library — the load-bearing, tested core
scripts/             the runnable CLIs and shell scripts (the "verbs" of each pipeline)
tests/               known-answer tests locking the math, parsers, and orchestration logic
data/                inputs: structures, ligands, and the early-warning snapshot store
work/                run outputs + every write-up (RESULTS_*.md) and pre-registration (PREREGISTRATION_*.md)
wiki/                you are here
```

The split to internalize: **`openafr/` is the deterministic, unit-tested library; `scripts/`
holds the entry points that call it** (and, for the parts that shell out to external tools or
the network, the un-reproducible glue that's deliberately kept *out* of the tested core). See
the [Codebase Reference](05-codebase-reference.md) for a file-by-file map.

---

## The recurring design patterns (learn these once)

These patterns show up in *both* tracks. Recognizing them makes every file easier to read.

1. **Hash-freezing / pre-registration.** Files that define "what counts as success" —
   `protocol.yaml`, the `PREREGISTRATION_*.md` docs, the ERG11 reference sequence — have their
   SHA-256 hash pinned in code. At run time the code re-checks the hash and loudly flags any
   result as "not credible" if the file was edited after freezing. This is anti-cheating,
   machine-enforced.

2. **"Uncalled is never the convenient default."** When the pipeline doesn't know something,
   it records *that it doesn't know* — it never silently fills in the flattering answer. An
   isolate whose resistance wasn't sequenced is **uncalled**, not "assumed susceptible." A
   molecule that produced no usable docking pose is ranked **last**, never dropped (dropping it
   would inflate every metric). This single principle prevents most ways of accidentally
   faking a good result.

3. **Separate "we measured it" from "we estimate it" from "we can't say."** Especially in
   track 2's structural stage, every verdict is tagged with an explicit confidence:
   `measured` (a real, pre-registered dock exists) vs. `estimate` (a direction, no magnitude)
   vs. `none` (an honest no-call). "No confident call" is a first-class, expected outcome, not
   a failure.

4. **One source of truth, read by everything.** The docking box parameters live in exactly
   one place (`protocol.yaml`) and are read by *both* the docking script and the grading gate,
   so what gets docked can never silently drift from what gets graded. Shared parsers live in
   exactly one place (`openafr/pdbqt.py`) instead of being copy-pasted across scripts.

5. **The big regenerable blobs stay out of git; the small durable memory is committed.**
   Multi-megabyte snapshot files and receptor files are gitignored and rebuilt on demand; the
   tiny provenance records of *what happened* (an `INDEX.tsv`, a run-log `.jsonl`, a
   `STATE.json`) are committed, so the git history itself is the log.

---

## Current status & the critical open dependencies

- **Track 1:** the validated core is complete and passing. The full novel-screening path
  (filter → prep → screen → rank → report) works end to end. **The one blocking dependency is
  a wet-lab collaborator** — nothing computational is validated as a *drug* until someone puts
  candidates on real fungus. A high rank is a hypothesis for a wet lab, never a hit.

- **Track 2:** the whole scaffold is built and tested. **The one blocker is running the ERG11
  re-caller on real data at scale** (it needs a Linux/x86 host with bioinformatics tools and
  large genome downloads). Until then the track honestly reports NOT-YET-VALIDATED. See
  [TODOS.md](../TODOS.md) for the exact run order.

---

## License, in brief

OpenAFR is released under the **PolyForm Noncommercial License 1.0.0** — *source-available*,
not OSI open-source. Free for noncommercial use (academic and public research, nonprofits,
education, health/government, personal). Commercial use requires a separate license from the
foundation. See [LICENSE](../LICENSE) and [Contributing](07-contributing.md).

---

Next, pick your deep dive:
- **[The Drug-Discovery Pipeline →](03-drug-discovery-pipeline.md)** (track 1)
- **[The Early-Warning Pipeline →](04-early-warning-pipeline.md)** (track 2)
