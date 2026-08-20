# OpenAFR Wiki

Welcome. This wiki is written to be picked up **cold** — by someone who knows nothing
about antifungal resistance, nothing about this codebase, or both. Read the pages in
order for a full onboarding, or jump straight to the one you need.

The repository's own [README](../README.md) tells you *what the pipeline does and how to
run it*, and [VISION.md](../VISION.md) covers *strategy and what's genuinely novel*. This
wiki sits underneath both: it **teaches the science and walks through the code** so you
can actually contribute.

## Start here

| # | Page | For you if… |
|---|------|-------------|
| 1 | [Antifungal Resistance — a primer](01-antifungal-resistance-primer.md) | You've never heard of *C. auris*, CYP51, or azoles and want the biology from zero. **No code.** |
| 2 | [Project Overview](02-project-overview.md) | You want the 10-minute map: what OpenAFR is, its two tracks, and the honesty philosophy that shapes every file. |
| 3 | [The Drug-Discovery Pipeline](03-drug-discovery-pipeline.md) | You want a deep, code-level walkthrough of track 1 — docking, the geometry criterion, and the self-validating gate. |
| 4 | [The Early-Warning Pipeline](04-early-warning-pipeline.md) | You want a deep, code-level walkthrough of track 2 — genomic surveillance for emerging resistance mutations. |
| 5 | [Codebase Reference](05-codebase-reference.md) | You want a module-by-module and script-by-script map to find the right file fast. |
| 6 | [Glossary](06-glossary.md) | You hit a term (heme, EF@1%, PDBQT, TR34…) and want a one-line definition. |
| 7 | [Contributing](07-contributing.md) | You're ready to set up the environment, run the tests, and pick a first task. |

## The one idea to take away first

Most of what looks unusual in this codebase — hash-frozen files, a "gate" that can make
the pipeline *fail on purpose*, honesty constraints repeated in every module — exists to
serve a single rule:

> **A computational prediction is a hypothesis, never a result, until it has proven it
> can re-discover something we already know to be true.**

If you understand *why* that rule matters (the [primer](01-antifungal-resistance-primer.md)
and [overview](02-project-overview.md) explain it), the rest of the code stops looking
paranoid and starts looking like the point.

## A note on honesty

This project deliberately documents its own failures and limitations in the same places it
documents its successes. If you find a page here that overclaims — that says something
"works" when the code only says it *might* — that's a bug in the wiki. Please fix it or
flag it. The credibility of the whole project rests on that discipline.
