# Screening strategy — the library ladder, the two paths, and the validation gate

Companion to [VISION.md](../VISION.md) (what the tool is for) and
[BUSINESS_PLAN.md](BUSINESS_PLAN.md) (how the org turns a hit into a medicine). This doc
answers one tactical question those two don't: **now that the tool works, what do we point
it at — and in what order?** Written to be picked up cold.

Date: 2026-09-05.

## TL;DR

- The first screen (the Drug Repurposing Hub, ~6–7k compounds) was the **smallest rung** on a
  ladder that runs up into the billions. It was the right *first* move — cheap, safe hits, fast
  path to clinic — but it is not where the tool's value lives.
- The tool's real power is **enrichment at scale**: screen a huge library, take the top fraction.
  A 10-molecule shortlist doesn't use that power at all.
- There are **two honest business paths** — *fast-to-clinic on constrained chemistry* (keep
  repurposing) vs. *best-possible molecule on make-on-demand space* (novel-at-scale). Different
  libraries, timelines, and funding stories.
- **Both paths run through the same gate:** one cheap wet-lab run that shows the geometry ranking
  produces a real hit. Until that gate is passed, scaling the library is spending compute and
  credibility we haven't earned.

## Where we are today

The current shortlist ([RESULTS_dossier.md](../work/RESULTS_dossier.md), snapshot 2026-09-03) is
**A=0, B=14, C=2** over 16 molecules — of which only **10 are novel candidates**; the other 6 are
azole reference controls. Of the 10 novel, by the tool's own validated criterion (an aromatic-ring
nitrogen reaching the heme iron): **4 in-mechanism, 1 borderline, 5 out-of-mechanism.** So half the
novel list is ranked for a reason the tool has not validated, and nothing clears automatic-A.

The **validated claim** is narrow and honest: broad **novel-chemotype triage** (AUC ≈ 0.81), *not*
razor-top ranking (the EF@1% figure collapsed under re-sampling — see
[RESULTS_reliability.md](../work/RESULTS_reliability.md)) and *not* within-class potency ranking
(AUC ≈ 0.65 — see [PREPRINT_geometry_ceiling.md](../work/PREPRINT_geometry_ceiling.md)). The tool
is a coarse geometric filter that says "this region of chemical space is worth a look," not a ruler
that measures which of two neighbors is better.

## The library ladder

The Repurposing Hub was chosen for *convenience and clinical safety*, not for being good chemical
space for CYP51. It is a rounding error of what exists. (Sizes are order-of-magnitude, not exact.)

| Tier | Rough size | What it is | Trade-off |
|---|---:|---|---|
| **Repurposing libraries** | ~1.5k–7k | Known drugs with clinical/development history (Drug Repurposing Hub, Prestwick, NCATS, FDA-approved sets) | Safe, fast to reposition — but tiny, constrained chemistry |
| **Buyable screening catalogs** | ~1–10 million | Real, in-stock molecules (Enamine ~4M, ChemBridge, ChemDiv, Life Chemicals) | Far broader/better chemistry — but *not* drugs; no safety data |
| **Make-on-demand virtual space** | ~1–50 billion | Synthesizable-on-order (Enamine REAL, ZINC) — a real molecule, just not made until ordered | Where modern big virtual screens live; the tool's scale sweet spot |
| **De novo / generative** | effectively infinite | Invented structures | **Out of scope** — the tool cannot score potency finely enough to *design* or optimize (see [INTERPRETING_A_SCORE.md](INTERPRETING_A_SCORE.md)) |

The Repurposing Hub is roughly 0.0001% of buyable space. There is enormously more we *can* scan,
and the engine was built for exactly that scale.

## The two business paths

Bigger is not automatically better — it trades safety for chemistry. The fork is real and should be
a deliberate choice, not a drift.

| | **Path A — Fast-reposition** | **Path B — Novel-at-scale** |
|---|---|---|
| Library | Larger repurposing sets (Prestwick, NCATS, full approved-drug catalogs) | Buyable millions → make-on-demand billions (Enamine-scale) |
| Hit is | A known drug, with safety history | A brand-new molecule, no safety data |
| Path to patients | Short — reposition, skip early safety | Long — full from-scratch development |
| Cost/time to clinic | Low–moderate, faster | High, years, licensed out |
| Chemistry quality | Constrained (drugs built for other targets) | Best available — the point of the engine |
| Uses the tool's scale power | Partly | Fully |
| Funding story | "Reposition a safe drug against a WHO-priority pathogen" | "Validated triage engine mines billion-scale space for novel antifungals" |

Neither is wrong. Path A is the cheaper, faster, more repurposing-native story; Path B is where the
tool's distinctive value actually pays off. They can run in sequence (A to validate and get an early
signal, B to build the real pipeline) — but the library choice, and therefore the timeline and the
pitch, should be made on purpose.

## The gate that unlocks both

Neither path is justified until **one cheap wet-lab run** shows the geometry ranking produces a real
hit. This is not bureaucracy — it is what earns the right to spend compute on billions and to ask a
partner (or a funder) to take the engine seriously. See
[BUSINESS_PLAN.md](BUSINESS_PLAN.md) checkpoint 1 and [VISION.md](../VISION.md) "the wet-lab hand-off."

**Scope the validation run honestly:**

- The real experiment is the **in-mechanism subset** — the molecules ranked for the reason the tool
  actually validated (aromatic-N reaching the iron) — with the azole controls as positive controls.
  Treat the out-of-mechanism ranks as exploration, not as the test of the criterion.
- The assay is **MIC on live *C. auris*, including the isogenic resistant strains** (Y132F/K143R).
  Those strains are the ground truth and are held by specific academic labs (see
  [../outreach/TARGETS.md](../outreach/TARGETS.md)) — an argument for the academic/reference-lab
  route over a generic CRO for *this* project.
- Cost is small: roughly the price of the compounds (mostly catalogued, ~$1–5k) plus either a
  collaboration (≈ free, co-authorship) or fee-for-service MIC (~$5–20k). The entry ticket to
  validating the entire thesis is low thousands.

**The result routes the whole project:**

- **Signal → scale the library** (choose Path A or B and commit compute), now with proof.
- **No signal → you have learned the method's real-world ceiling cheaply and honestly**, before a
  company is built on it. A ~$10k answer to a question that otherwise costs years.

## The compute reality (why "just scan everything" isn't free)

Docking cost scales with library size, and the pipeline has already hit a throughput wall (the
exhaustiveness-32 "no free speedup" gotcha). Billion-scale screens are tractable only with **staged
filtering** — a cheap first pass over the full set, expensive re-docking under the frozen protocol
only on survivors. The geometry-triage angle is relatively cheap per pose, which helps, but a
full-rigor scan of millions is a real engineering line item, not a button press. Budget it as its
own workstream when Path B is chosen.

## What this does not change

- **The mission stays broad, the active project stays narrow.** Scaling the library is still one
  target (CYP51), one method. New targets and diagnostics remain later projects (see
  [VISION.md](../VISION.md)).
- **The honesty stays load-bearing.** Every larger screen runs under the same hash-frozen protocol
  and the same self-validating gate. A bigger library does not lower the bar for what counts as a
  trustworthy ranking.
- **De novo generation stays out** until the tool can score potency well enough to optimize against —
  which the geometry ceiling says it currently cannot.
