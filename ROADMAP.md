# OpenAFR — Roadmap

The **ranked, living** list of what to work on next. Not a sprint tracker: 4–6 bets, each a
move toward the north star, re-ranked as work lands. Pick the top live bet, ship one PR,
re-rank. For *why* the project exists see [VISION.md](VISION.md); for the archive of
deferred/finished work with cold-pickup context see [TODOS.md](TODOS.md).

**North star (VISION.md, 2026-10-06):** a *live*, auto-updating **OpenAFR Resistance Weather**
page — emerging echinocandin (FKS1) resistance in *C. auris*, backed by the benchmark-certified
FKS1 caller. The page is live (release 710). The work now is making it deeper and keeping it current.

_Last re-ranked: 2026-10-26._

**Operational board:** the [OpenAFR GitHub Project](https://github.com/users/jajjer/projects/1)
is the live tracker — each active bet is an issue with a `Rank` and `Status` (Backlog / Next /
In progress / Parked / Done). Rank and status live on the board; this file is the narrative
*why* and is refreshed when a bet lands. Keep the two in sync: re-rank on the board, re-tell the
story here.

## Active bets (ranked)

1. **Weather depth — keep extending what the page tells you.** ([#222](https://github.com/antifungal-resistance-project/OpenAFR/issues/222))
   The winning thread: each
   recent PR adds a layer (#215 panel breakdown, #216 geo prevalence, #217 plain language,
   #218 claims-not-made, #219 QA, then **marker × region** — each geo-prevalence hotspot now
   names which FKS1 marker drives it; the data shows Mid-Atlantic is clonal D642Y and
   Guangdong is M690I, i.e. independent emergence; schema v5→**v6**). Next increments are
   judgment calls — pick the one that most changes what a surveillance reader learns:
   - **Fix the per-release trend honesty.** The "Resistance trend" section plots 709
     (300-isolate pilot) vs 710 (3375 full) as "prevalence up 0.8 pts since the previous
     batch" — but both cover the ~same April window, so that delta is a sampling-depth
     artifact, not emergence. Either suppress the delta when two fills share a window, or
     label batches by sampling depth rather than time. (Good next pick: a real honesty fix.)
   - Finer geo resolution; clearer forecast surfacing.

2. **Live-loop durability — don't let the page go dark.** ([#223](https://github.com/antifungal-resistance-project/OpenAFR/issues/223))
   The display window anchors to the
   committed fill; on the current fill it goes dark ~2027-04-10. The non-billable monthly
   watchdog measures drift and rings for a human-gated billable refill. Bet: make sure that
   watchdog→refill path is solid and the dark-date is never a surprise.

3. **Raise the echinocandin-call coverage ceiling (FKS1 HS3 panel extension).** ([#224](https://github.com/antifungal-resistance-project/OpenAFR/issues/224))
   The #137
   diagnostic-accuracy arm read FAIL (VME 10% vs ≤3% bar) — a *coverage* ceiling, not a caller
   defect. The fix is the HS3 panel extension (v2 prereg; D642Y/M690I/R1354S path via #154,
   already on main). Billable GCP re-measure to confirm it lifts coverage.

4. **FKS1 structural so-what — optional deepening.** ([#225](https://github.com/antifungal-resistance-project/OpenAFR/issues/225))
   Rungs 1–3 shipped (position_context,
   structural_proximity, proximity_label verdict). Honesty boundary holds: no metal → no fit
   claim, never ACT-NOW. Any further structural work is a judgment call, not a committed need.

## Parked / blocked (do not start without the unblock)

- **Calibration engine (#136)** — blocked on a *C. albicans* continuous-score panel that does
  not exist. Novel regime has no continuous score by design. Leave parked.
- **External-truth backtests (ERG11 & FKS1 concordance beyond the 98-isolate PASS)** — the
  PMC12323592 truth set is not open-access and must be hand-transcribed from PDF before the
  grader can run. Optional; the deliverable prevalence numbers are already measured.

## Out of scope (by standing decision)

- **Outreach of any kind** (wet-lab, orgs, public-health agencies, collaborators) — software
  first. Revisit only if the user raises it.
- **Within-azole-class ranking** — closed for good (TODOS.md): pose geometry over a fixed
  CYP51 target cannot rank working azoles above failed analogues to the 0.70 bar.
- **Drug-discovery CYP51 triage** is the research/credibility arm, not the product.
