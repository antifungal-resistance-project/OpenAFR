# OpenAFR — Vision & Strategy

Strategic framing for the project. The [README](README.md) covers *what the pipeline
does and how to run it*; this covers *what it's for, what's genuinely new, and where it
goes next.* Written to be picked up cold.

## What the tool is

A **triage tool for antifungal drug discovery**, aimed at CYP51 — the fungal enzyme the
azole drugs jam and that drug-resistant *Candida auris* is defeating.

Its one distinctive move: **it does not trust the docking program's affinity score.**
Instead it ranks molecules by a *geometric, mechanistic* criterion — how closely a
nitrogen atom reaches the heme iron at the enzyme's core, which is the actual mechanism
by which azoles kill. On a blind held-out test that criterion hit AUC 0.794 (p=0.0028),
while AutoDock Vina's own score ranked the real drugs *worse than random* (AUC 0.471) on
the identical poses. **Geometry beats the scoring function on this system.**

**One caveat, established later and carried honestly.** The same single-search run also
reported EF@1% 12.68x, but a pre-registered reliability check
([RESULTS_reliability.md](work/RESULTS_reliability.md)) found that *top-1%* figure was an
artifact of docking once — it collapses to 0.00x under 5-seed sampling, under every
aggregation tried. What is robust is *broad* enrichment (AUC ~0.79–0.85, EF@5% ~5.6–8.5x).
So the defensible headline is a **trustworthy top-~5% shortlist, not razor-top ranking**;
the AUC is the durable number, EF@1% is not.

**Two later rounds of pre-registered stress-testing sharpened the product.** Re-run against
*measured* (not presumed) inactives the criterion passes — and locates a ceiling: it separates
**novel chemotypes** well (AUC ≈ 0.81) but cannot rank a *working* azole above a *failed* azole
analogue (AUC ≈ 0.65), a limit confirmed on both the rigid and the crystallographic-ensemble
receptor. So the honest product is **novel-chemotype triage**, and the frozen criterion also
generalizes to never-scored and newest-literature holdouts (temporal AUC 0.754). On top of the
raw rank the tool now has a **front door** (score your own SMILES) and a **candidate-confidence
layer** — precedent (already measured and failed?), ADMET, human-CYP off-target liability,
synthesizability, and a *coordinator-identity* check that catches out-of-mechanism iron-reach —
fused into a per-candidate A/B/C dossier. Every axis annotates a rank; none re-ranks or invents.

Everything else — pre-registration, hash-frozen protocol, a validation gate that blocks
the screen if it can't re-discover known drugs it has never seen — exists to *prove the
ranking is trustworthy before it is allowed to rank anything new.* That discipline is the
product today, more than any single molecule.

## What's actually novel (be honest)

- **The ingredients are not new.** Mechanism/geometry-based rescoring of docking poses
  (pharmacophore constraints, metal-coordination constraints) is a known, decades-old
  idea. We are not the first to distrust a docking score.
- **The discipline and the specific validated result are the contribution.** A
  pre-registered, hash-frozen, blind-holdout, permutation-tested screen — with a gate
  that documents a run that *failed* first — is rare, especially in open/amateur docking
  repos that publish ranked lists with zero evidence the ranking means anything.

Frame it as *"a known idea, actually proven, honestly, on a target that matters"* — not
as a novel algorithm. The honesty is the edge; overselling it forfeits that.

## End goal & how it's used

- **Not a run-forever cron job.** It's hypothesis generation / triage: point it at a
  candidate library → dock once under the frozen protocol → get a short ranked shortlist
  → hand the top of that list to a wet lab. Re-run when there's a *new library*, not on a
  schedule.
- **A high rank is a hypothesis, never a hit.** The tool makes a wet-lab yes/no *more
  likely to succeed than random* and cheaper to reach. It does not make a drug.

## New drug vs. better drug vs. new mechanism

- Aimed at **CYP51 (the azole target)**, the default output is *a new molecule against
  the same target* — ideally one whose geometry still coordinates the iron **in resistant
  strains** where mutations reshaped the pocket. Valuable (resistance is the problem), but
  it's the same **mechanism of action** azoles use.
- A genuinely **new mechanism** requires re-aiming the pipeline at a *different* essential
  fungal protein. The method transfers; the target is a deliberate choice. **Strategic
  fork:** "beat resistance on the proven target" vs. "open a new target." The tool is
  currently set up for the former.

## The wet-lab hand-off (the critical open dependency)

No wet-lab collaborator yet — **this is the one blocking dependency.** Once a shortlist
meets a partner, the funnel is roughly:

1. Acquire / synthesize the top candidates.
2. **MIC assays** on live *C. auris*, including resistant strains — the first real yes/no.
   Most computational hits die here; that's expected.
3. **Selectivity / cytotoxicity** — kill the fungus without killing human cells (both are
   eukaryotes; this is the hard part of antifungals specifically).
4. Survivors → hit-to-lead → lead optimization → ADME/tox → animal models.
5. IND → Phase I/II/III. Years and tens-to-hundreds of millions — so this stage almost
   always happens via **licensing / partnership**, not the discoverer going it alone.

## Portability — is the tool only good for antifungals?

No. Two separable, reusable assets:

- **Mechanism-anchored geometric rescoring** — works for any target with a knowable
  structural mechanism (metalloenzymes, clear pharmacophores). Antibacterials, other
  enzyme inhibitors, etc.
- **The self-validating gate** — *fully domain-independent.* "Don't trust a screening
  ranking until it re-discovers what you already know." Applies to essentially any virtual
  screen. Arguably the most portable and most novel-*feeling* contribution.

The antifungal target is the **demonstration case, not the ceiling.**

## Second track: genomic early-warning surveillance

A second track now exists in the repo (issues #20–27): watch **NCBI Pathogen Detection**
for emerging *C. auris* azole-resistance mutations and give each flagged mutation a
**structural so-what** via the same CYP51 pocket tooling. Track 1 finds new drugs; track 2
watches resistance evolve against the drugs we have. The structural verdict is what keeps it
from being a bare novelty feed — it reuses track 1's moat rather than starting a new one.

**Honest status — scaffold complete, re-caller run, track redirected.** The full chain
(snapshot → emergence → mapping → structural so-what → alert → scheduled delivery) is built,
tested, and pulls real NCBI data. The killer finding still holds: **NCBI runs no AMR pipeline on
*C. auris*** — it's a metadata + genome-pointer feed, not a resistance feed — so the resistance
signal must be manufactured by *us* from SRA reads, via an **ERG11 re-caller**. That re-caller
has now been **built and run** on a representative real sample, measuring the number the track
was gated on: the **azole event frequency = 80.4%** (Wilson 95% CI 74.3–85.3%, n=199;
[work/RESULTS_prevalence.md](work/RESULTS_prevalence.md)). That result *redirected* the track:
azole resistance is near-saturated at baseline, so azole *emergence* is a weak early-warning
signal — and a parallel **FKS1/echinocandin detection track** was built (honestly detection-only;
echinocandins coordinate no metal, so the CYP51 structural moat does not transfer). Fully
*validating* a warning still needs a wet-lab-anchored backtest truth set, so nothing over-claims
— but the scaffold is real, the deliverable is measured, and the gaps are named precisely.

## Foundation scope — diagnostics and the rest of resistance?

**Recommendation: broad mission, focused execution.**

- Antifungal resistance is much larger than drug discovery, and **diagnostics is arguably
  the most neglected piece** — fungal disease is diagnosed slowly and poorly, and people
  die of *C. auris* partly because ID and susceptibility testing lag. A foundation named
  "The Antifungal Resistance Project" legitimately spans drugs, diagnostics, surveillance,
  and stewardship, and framing the *mission* that broadly is honest and compelling.
- **But** the asset today is a computational drug-discovery tool. Diagnostics is a
  different competency (wet lab, clinical validation, regulatory, sometimes hardware).
  Spreading a small effort across both before either has a validated result or a partner
  is how good projects stall.
- **So:** let the *mission* be broad (all of antifungal resistance); keep the *active
  project* narrow until it reaches (a) a validated wet-lab hit or (b) a wet-lab
  collaborator. Diagnostics becomes the natural *second* project afterward.
- **Bridge:** discovery and diagnostics share a natural partner — a clinical/mycology wet
  lab. The relationship that tests your candidates could be the one that later opens the
  diagnostics work. Sequential, not competing.
- **Where the early-warning track fits:** surveillance (the second track above) is the one
  adjacency taken on before a validated hit, defensible *because it reuses track 1's structural
  moat* rather than opening a new competency — not the diagnostics expansion this section
  cautions against. It has taken on exactly **one** data-justified widening: **FKS1/echinocandin
  detection**, opened only *after* the ERG11 re-caller measured an ~80% azole baseline that makes
  azole emergence a weak signal. That FKS1 half is kept honest as **detection-only** — no
  structural verdict, since echinocandins coordinate no metal and the moat does not transfer.
  Everything past that (a structural FKS1 verdict, TAC1B/ERG3, diagnostics) stays out of scope
  until the current work is validated (see [TODOS.md](TODOS.md)).
