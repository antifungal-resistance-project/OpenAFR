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
by which azoles kill. On a blind held-out test that criterion hit AUC 0.794 / EF 12.68x
(p=0.0028), while AutoDock Vina's own score ranked the real drugs *worse than random*
(AUC 0.471) on the identical poses. **Geometry beats the scoring function on this system.**

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
