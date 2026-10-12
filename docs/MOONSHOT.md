# OpenAFR Resistance Weather — the moonshot

*Aspirational north star, not the active backlog. See [VISION.md](../VISION.md) for the
real product bet and the current state. This file answers one question: if we had all the
money and resources in the world, what would the biggest, most successful version of this
software do that it doesn't today?*

## The one-line moonshot

Today OpenAFR Resistance Weather is a **hindcast** of **one mechanism** (echinocandin /
FKS1) in **one bug** (*Candida auris*) from **one lagging feed** (public NCBI), shown on
**one page**. The biggest version is a **real-time, forecasting, global nowcast of all
antifungal resistance** that both *warns* public health and *feeds a drug-design loop that
beats the resistance it detects* — the authoritative "weather service" for antimicrobial
resistance.

## What it would do that it doesn't now

**1. Every bug × every drug — a resistance matrix, not a single cell.**
Expand from *C. auris*/FKS1 to the whole WHO fungal-priority list (*A. fumigatus*,
*C. glabrata*, *C. neoformans*, …) × every mechanism (azole/ERG11-CYP51, echinocandin/FKS,
polyene, flucytosine). Each cell is a separately benchmark-certified caller. Same
discipline, N targets.

**2. Near-zero-latency ingestion.**
Our very first fill already drifted out of the live window because public deposits lag
months. The max version pipes directly from sentinel labs, hospital sequencers, and the
CDC AR Lab Network — the picture is hours old, not quarters old.

**3. Actual *forecasting*, not just reporting.**
Real weather predicts. Layer phylogeographic spread + transmission modeling on the mutation
stream to project *where resistance appears next, which clade is moving, months out* —
"resistance front" maps and alerts, every forecast backtested.

**4. Structural verdicts on never-seen mutations.**
Today a caller matches *known* resistance alleles; surveillance goes blind on novel
variants. The pocket-verdict engine would structurally predict whether a brand-new mutation
confers resistance — so emergence is caught on first sighting, not after it's characterized.

**5. Cross the clinical line (deliberately deferred today).**
VISION.md parks per-patient use behind the regulatory bar. The maximal version clears it:
an FDA-cleared genotype→calibrated-treatment verdict wired into hospital systems — "this
isolate → echinocandin resistance predicted, switch therapy."

**6. Close the drug loop — and make the two arms a flywheel.**
Today the CYP51 triage tool is a one-shot shortlist with no wet lab. The big version is a
closed design→make→test loop (generative chemistry against resistant pockets + a robotic
assay partner feeding results back), expanded past CYP51 to new targets/mechanisms. **The
surveillance arm tells you which resistance is emerging; the discovery arm designs the
molecule to beat it; wet-lab results re-certify the callers.** That mutual feed is the real
moat.

**7. Become the reference layer.**
The "USGS earthquake feed" of resistance: a public API that WHO/CDC dashboards, stewardship
programs, and pharma pipeline teams all consume and cite.

**8. Generalize beyond fungi to all AMR.**
The pattern — certified genotype→resistance caller + live feed + backtested forecast —
transfers to bacterial AMR, a far larger killer. "Resistance Weather" for all of
antimicrobial resistance.

## The through-line

The edge stays the same at every scale: the **honesty discipline** — pre-registered,
benchmark-certified, calibrated-uncertainty, a gate that documents runs that *failed*
first. The moonshot is that same proof-before-claim rigor applied to *every* pathogen, in
*real time*, with a drug-design loop on the back end.
