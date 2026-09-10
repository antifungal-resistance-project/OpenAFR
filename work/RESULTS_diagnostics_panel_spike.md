# Findings — the calibration-grade panel hunt (issue #132), pass 1

> **RESOLVED 2026-09-09 → FAIL, fallback (a): concordance-only.** Per the frozen prereg's
> pre-committed FAIL path, no obtainable public *Candida* collection meets the calibration bar,
> and the failure is a *near-miss* (real paired data exists but is class-imbalanced /
> hotspot-skewed for C. auris, and too small for C. albicans). Decision (user, 2026-09-09):
> **narrow the day-one claim to concordance-only** — keep the hard S/R call + the honest
> "uncharacterized" verdict on the existing C. auris caller, and **drop "calibrated probability"
> from the #133 contract.** The calibrated-probability product (option C, a pooled C. albicans
> panel) is deferred to a separate track, not killed. Panel size that *would* unlock it is
> recorded below (N≥150, ≥40/class, ≥5 variants ≥5 carriers). The frozen bar was NOT loosened.

> Run against the frozen bar in `work/PREREGISTRATION_diagnostics_panel.md`
> (sha `1211e973…`). This is an honest interim finding, not a certification. **No panel is
> certified PASS in this pass, and none is declared a final FAIL.** The verdict is *not yet
> resolvable without obtaining specific paywalled supplementary tables* — which, per the
> prereg's own criterion #5, is a manual hashed-transcription step, not an LLM scrape.

## The landscape, mapped against the 5 criteria

| Candidate | N (WGS) | Class balance | ERG11 variant diversity | Genotype+MIC paired | Public reads | Read vs bar |
|---|---|---|---|---|---|---|
| **C. auris — Chowdhary 2018 JAC** (350 isolates, India, 2009–17; ERG11+FKS1+AST) | ~350 | **Fails** — clinical C. auris is overwhelmingly FLC-resistant; susceptible ≪ 40 | **Likely fails** — Clade I, dominated by Y132F/K143R hotspots | Reported in paper tables | Some BioProjects exist | **Near-miss (crit #2, #3)** |
| **C. auris — JCM 2025 transmission/AR set** (~188; "181/188 with FLC MIC >32 had ERG11 mut") | ~188 | **Fails** — almost all resistant | **Likely fails** — VF125AL (clade III) / K143R clade-linked | WGS + FLC MIC | Likely SRA | **Near-miss (crit #2, #3)** |
| **Multi-species — Canadian reference lab, AAC 2020 (McTaggart)** | **4,715 AST'd, but WGS on only ~11** in-host-evolution cases | n/a | n/a | **No** — it is a *susceptibility survey*; WGS applied to a handful, not a paired panel | partial | **FAIL (crit #2 N) — collapses on inspection; not a calibration panel** |
| **C. albicans — EUCAST-above-breakpoint set, PRJNA592373** (40 isolates; ERG11/ERG3/TAC1/FKS1; 54 missense) | 40 | **Fails** — all above breakpoint (resistant-only by design) | **Diverse** (54 missense) | Yes + **public SRA** | **PRJNA592373** | **Near-miss (crit #2 N & class)** — diverse + public but small & resistant-only |
| **C. albicans — pooled ERG11 literature** (Flowers 2015; Morio; Iranian/Chinese cohorts) | per-study 4–100 | Balanced within studies | **Diverse** — ~20 substitutions (Y132F, K143R, Y123H, G448E, G464S, S405F, S154F…) | Per study | Fragmented | **Fails as one set (crit #2 N, #5 single obtainable collection)** |

## The structural tension (the real result of this pass)

1. **The default target (C. auris ERG11) is trending toward FAIL on the frozen bar — but as a
   *near-miss*, not "no data."** Big-N collections exist (350, 188), so the data is real; they
   fail on **class balance** (clinical C. auris is almost all azole-resistant — too few
   susceptibles to anchor the "S" end of a calibration curve) and on **variant diversity**
   (ERG11 resistance is a few clade-linked hotspots). This is exactly the prereg's fallback-(a)
   trigger: *narrow the day-one claim to concordance-only*, not force-fit a skewed panel.

2. **The organism that best clears the bar is not the organism the callers were built for.**
   *C. albicans* ERG11 has the diversity (~20 documented substitutions) and the susceptible
   isolates that calibration needs — but our re-callers and the structural port are
   *auris*-centric ([[openafr-auris-port]], 5TZ1). So "PASS the panel bar" and "reuse the
   existing caller" pull toward different organisms. That tension must be decided, not glossed.

3. **Certification is paywall-gated, same wall as FKS1.** The single best PASS candidate (the
   Canadian reference-lab WGS+MIC collection) is behind a 403; verifying per-isolate
   genotype+MIC, counting per-variant carriers, and confirming SRA accessions needs the actual
   supplement. Per prereg crit #5 this is a **verified manual transcription**, mirroring the
   FKS1 blocker ([[fks1-caller-concordance]]) — not something to scrape or infer.

## Verdict of this pass — leaning FAIL-on-a-single-collection; decision required

After pass 1 no existing public dataset clears the frozen bar as a **single obtainable
collection**, and the pattern is structural, not bad luck:

- **Big-N paired sets are C. auris** (Chowdhary 350, JCM 188) → **fail class balance + diversity**
  (mostly resistant, few clade-linked ERG11 hotspots).
- **Diverse, class-relevant sets are C. albicans but small** (PRJNA592373 n=40 resistant-only;
  APECED n=14; pooled literature) → **fail N≥150 and/or class balance**.
- **The one "large" collection (AAC 2020) is a susceptibility survey**, not a paired WGS+MIC
  panel — WGS on ~11, not 150+ → **fail N**.

So the frozen bar, as written (single collection, N≥150, ≥40/class, ≥5 variants ≥5 carriers,
paired genotype+MIC), is **not met by any public Candida dataset found**. The only routes to a
genuine PASS are: **pool** several public C. albicans studies into one hashed fixture (needs a
prereg *amendment* — pooling mixes CLSI/EUCAST MIC methods across labs, a real heterogeneity
problem crit #5's "single collection" rule was guarding against), or obtain a not-yet-found
larger paired set behind access we don't have.

## Next actions (a target-strategy fork — pick before any build)

Candidate (A) from pass 1 (the AAC 2020 collection) is **withdrawn** — it is a survey, not a
panel. The live options:

- **(B) Accept the near-miss → concordance-only on C. auris.** Keep the existing auris caller,
  drop "calibrated probability" from the #133 contract, ship the honest S/R + "uncharacterized"
  verdict validated against the C. auris benchmarks as a *concordance* check. Cheapest; smallest
  claim; buildable now on the callers we have.
- **(C) Pool a C. albicans calibration panel (prereg amendment).** Assemble PRJNA592373 + other
  public C. albicans ERG11 studies into one hashed fixture to reach N≥150 with diversity. Best
  shot at a real *calibrated-probability* product — but needs a frozen amendment addressing
  cross-study MIC-method heterogeneity, and new C. albicans caller/structure work (our port is
  auris/5TZ1).
- **(D) Keep hunting** for a single larger paired set (incl. Aspergillus CYP51A per the flip
  rule) before deciding — costs another research pass, no build.

Recommendation: **(B) as the shippable MVP now, with (C) as the flagship track in parallel.**
The data says the default target can carry an *honest concordance* claim today but not the
*calibrated-probability* claim; splitting them keeps #132's kill gate honest instead of quietly
loosening the bar to force a PASS.
