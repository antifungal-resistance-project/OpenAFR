# Customer-discovery brief — track B (issue #165)

> **STATUS: PRE-INTERVIEW BRIEF — HYPOTHESES, NOT YET VALIDATED (2026-10-11).**
> **No labs have been interviewed yet.** Everything below is a plan and a set of
> falsifiable guesses, not findings. There are no real quotes, no named commitments, and
> no demand signal in this document. Its job is to make the eventual 3–5 conversations
> fast, non-leading, and decision-grade — and to record the regulatory boundary the
> product already commits to. When interviews happen, their notes go in a **separate**
> "Findings" file so this brief and the evidence never blur.

Child of the surveillance-product epic (#170). This is the **track-B gate**
([[direction-and-workflow]]): per the charter, *no B-side build or commercialization work
starts until this gate clears.* The question it must answer is the plain one from #165 —
**who would pay for (or fund) a research-grade resistance-interpretation tool, for what, and
what is the regulatory reality** (non-regulated surveillance software vs. a clinical
diagnostic).

See also: [`../VISION.md`](../VISION.md) (the product bet), [`BUSINESS_PLAN.md`](BUSINESS_PLAN.md)
(the foundation/mission structure), [`DIAGNOSTIC_GO_NO_GO.md`](DIAGNOSTIC_GO_NO_GO.md) and
[`DIAGNOSTIC_VERDICT_CONTRACT.md`](DIAGNOSTIC_VERDICT_CONTRACT.md) (the RUO boundary this brief
must not contradict).

---

## 1. What we are actually selling (so the interviews stay honest)

Two distinct things, at two very different regulatory bars. Keep them separate in every
conversation — conflating them is the fastest way to a false positive.

1. **OpenAFR Resistance Weather** — a live, public, auto-updating report on *emerging
   echinocandin resistance in C. auris*, backed by a benchmark-certified FKS1 caller. This is
   **population-level surveillance**, not a per-patient result. It is the beachhead
   (`VISION.md`, product bet 2026-10-06).
2. **The resistance-interpretation engine** — takes a *typed C. auris genotype* and emits a
   per-drug-class categorical verdict (`RESISTANCE_MARKER_DETECTED` / `UNCHARACTERIZED_VARIANT`
   / `NO_KNOWN_MARKER` / `UNRESOLVED`). It is **RUO** — research use only; it emits **no
   probability, no MIC, no susceptibility call, and no clinical recommendation**
   (`DIAGNOSTIC_VERDICT_CONTRACT.md`, `DIAGNOSTIC_GO_NO_GO.md` §3). The honest boundary *is*
   the product.

The hardest truth to hold in discovery: **the beachhead may be non-paying public-health
infrastructure.** "Who uses it" (a PH-lab epidemiologist) and "who funds it" (a CDC contract,
a grant, a philanthropy) can be different people. The interviews must surface *both* — usage
pain **and** who controls the budget that would keep a surveillance feed running.

---

## 2. Segment map

For each segment: who, the job they'd hire this for, why they might pay or fund (vs. use it
free), and the regulatory bar they operate under. **All demand below is hypothesized.**

### A. Public-health / genomic-surveillance labs — **beachhead**
- **Who:** CDC Antimicrobial Resistance (AR) Lab Network and its regional labs; state public-
  health labs; dedicated *C. auris* genomics groups; APHL-member labs.
- **Job:** "Tell me whether echinocandin resistance is emerging in the *C. auris* genomes
  flowing through NCBI/our pipeline, without my team hand-curating FKS1 calls."
- **Pay / fund:** Unlikely to "buy" a tool; likelier path is **free public infrastructure**
  funded by a CDC cooperative agreement, an AMR surveillance grant, or philanthropy, or
  adopted as a method/validation reference. WTP ≈ low; **willingness-to-fund / willingness-to-
  adopt** is the real axis.
- **Regulatory bar:** **Low** — population surveillance is not a diagnostic; see §5. This is
  why it's the beachhead.

### B. Clinical microbiology labs — **later, higher bar**
- **Who:** hospital/reference clinical micro labs making per-patient susceptibility calls.
- **Job:** "Given this patient's isolate genotype, is it echinocandin-resistant?" — a
  per-patient actionable result.
- **Pay / fund:** Real budgets exist for diagnostics, but this is where the product would
  become a **clinical test** and cross into CLIA/validation territory (§5). The current engine
  **explicitly refuses this claim** (RUO; no S call, no probability — `DIAGNOSTIC_GO_NO_GO.md`
  §3). Interview them to *learn the pull and the bar*, **not** to sell today.
- **Regulatory bar:** **High** — see §5.

### C. Pharma / biotech surveillance & epidemiology teams
- **Who:** antifungal developers, clinical-trial epidemiology groups, AMR market-access teams.
- **Job:** "Track the resistance landscape for a target/compound; size an emerging-resistance
  population."
- **Pay / fund:** **Plausible WTP** — commercial budgets, data/report subscriptions. This is
  the clearest "who pays" candidate even though it's not the mission beachhead.
- **Regulatory bar:** Low (RUO/informational).

### D. Academic *C. auris* / AMR researchers
- **Who:** mycology and AMR genomics labs, method developers.
- **Job:** "A certified FKS1 caller and an honest genotype→verdict reference I can cite and
  build on."
- **Pay / fund:** ≈ free/open; value is **credibility, citations, and collaborators**, not
  revenue. Feeds the research/credibility arm (`VISION.md`).
- **Regulatory bar:** None (research).

---

## 3. Core hypotheses to test

Numbered so interview notes can map evidence back to each. Each has a **confirm** and a
**refute** signal. A hypothesis is only "supported" by *past behavior or committed budget*,
never by "that sounds useful."

- **H1 — Pain exists.** Surveillance/PH labs currently find emerging echinocandin resistance
  in *C. auris* hard, slow, or manual.
  - *Confirm:* they describe a real recent instance of hand-curating FKS1 calls or missing an
    emergence signal. *Refute:* "an existing pipeline/tool already handles this fine."
- **H2 — Our seam is the bottleneck.** The specific hard part is *trustworthy FKS1 calling at
  scale over the live feed*, not data access or phenotype.
  - *Confirm:* they name caller trust / curation burden as the limiter. *Refute:* the real
    blocker is sequencing throughput, metadata, or MIC data — things we don't provide.
- **H3 — A certified-caller surveillance feed is trusted/adoptable.** A benchmark-certified
  caller + an honest freshness/coverage label is enough to use or cite.
  - *Confirm:* they state what certification/validation evidence they'd need, and it's within
    reach (we already have a 98-isolate concordance cert — [[fks1-caller-concordance]]).
    *Refute:* they'd only trust an in-house or CDC-run pipeline regardless of evidence.
- **H4 — The RUO / surveillance boundary is acceptable.** "No per-patient call, no
  probability, no MIC" is fine for their use.
  - *Confirm:* surveillance use needs no clinical actionability. *Refute:* the only valuable
    version is the per-patient clinical call — i.e., they want segment B, the high-bar product.
- **H5 — Someone will fund or pay.** There is a budget line (grant, CDC agreement,
  philanthropy, or commercial subscription) that would sustain a live feed.
  - *Confirm:* they can name the budget and who controls it. *Refute:* "great, but there's no
    money for it / it would have to be free forever with no funder."
- **H6 — Freshness cadence matches need.** Monthly (or similar) updates meet the surveillance
  need given NCBI's ~97-day median arrival lag ([[weather-fks1-render-wired]]).
  - *Confirm:* monthly is "plenty." *Refute:* they need near-real-time, which our feed can't
    honestly promise.

---

## 4. Interview guide

Mom-Test discipline: ask about **past behavior and real instances**, never "would you use /
pay for." Open questions, no pitching until the end. 15–20 minutes.

**Warm-up / context (all segments)**
1. Walk me through what your lab does with *C. auris* genomic data today.
2. The last time you needed to know about echinocandin resistance in a set of isolates — what
   did you actually do, step by step?
3. Who on the team does the resistance calling, and what tools do they use?

**Pain (H1, H2)**
4. Where does that process break down or eat the most time?
5. Tell me about a time a resistance signal was caught late, or missed. What happened?
6. When you get an FKS1 genotype, how do you decide what it means? Who do you trust for that?

**Trust & boundary (H3, H4)**
7. If a tool handed you FKS1 resistance calls over the live NCBI feed, what would you need to
   see before you trusted or cited it?
8. How do you feel about a tool that says "I detect these named markers and explicitly abstain
   on everything else," rather than giving a yes/no for every isolate?
9. Is population-level surveillance useful to you on its own, or do you only care about
   per-patient answers? (branch: if per-patient → segment B, probe the clinical bar.)

**Cadence & funding (H5, H6)**
10. How fresh would such a feed need to be to matter — daily, monthly, quarterly?
11. If this existed and worked, who would own it in your world — and whose budget keeps a
    service like it running? Have you funded anything similar?
12. What have you (or peers) actually paid for in this space — data, reports, software,
    validation services?

**Clinical-lab branch (segment B only)**
13. What would it take — validation, regulatory, operational — to use a genotype-based
    resistance result for a patient?
14. Who signs off, and against what standard (CLIA validation, in-house LDT process)?

**Close (all)**
15. Who else should I talk to about this?
16. If I built a first version, would you want to see it? (interest signal, not a sale)

---

## 5. Regulatory reality (verified 2026-10-11)

The decision-critical section. It separates the **low-bar surveillance** product from the
**high-bar clinical** one, and it is deliberately consistent with the RUO posture already
fixed in `DIAGNOSTIC_GO_NO_GO.md` §3.

**Two regimes:**

1. **Non-regulated public-health surveillance software (the beachhead).** A tool that reports
   *population-level* resistance trends from public genomic data, makes **no per-patient
   claim**, and carries **no clinical actionability**, is **not** an in-vitro diagnostic (IVD)
   and is not a medical device. This is the regime OpenAFR Resistance Weather operates in by
   design. It is also exactly why the engine's contract forbids probability, MIC, S-calls, and
   treatment guidance (`DIAGNOSTIC_VERDICT_CONTRACT.md`): staying RUO/surveillance keeps the
   product out of device regulation.

2. **Clinical diagnostic (IVD / laboratory-developed test).** The moment a result is used to
   guide the care of an **identified patient**, it becomes a clinical test. In the U.S. that
   means **CLIA** oversight of the performing laboratory, administered by **CMS** — the
   lab must validate the test and hold appropriate CLIA certification. A test designed,
   validated, and run within a single CLIA lab is a **laboratory-developed test (LDT)**.

**FDA LDT status — time-sensitive, confirmed current as of 2026-10-11:**
- FDA issued a final rule in **May 2024** asserting authority to regulate LDTs as medical
  devices (phased enforcement).
- On **March 31, 2025**, the U.S. District Court for the Eastern District of Texas **vacated**
  that rule, holding FDA lacks statutory authority over LDTs — Congress assigned laboratory-
  test oversight to **CMS under CLIA (1988)**, not FDA.
- HHS **did not appeal**, and on **September 19, 2025**, FDA formally **reverted** the
  regulatory text to its pre-May-2024 state.
- **Net, 2026:** LDT oversight sits with **CLIA/CMS**, not FDA device regulation. FDA
  continues to regulate commercially distributed IVD *kits* and retains **RUO labeling** rules
  (an RUO product must not be marketed for clinical/diagnostic use).

**What this means for OpenAFR (and for the interviews):**
- Staying in regime 1 (surveillance, RUO) keeps the regulatory bar **low** — the right place
  to launch and the reason segment A is the beachhead.
- Any move toward regime 2 (segment B, per-patient calls) brings **CLIA validation** into
  scope for whichever lab runs it, plus the calibrated-accuracy evidence the engine does not
  yet have (`DIAGNOSTIC_GO_NO_GO.md` §4 gate ladder; calibration track #136 is blocked,
  [[calibration-track-blocked]]). Discovery with clinical labs is to **map that bar**, not to
  promise crossing it.
- RUO labeling is a genuine constraint, not a formality: the product must not be marketed as
  clinically actionable. The contract already enforces this.

---

## 6. Go / no-go for any track-B build work

Mirrors the gate discipline of `DIAGNOSTIC_GO_NO_GO.md`: the gate is **NO-GO until the
evidence below is in hand.** Meeting it is necessary, not sufficient.

**GATE (what we must hear from ≥3 of the 3–5 interviews before committing B-build effort):**
- **H1 + H2 supported:** at least 3 contacts describe a *real, recent* instance where
  trustworthy FKS1 calling over the feed was the bottleneck (past behavior, not speculation).
- **H4 supported:** the RUO / surveillance boundary is acceptable to the beachhead segment —
  i.e., surveillance value stands on its own without per-patient claims.
- **H5 directional:** at least one credibly-named funding or payment path (a budget line, a
  grant mechanism, a CDC agreement, or a commercial subscriber) — not a vague "someone should
  fund this."

**Current status: NO-GO (gate not yet attempted — zero interviews conducted).** This brief
makes the gate runnable; it does not clear it.

**If the gate fails** (pain is thin, or the only value is the high-bar clinical call with no
funder): track B stays deferred and effort returns to the research/credibility arm and the
live-loop operations (#171), per `VISION.md`.

---

## 7. How to run this (operational)

1. Build the target list: name 3–5 specific contacts across segments A (≥3, the beachhead) and
   B (1–2, to map the clinical bar). Channels: CDC AR Lab Network / Mycotic Diseases Branch,
   APHL membership, authors on recent *C. auris* FKS1 surveillance papers, state PH-lab
   directors.
2. Outreach template (adapt per contact):
   > *Subject: 15 min on C. auris echinocandin-resistance surveillance?*
   > I'm building an open, certified FKS1 caller that reports emerging echinocandin resistance
   > in C. auris from public genomic data (research use, population-level — no patient calls).
   > Before I build further I'm trying to learn how labs like yours actually handle this today.
   > Could I ask you 15 minutes about your current workflow? Not selling anything — just
   > learning where it breaks.
3. Run interviews with §4 guide; log notes in a **new** `docs/CUSTOMER_DISCOVERY_FINDINGS.md`
   (one entry per conversation, dated, attributed only with consent), mapping each note to the
   H# it bears on.
4. Re-evaluate §6 gate once ≥3 interviews are done.
