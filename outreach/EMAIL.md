# Cold outreach — template, targeting, and logistics

Companion to [ONE_PAGER.md](ONE_PAGER.md). The one-pager is the attachment; this is how it gets sent
and to whom.

---

## The email

Keep it short. The one-pager does the work; the email's only job is to get it opened. Do **not**
paste the shortlist into the email body — it reads as spam and it front-loads the claim instead of
the evidence.

> **Subject:** 10 approved drugs, pre-registered CYP51 shortlist — would you run MICs?
>
> Dear Dr. ——,
>
> I run a small open-source computational project on antifungal resistance, and I have a shortlist I
> can't test. I'm writing to ask whether your lab would run MIC panels on ten approved or
> clinical-stage drugs against *C. albicans* and *C. auris* — about a week of bench work.
>
> The shortlist comes from a virtual screen of the Broad Repurposing Hub against *Candida* CYP51,
> ranked by heme-iron coordination geometry rather than by the docking score. Two things make it
> worth a look rather than being another ranked list:
>
> - On a blind, pre-registered held-out set the docking score ranked known azoles **worse than
>   random** (AUC 0.471); the geometric criterion on the same poses gave AUC 0.794, p = 0.0028. And
>   screening 2,776 unlabeled drugs blind, it put seven known azoles in the top 3.2%.
> - Four pre-registered follow-up tests **failed**, and I've published all four — including the one
>   showing the top-1% enrichment is an artifact of docking once. What survives replication is broad
>   enrichment (~8x over random in the top 5%), not razor-top ranking. I'd rather send you the
>   honest version of the claim than the impressive one.
>
> Attached is a one-page summary with the compounds, the suggested panel, and — at more length than
> is probably normal — everything I already know is wrong with it, including that the method says
> nothing about fungal-vs-human selectivity.
>
> I have no commercial interest here; the repository is noncommercial-licensed and there's no
> company behind it. The result gets published either way, a null result included, and I'd write the
> pre-registration with you before anything is run. The experimental work and the authorship
> position that goes with it are yours.
>
> If this isn't a fit for your lab, I'd be grateful for a pointer to someone it might suit — or for
> a sentence telling me the shortlist is uninteresting and why, which would genuinely be useful.
>
> Best regards,
> Jacob Jensen
> The Antifungal Resistance Project · https://github.com/antifungal-resistance-project/OpenAFR
> `<CONTACT EMAIL>`

**Attachment:** `ONE_PAGER.md` exported to PDF (`/make-pdf outreach/ONE_PAGER.md`). Attach the PDF,
link the repo. Don't attach the preprint unless asked — offer it if they reply.

### Why it's written this way

- **The negatives are in the second bullet, not buried.** Leading with four published failures is
  the single most credible thing available, and it's what separates this from the mass of unsolicited
  docking lists that mycology labs already ignore.
- **The ask is bounded and concrete** — a week, a standard panel, ten catalog compounds. Vague
  "would you like to collaborate" emails don't get answered.
- **Authorship and funding are addressed before they have to ask.** The most common silent objection
  to an unsolicited computational collaboration is "what does this person want from me."
- **The exit ramp at the end** ("tell me it's uninteresting and why") converts a large fraction of
  otherwise-silent nos into useful information. Take those replies seriously — three people
  independently saying the shortlist is chemically implausible is a result.

## Who to send it to

Send **~25–30**, individually, over about two weeks. Expect a low single-digit number of replies and
one yes if it goes well. Personalize one sentence per email — reference an actual paper of theirs;
generic mail-merge is visible and fatal.

Priority tiers:

1. **Academic mycology labs already doing *C. auris* susceptibility work.** The best filter is
   recent papers reporting MICs on *C. auris* clinical isolates — those labs have the strains, the
   protocol, and the CLSI plates already running. Search PubMed for *C. auris* + MIC + ERG11 from the
   last three years and write to corresponding authors.
2. **Labs that have published on ERG11/Y132F specifically.** They will immediately understand the
   transfer-test negative and are the most likely to find it interesting rather than disqualifying.
3. **Antifungal drug-repurposing groups.** They already believe the premise; the sell is the
   validation discipline, not the idea.
4. **Public-health and reference labs** — CDC Mycotic Diseases Branch, UK HSA Mycology Reference
   Laboratory, and the equivalents. Slower, more bureaucratic, but they have the isolate collections
   and a mandate that matches the mission.
5. **Contract labs** as the paid fallback if no academic partner materializes in ~2 months. A
   ten-compound MIC panel is a routine quoted service. This costs money but removes the dependency
   entirely, and it's a small enough number that it may be worth just paying for.

Also worth a separate, differently-framed email: **MMV / Medicines for Malaria Venture-style open
compound-screening consortia** and **the Gates Foundation-adjacent AMR programs** — they run open
screening cascades and sometimes accept external nominations.

## Logistics to have settled before the first email goes out

So that a yes doesn't stall:

- **Strains.** CDC & FDA Antibiotic Resistance Isolate Bank distributes a *C. auris* panel (multiple
  clades, including fluconazole-resistant) free to qualifying labs. Confirm the current request
  process and be able to say so in the reply.
- **Compounds.** Price the ten at MedChemExpress/Selleck/Sigma and have a real number, not the
  "$100–300 each" estimate. Check availability — some clinical-stage compounds (taranabant,
  R-1479/balapiravir, LDN-27219, nolatrexed) may be harder to source than the approved ones, and if
  two are unobtainable it's better to know before proposing them. Consider whether to offer to buy
  the compounds yourself; it is the cheapest possible way to remove friction from a yes.
- **The pre-registration.** Draft the MIC pre-registration in advance — endpoints, breakpoints,
  what counts as a hit, what gets published on a null — so it can be sent within a day of interest.
  A collaborator who sees the pre-registration discipline applied to *their* experiment, before they
  commit, is far more likely to say yes. Same hash-frozen format as `work/PREREGISTRATION_*.md`.
- **What you'll do with a partial yes.** A lab may offer *C. albicans* only, or three compounds
  instead of ten. Take it. In-domain *C. albicans* data is the more defensible test anyway, given
  the Y132F transfer failure.

## Sequencing note

Post the preprint (`work/PREPRINT_geometry_ceiling.md`) publicly **before or alongside** the
outreach, not after. A citable, timestamped writeup that includes four published failures is what
makes a cold email from an unknown project legible to an academic — it moves you from "someone with
a spreadsheet" to "someone with a preprint and a stopping rule." bioRxiv or Zenodo; either gives a
DOI to put in the signature line.

## Tracking

Keep a simple table in this directory — recipient, institution, date sent, personalization used,
reply, outcome. Not for optimization, but so that at 30 sent with 0 yeses you can look at the
pattern and tell whether the problem is the list, the framing, or the premise. Set that review point
now, in advance: **if 30 emails produce no wet-lab interest and no useful critique, that is
information about the project, not about the emails.**
