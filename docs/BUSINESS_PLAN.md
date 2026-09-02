# The Antifungal Resistance Project — Business Plan

How the software tool turns into a drug company that actually gets a new
antifungal to patients, while the foundation stays in control and the mission
stays first.

Read alongside [VISION.md](../VISION.md) (what the tool is) and
[CHARTER.md](CHARTER.md) (the rules that keep us honest).

> **Where we really are.** Today all we have is a computational triage tool with
> one honest result: on CYP51, geometry beats the docking score (AUC ~0.79). We
> have no wet-lab partner, no validated hit, and no molecule. This plan is the
> path from here to a real drug. Almost every step has a point where most projects
> die. I've written it so that failing a step is survivable and we say so out
> loud. The dollar amounts and timelines are rough planning ranges, not promises.

---

## 1. The short version

The idea in one breath: the software finds a promising molecule, a small company
develops it with the help of outside labs, and the foundation owns enough of that
company that nobody can ever sell the mission out from under us. Money that comes
back gets spent on the next drug, on diagnostics, on surveillance — not on cashing
anyone out.

Three things I care about, and everything in the plan follows from them:

1. **The mission stays in control.** We copy the Novo Nordisk Foundation setup: a
   foundation that owns the drug company through special shares, so the company
   answers to the mission instead of to whoever wants to buy it.
2. **The team stays tiny.** We own the science and the IP and rent everything else
   from contract research organizations (CROs). A big payroll is what forces bad
   decisions, so we keep it small on purpose.
3. **We use money that doesn't cost us ownership.** Grants, prizes, and
   mission-friendly funders first. We only take investment money when there's no
   other way, and even then we structure it so the foundation keeps the wheel.

The loop I'm building: the tool finds something, the foundation owns it, outside
labs do the lab work, non-dilutive money pays for it, a drug reaches patients, and
the money that comes back refills the foundation. Then we do it again.

**What I'm asking for right now is small.** Not a big raise. Just the one thing
blocking everything else: a wet lab willing to test a short list of candidates for
us, plus a little grant money to pay for that first real yes-or-no. Everything
after that depends on how that test goes.

---

## 2. What we're trying to do, and why it might work

**The mission is broad:** a world without antifungal resistance — drugs,
diagnostics, surveillance, stewardship. Fungal disease kills roughly 2.5 million
people a year, it gets diagnosed too slowly, and *Candida auris* is on the WHO's
critical priority list because it's beating the azole drugs.

**The work stays narrow on purpose:** we have exactly one proven asset — a tool
that triages molecules for CYP51, the target the azoles hit. Broad mission, narrow
execution. Diagnostics and the other tracks come later, one at a time, not all at
once.

Here's the chain we're betting on:

- The tool ranks a library and hands the top ~5% to a wet lab.
- The wet lab finds a real hit against resistant *C. auris* (this is the first
  honest yes/no — most computational hits die here).
- That hit becomes a lead compound whose IP the foundation owns.
- CROs take it through hit-to-lead, lead optimization, and the IND-enabling work.
- A pharma partner takes it through late trials and to market.
- The license money comes back to the foundation and pays for the next thing.

Every one of those arrows is a gate we can fail. The whole point of the plan is to
make each gate cheap to reach, honest to report, and non-fatal if the answer is
no.

---

## 3. An honest look at what we have

| Thing | Where it stands |
|---|---|
| The triage tool (geometry beats the docking score) | Built and honestly validated. AUC ~0.79 on a blind holdout; broad enrichment holds up; we retracted the flashier top-1% number when it didn't survive multiple seeds. |
| The self-checking gate (won't rank anything until it re-finds known drugs) | Built. Probably our most reusable, most defensible piece. |
| Pre-registration and hash-freezing | In place across ~25 analyses. The discipline is the moat, not some clever algorithm. |
| Candidate shortlists | Generated (see `work/RESULTS_repurposing.md`, `work/RESULTS_candidate_shortlist.md`). |
| Genomic early-warning track | Scaffold done, one blocker (the ERG11 re-caller isn't built; NCBI doesn't run an AMR pipeline on *C. auris*). |
| A wet-lab partner | None. This is the one thing blocking us. |
| A validated hit | None yet. |
| The legal entities | Not formed yet. We do have the website, the `.org` domain, the org GitHub, and the tool under a noncommercial license set up so the foundation can license it commercially later. |

The takeaway: we're before the hit. So the near-term plan is really two things —
get that first wet-lab answer cheaply, and set up the ownership structure now
while it's cheap and nobody's fighting over it. It's much easier to make something
incorruptible before there's money on the table.

---

## 4. Who owns what

Three layers, kept separate on purpose so that control and money don't sit in the
same hands.

```
  The Antifungal Resistance Project Foundation   (nonprofit — no owners)
        owns the mission, the charter, the IP, and the control shares
                              │  controls
                              ▼
  DrugCo   (a public benefit corporation)
        tiny team, owns the drug program, hires CROs for everything
                              │  licenses / partners
                              ▼
  CROs, clinical partners, and later a pharma company for late trials + selling
```

**The Foundation** is a nonprofit. It has no shareholders, so it can't be bought.
That's the deepest protection there is — you can't acquire something nobody owns.
It holds the control shares in DrugCo, the core IP, the brand, and the charter. If
it ever shuts down, everything goes to another mission-aligned nonprofit, never to
a person.

**DrugCo** is the company that actually does the drug development. It's a public
benefit corporation, which means its directors are legally required to weigh the
mission (getting antifungals to people, fighting resistance) against profit. We
need a for-profit here because clinical work means signing deals, taking milestone
payments, and maybe taking investment — a pure nonprofit is bad at that. DrugCo is
the hands; the foundation is the owner and the conscience.

**The control shares are the key trick.** The foundation holds a class of shares
with extra votes, so even after investors put in money and take most of the upside,
the foundation still controls the votes and can veto anything that betrays the
mission — a sale, a change of purpose, dropping our access commitments. This is
exactly how the Novo Nordisk Foundation controls Novo Nordisk while owning a
minority of the actual economics.

**The CROs** do all the wet-lab, chemistry, tox, manufacturing, and clinical work
under contracts that assign the IP back to DrugCo.

---

## 5. How we keep it incorruptible

The charter ([CHARTER.md](CHARTER.md)) is the real instrument; this is the
thinking behind it.

The Eric Ries point — from *The Startup Way* and his Long-Term Stock Exchange — is
that you have to build mission protection into the structure at the very start,
because it's almost impossible to bolt on later once money and incentives show up.
I agree, so we're doing it now.

We're borrowing from a few places that already pulled this off:

- **Novo Nordisk Foundation** — a foundation controlling a pharma company through
  special shares. This is our closest model.
- **Patagonia** — voting stock in a purpose trust, a locked mission, "Earth is our
  only shareholder."
- **Anthropic's Long-Term Benefit Trust** — an independent trust that holds a
  special share class insulated from money pressure.
- **Steward-ownership / the golden share** — separate control from cashing out;
  keep a veto share that blocks anything that breaks the mission.
- **Public benefit corporation** — a legal duty to weigh the mission against
  profit.
- **OpenAI** — we keep this one as a warning. A nonprofit-controls-a-for-profit
  setup can get strained when the capital demands get big, so we harden the parts
  that bend (see the risks section).

Five locks, spelled out in the charter:

1. **Ownership lock** — the top of the stack is a nonprofit with no owners. The
   mission isn't for sale.
2. **Control lock** — the foundation's golden share can veto a sale, a change of
   purpose, dropping access commitments, diluting the foundation below control, or
   changing the charter.
3. **Purpose lock** — the mission clause is hard to change: it takes a supermajority
   and the foundation's consent.
4. **Access lock** — real, binding affordability and access commitments that
   survive even after we license the drug to someone else.
5. **Transparency lock** — the honest reporting we already do (pre-registration,
   publishing the failures) becomes a governance duty, not just a habit.

And it's people, not only paper. A small board of mission trustees — independent,
term-limited, chosen for the mission, including scientific and patient-access
voices — holds the golden share's powers. I'm a trustee and the steward, not an
owner who can cash out or sell. That's the whole Ries idea: my leverage is vision
and stewardship, on purpose not capital I can flip.

---

## 6. The smallest possible drug company

The plan is a "virtual biotech." Our few employees own the decisions, the IP, and
the coordination. The CROs do the actual work. I'd like to reach IND with roughly
3 to 6 people.

| Work | Us (small) | Outside (CRO / partner) |
|---|---|---|
| The tool and the hypothesis | us | — |
| Chemistry and synthesis | design and oversight | CRO |
| Microbiology (MIC, resistant strains) | protocol and reading the data | academic/clinical lab or CRO |
| Selectivity, tox, ADME | study design | CRO |
| Manufacturing | — | CDMO |
| Regulatory / IND | a fractional consultant | regulatory CRO |
| Clinical trials | strategy and oversight | clinical CRO |
| Late trials + selling | license it out | pharma partner |

We hire in-house only as the gates clear: me as steward, the computational lead
who owns the tool, then a head of chemistry and biology who can turn a shortlist
into real molecules and read CRO data with a skeptical eye. Later a development/
regulatory person and someone to keep the CROs on schedule.

Staying small isn't just cheaper — it's the right shape here. The tool is a
handoff instrument by design, late trials happen through a partner anyway, and a
small team keeps the burn low. Low burn is itself a protection: it's what keeps us
from taking bad money to make payroll.

---

## 7. The path, gate by gate

Each stage has a way in, an honest yes/no to get out, and a point where we stop
and say so.

| # | Stage | To start | The real yes/no | Rough $ | Rough time |
|---|---|---|---|---|---|
| 0 | Structure + shortlist | — | Entities formed, charter signed, a defensible shortlist | $10–50k legal | 1–3 mo |
| 1 | First wet-lab answer | a wet-lab partner | a compound with real activity against resistant *C. auris* and OK selectivity | $25–150k | 3–9 mo |
| 2 | Hit-to-lead | a validated hit | a lead series with SAR, potency, early selectivity, clear IP | $0.5–3M | 12–24 mo |
| 3 | IND-enabling | a lead compound | clean tox/ADME, IND accepted | $5–15M | 12–24 mo |
| 4 | Phase I/II | an IND | human safety plus a signal | $20–80M | 2–4 yr |
| 5 | Phase III + market | a Phase II signal | approval, licensed to a pharma partner | $100M+ (partner) | 3–5 yr |

The discipline: we don't raise or spend Stage-3 money until Stage 1 says yes. Most
computational hits die at Stage 1, and we've planned for that. A "no" there sends
us back to the tool with a new library or a new target — not into pretending.

---

## 8. Money — non-dilutive first, control-preserving always

The rule is simple: use money that doesn't cost us ownership before money that
does, and when we can't avoid the second kind, structure it so the foundation
keeps control and the access promises stay in place.

**First choice — non-dilutive.** This is where I want to fund Stages 0 through 3
if at all possible:

- Government and public health: NIH/NIAID grants, BARDA, and the AMR-focused
  programs. SBIR/STTR is a great fit — it's non-dilutive, we keep the IP, and it
  pairs a small company with an academic lab, which also solves our wet-lab
  problem.
- AMR-specific funders: CARB-X (funds early AMR work including antifungals),
  GARDP, the Wellcome Trust, the Gates Foundation, the Novo Nordisk Foundation.
  This is mission-friendly money that doesn't want to own or flip us.
- Prizes and challenges in AMR and drug discovery.
- Academic and in-kind: university mycology labs giving us strains, assays, and
  co-authored grants; open-science collaborations; compute credits.
- Philanthropy: recoverable grants and program-related investments from health
  foundations — money that wants the mission, not control.

**Second choice — mission-aligned catalytic money** as a bridge: the AMR Action
Fund and impact investors, recoverable grants, revenue- or milestone-based
instruments that avoid handing over equity. Any equity here is non-voting or
economic-only.

**Last resort — investment money, only at Stage 3 and beyond, always
control-preserving.** Investors buy ordinary shares with real upside; the
foundation keeps the control shares. They get the money side; we keep the wheel.
The big one is the Phase III license: instead of raising and burning $100M+
ourselves, we license late trials and selling to a pharma partner and take
upfront + milestones + royalties back to the foundation.

Then it loops: license money and any distributions go to the foundation and fund
the next target, the diagnostics project, and the surveillance work. The mission
compounds instead of the cap table.

---

## 9. Money math, roughly

These are planning sketches. The real numbers depend entirely on how Stage 1 goes.

- **To get the first real answer:** about $25–150k plus a wet-lab partner. This is
  the only number that matters right now, and it's reachable almost entirely with
  non-dilutive money — one SBIR grant, one foundation grant, or an in-kind academic
  collaboration.
- **Lead through IND:** somewhere around $6–18M total, aiming for most of it
  non-dilutive, the rest catalytic.
- **Trials and market:** hundreds of millions — not on our books, carried by a
  licensee. Our return is upfront, milestones, and royalties.
- **Keeping the foundation running:** grant-funded at first, then an endowment
  seeded by the first license deal, plus dividends from the control shares if
  DrugCo ever turns a profit (again, the Novo pattern).

The measure of success isn't a valuation. It's a real antifungal that resistant
patients can afford, owned by a mission that can't be sold. The money is the fuel,
not the point, and the structure is built to keep it that way.

---

## 10. What could go wrong

| Risk | What we do about it | If it happens anyway |
|---|---|---|
| No wet-lab partner (today's blocker) | SBIR/STTR and CARB-X both pair you with academic labs; plus direct outreach | We stall at Stage 0 — and we say so. The tool stays a public good either way. |
| Stage-1 no (hits fail the assays) | We expect this and budget for it; test several candidates; fork to a new library or target | Retool, re-run, re-target, and publish the null. |
| Non-dilutive money is slow | Apply in parallel; keep burn near zero so slow isn't fatal | Stretch the timeline instead of taking money that costs us control. |
| Money pressure erodes the mission (the OpenAI failure) | Golden share + entrenched purpose + low burn + PBC duties | The veto blocks the sale or pivot; access promises survive licensing. |
| Too dependent on me | Term-limited trustees, everything documented, the tool is open and reproducible | The trustees carry on; the charter outlives any one person. |
| The noncommercial license annoys partners | Free for all academic and noncommercial use; commercial license funds the mission | If it blocks a real partnership, we revisit it — mission beats leverage. |
| Scope creep (diagnostics or FKS1 before a hit) | The charter and VISION both say: broad mission, narrow execution, one thing at a time | We say no and keep the second project sequential. |

---

## 11. The next twelve months

**Now to 3 months:**
1. Lock down the candidate shortlist from the tool.
2. Set up the ownership structure while it's cheap: incorporate the foundation,
   adopt the charter, reserve the DrugCo shell and the control/golden share
   classes — even before there's value in them.
3. Get low-cost or pro-bono legal review of the structure (a lot of firms will do
   this for mission-driven AMR work).
4. Line up wet-lab partners through SBIR/STTR and CARB-X applications that bundle a
   lab.

**3 to 9 months:**
5. Land one non-dilutive grant or in-kind collaboration.
6. Run MIC and selectivity on the top candidates against resistant *C. auris*.
7. Report the result honestly, yes or no, in `work/`.

**9 to 12 months:**
8. If yes: file provisional IP (assigned to the foundation/DrugCo), scope
   hit-to-lead, and start raising non-dilutive money for Stage 2.
9. If no: run the fork — new library or new target — and go back into Stage 1. The
   structure and the tool stick around either way.

---

## 12. The ask

- **Most important:** a wet-lab partner (a clinical or academic mycology lab) to
  run that first MIC and cytotoxicity panel. This one relationship unblocks the
  whole plan, and it's the natural bridge to the diagnostics project down the road.
- **Next:** about $25–150k in non-dilutive money to pay for that first answer —
  one SBIR Phase I, one foundation grant, or one prize.
- **Structural:** the go-ahead to incorporate the foundation and adopt the charter
  now, so the mission is locked before there's anything worth corrupting.

Written the way we run the science: the scaffold is real, the gaps are named, and
nothing here claims a win — scientific or financial — that we haven't earned yet.
