# A pre-validated 10-compound repurposing shortlist for *Candida* CYP51 — looking for a wet lab

**OpenAFR / The Antifungal Resistance Project** · 2026-08 · [repository + all data](../README.md)

We are looking for one mycology or clinical-microbiology lab willing to run **MIC panels on ten
already-approved or clinical-stage drugs** against *Candida albicans* and *C. auris*. Roughly a
week of bench time. We supply the compound list, the rationale, the controls, and every negative
result behind it. We have no financial interest in the outcome; everything is published either way.

---

## Why these ten and not another computational shortlist

Virtual screens produce ranked lists constantly and almost none of them come with evidence the
ranking means anything. Ours does, and it also comes with the evidence of where it stops working.

**1. The docking score was thrown away, because on this target it is anti-informative.**
On a blind, pre-registered held-out set — 7 azole antifungals hidden among 348 property-matched
decoys — AutoDock Vina's own affinity score ranked the real drugs **worse than random** (AUC
0.471). Sorting the *identical poses* by a mechanistic geometric criterion instead — how closely
the molecule brings a nitrogen to the heme iron, which is how azoles actually work — gave **AUC
0.794, permutation p = 0.0028**.

**2. The criterion re-discovers real antifungals in a library it was never told the answer to.**
Pointed blind at 2,776 drugs from the Broad Drug Repurposing Hub, it placed seven known azoles in
the top 3.2%: flutrimazole #1, ravuconazole #9, fosfluconazole #22, fluconazole #31, voriconazole
#57, luliconazole #73, croconazole #89.

**3. We know exactly how good it is not.** Four pre-registered follow-up tests failed, and we
published all four. Under 5-seed replication the top-1% enrichment collapses to zero — the
property-matched decoys reach the iron just as closely (2.69–2.84 Å) and just as axially (2.7–7.5°)
as the real drugs do. What survives replication is **broad** enrichment: AUC ~0.86, EF@5% ~8.5x.

**So the honest claim is narrow and, we hope, credible: the top ~5% of this ranking is enriched for
real iron-coordinating molecules by roughly 8-fold over random. It is not a prediction that any of
these ten kills a fungus.** We are asking you to test a shortlist that is better than random, not a
drug candidate.

## The shortlist

Ranked by the validated criterion out of 2,776 docked compounds. All are approved or have reached
clinical development, so they have existing human safety and PK data, and all are catalog items.

| rank | drug | usual indication / class | N–Fe (Å) | Vina (demoted) |
|---:|---|---|---:|---:|
| 2 | **verinurad** | URAT1 inhibitor (gout) | 2.417 | −7.3 |
| 3 | **vatalanib** | VEGFR tyrosine-kinase inhibitor | 2.479 | −11.0 |
| 4 | **flucloxacillin** | β-lactam antibiotic (isoxazolyl) | 2.558 | −7.8 |
| 5 | **taranabant** | CB1 inverse agonist | 2.586 | −10.5 |
| 6 | **lersivirine** | HIV NNRTI | 2.599 | −7.1 |
| 7 | **LDN-27219** | transglutaminase-2 inhibitor | 2.603 | −9.3 |
| 8 | **R-1479 (balapiravir)** | HCV polymerase inhibitor | 2.621 | −7.4 |
| 12 | **epacadostat** | IDO1 inhibitor (oncology) | 2.674 | −8.5 |
| 14 | **nolatrexed** | thymidylate synthase inhibitor | 2.692 | −9.4 |
| 19 | **apatinib** | VEGFR2 inhibitor | 2.711 | −10.7 |

For calibration, the real azoles in the same run sit at 2.37–2.84 Å. Several of these are
N-heterocyclic kinase inhibitors — chemistry that can coordinate a heme iron, which is exactly what
the criterion measures and exactly why coordination is necessary but not sufficient.

## What we would ask you to run

Nothing exotic — this is a standard panel, and we have tried to make the ask as small as possible.

- **Assay:** CLSI M27 broth microdilution, MIC at 24/48 h. Two-fold dilution series to 64 µg/mL.
- **Strains:** at minimum one azole-susceptible *C. albicans* (the criterion was validated on a
  *C. albicans* structure, 5TZ1 — this is the in-domain test) and one fluconazole-resistant
  *C. auris* (the mission-relevant stretch; CDC & FDA AR Isolate Bank supplies these panels free to
  qualifying labs). More clades if convenient, not required.
- **Controls:** fluconazole and voriconazole on the same plates.
- **Optional and, we suspect, the more interesting experiment:** checkerboard synergy with
  fluconazole. A weak CYP51 binder that potentiates an azole against a resistant isolate would be a
  more plausible outcome than standalone activity, and it is the same plates.

**Cost estimate:** the ten compounds are catalog items at roughly $100–300 each (MedChemExpress,
Selleck, Sigma); with plates and media the whole panel lands in the low four figures. We would
rather help you find funding for it than have you absorb it — see below.

## What we already know is wrong with this, stated up front

We would rather you hear these from us.

1. **Selectivity is not addressed at all, and this is the serious one.** The criterion measures heme
   coordination, and human CYPs have heme too. The same screen ranked **letrozole — a human
   aromatase inhibitor — at #23.** A compound that binds fungal CYP51 well and human CYP3A4 or
   CYP51 equally well is not a drug lead. Any hit here needs a human-cell cytotoxicity and CYP
   counter-screen immediately, and several of these compounds are plausible human CYP inhibitors on
   their face.
2. **The receptor is *C. albicans*, not *C. auris*, and it is rigid.** We ported it to the dominant
   *C. auris* Y132F resistance substitution and pre-registered the transfer test — **it failed** at
   the top of the ranking, though broad separation held (AUC 0.805). So the *auris* claim is weaker
   than the *albicans* claim, and we are not hiding that.
3. **This tool cannot, in principle, tell you a compound *beats* resistance — and we can show why.**
   We tried to validate resistance-awareness directly: assemble compounds with measured activity
   against azole-*resistant* *C. auris* and check whether the criterion ranks them above the azoles
   the mutant defeats. It doesn't work, for a structural reason rather than a data-gathering one. The
   compounds that retain activity against the resistant strain are all large, long-tailed molecules
   (VT-1598, 43 heavy atoms; opelconazole, 50; the oteseconazole/posaconazole class) — precisely the
   size range where a rigid receptor fails, because they need induced fit our model can't provide.
   The small molecules the tool docks reliably are exactly the ones resistance defeats. Resistance-
   retention is confounded with molecular size, so **the honest ceiling on this method is broad
   azole-like enrichment on a wild-type pocket — not resistance-breaker discovery.** We would rather
   state that plainly than let the *auris* framing imply more than the tool can carry.
4. **The per-molecule noise floor is real.** In the same 2,776-compound run, clotrimazole — an
   actual azole antifungal — ranked #2774 of 2776, because docking failed to find it a coordinating
   pose. A criterion that puts a real drug dead last on one molecule has genuine per-molecule
   variance. The *ranking* is enriched; individual ranks are not reliable.
5. **A single rigid pose says nothing about permeability, efflux, or metabolic stability**, all of
   which kill antifungal candidates routinely. *C. auris* in particular has formidable efflux.
6. **Nothing here has touched a living organism.** Most computational hits die at exactly this step.
   We expect that. It is why we want the experiment run rather than the list published as if it
   meant more.

## What we are pre-committing to

- **We publish the result either way, including "all ten inactive at ≤64 µg/mL."** That is written
  down before the experiment, in the same pre-registration style as everything above. A clean
  negative on ten repurposing candidates is a genuinely useful publication and we will treat it as
  the success case, not the failure case.
- **You own the experimental work and the authorship position that goes with it.** We are bringing
  a computational shortlist and a documentation trail, not a claim on your lab's output.
- **We will write the pre-registration with you before you run anything**, so the criteria for
  "this worked" are fixed in advance and neither of us can move them afterward.
- **No commercial interest.** The repository is PolyForm Noncommercial; there is no company, no
  patent position, and no funding contingent on a positive result.

## What we can bring besides the list

Full computational support at no cost: re-screening any library you care about under the frozen
protocol, docking any strain-specific ERG11 variant you are working with, and the pose files and
per-compound coordination evidence for anything on the list. If a specific resistance genotype from
your collection is of interest, we can model that pocket rather than the generic one.

---

**Everything above is reproducible.** All code, all five pre-registrations with their sha256 hashes
(one pass, four failures), the frozen protocol, the actives and decoy sets, and the complete
2,776-compound ranking are in the public repository. The full methodological writeup, including the
mechanistic explanation of why the criterion stops working at the top of the ranking, is
[here](../work/PREPRINT_geometry_ceiling.md).

**Contact:** `<CONTACT EMAIL>` — fill in before sending. A project-domain address reads better to an
academic recipient than a personal one.
