# Mechanism-anchored geometric rescoring on fungal CYP51: broad enrichment is robust, top-rank enrichment is not, and property-matched decoys explain why

**Draft manuscript — OpenAFR / The Antifungal Resistance Project.** Assembled 2026-08-11,
revised 2026-09-06, from ten pre-registered runs in this repository. Every number below is
traceable to a `work/RESULTS_*.md` file graded against a hash-frozen pre-registration; nothing
here is a new analysis, and no result was re-graded for this writeup. The dated update banners
above and §§4.3–4.5 fold in five later pre-registered results — a measured-inactive decoy set
(look #6), an ensemble receptor (looks #8–#9), a 2026-08-22 compound-curation finding that
bounds resistance claims (§4.4/§6.9), a 2026-09-02 coordinator-identity caveat (§4.5), and a
2026-09-06 power run that widened the primary gate to a blind N = 300 and retired the n = 7
caveat (look #10, §3.1/§6.1) — each graded against its own frozen hash and each re-grading
nothing that preceded it.


---

> **UPDATE 2026-08-23 — look #6, on an independent dataset: the decoy ceiling is REAL, not an
> artifact of decoy construction.** §4.1 below names experimentally confirmed inactives as "the
> most useful next contribution anyone could make to this benchmark", and §5's stopping rule
> permits exactly one thing: a new decoy construction. Both have now been done
> ([RESULTS_verified_inactives.md](RESULTS_verified_inactives.md), pre-registered and graded).
> The 7 held-out azoles were re-ranked, under the identical protocol, criterion and seed,
> against **279 compounds MEASURED not to inhibit *C. albicans*** (ChEMBL, consistent-inactivity
> rule). Result: **AUC 0.716 — the gate passes**, so the separation reported below is *not* a
> product of how the decoys were built, and Limitation 2 is discharged. But the tip is worse,
> not better: **7 measured inactives outrank the best real drug**, and of the 39 molecules
> inside the validated iron-bound band, **33 are compounds measured not to work and only 2 are
> actives**. A pre-specified subgroup split locates the effect precisely — AUC **0.810** against
> non-azole chemotypes versus **0.650** against measured-inactive *azole analogues*, below the
> bar. What the criterion mostly measures is chemotype, not potency. §4.1's bind is therefore a
> fact about the chemistry rather than about our decoy generator, and §4.4's "trustworthy
> top-~5% shortlist" should be read with the added caveat that such a shortlist will be crowded
> with azole-like molecules that have already been tried and failed.

---

> **UPDATE 2026-08-26 — the ceiling holds under an ensemble receptor, closing §4.3's leading
> escape route.** §4.3 names "induced-fit or ensemble receptor treatment" as the most promising
> signal a rigid single-conformer pose omits. It has now been tested. An exhaustive 3-conformer
> experimental *C. albicans* CYP51 ensemble (5TZ1 + 5FSA + 5V5Z, hash-pinned, mode C unchanged,
> aggregated by ensemble-minimum) was pre-registered and graded twice on the within-azole
> question. On the 7 held-out azoles it **passed** (AUC 0.750, look #8,
> [RESULTS_ensemble.md](RESULTS_ensemble.md)) — but with a wide n=7 CI (0.683–0.825) and a
> pre-committed demand for confirmation on an independent active set. That confirmation (look #9,
> [RESULTS_ensemble_confirm.md](RESULTS_ensemble_confirm.md), pre-reg sha `7754050e…`) ran on an
> **independent N=200** sample of ChEMBL measured-active azoles: **ensemble-min mode C AUC 0.682 <
> 0.70, CI 0.647–0.715 — FAIL.** The ensemble adds a genuine lift (5TZ1 alone 0.638 → ensemble
> 0.682, permutation p<0.0001) but lands below the usefulness bar, and the n=7 pass did not
> replicate. The within-azole ranking line is now closed on **both** the rigid and the ensemble
> receptor — a property of the question, not of any one conformer. The novel-chemotype guard is
> unaffected (AUC 0.823 vs 0.810). See revised §4.3 and ledger rows 7–9.

---

> **UPDATE 2026-09-02 — which nitrogen? A criterion-identity caveat on novel libraries.** The
> validated criterion (`min_nitrogen_iron_distance`) measures the closest approach of the nearest
> ligand nitrogen of *any* type. On the azole controls that is always the aromatic ring N the
> mechanism was validated on (6/6 in-mechanism, [RESULTS_coordinator.md](RESULTS_coordinator.md)).
> On a *novel* library it need not be: two of the strongest, cleanly-selective repurposing hits
> (taranabant, lersivirine) reach the iron through a **nitrile** N whose nearest aromatic ring N is
> 6.87 / 4.05 Å away — out-of-mechanism coordination the method never validated. The principled fix
> — rank on the aromatic ring N only — was pre-registered and tested
> ([RESULTS_aromatic_criterion.md](RESULTS_aromatic_criterion.md)): it **lowers** the validated gate
> AUC from 0.794 to 0.729, because a real held-out active, **luliconazole, itself coordinates via
> its nitrile (2.76 Å) not its imidazole (5.28 Å)**. A non-aromatic coordinator is therefore
> *unvalidated, not disproven*; the frozen criterion stands, and coordinator identity is reported
> as a per-candidate flag. See new §4.5.
> *(The PDF in this directory predates these updates.)*

---

## Abstract

Docking scoring functions rank known azole antifungals against fungal CYP51 **worse than
random**. A mechanism-anchored geometric criterion — the minimum distance from a ligand
nitrogen to the heme iron — ranks the same poses at AUC 0.794 on a pre-registered, blind
held-out set (permutation p = 0.0028), and reproduces this on a 2,776-compound unlabeled
library by pulling seven known azoles into the top 3.2%. That part holds up.

The same single-search run reported EF@1% = 12.68x. It does not survive replication. Under
five-seed sampling, EF@1% collapses to 0.00x under both a best-pose (min) and a
pool-size-neutral (median) aggregator, and again under an orthogonal *axial-direction*
feature — while AUC *improves* to 0.862. We show the cause is not estimator noise and not the
choice of geometric feature: **property-matched decoys coordinate the heme iron both closely
(2.69–2.84 Å) and on-axis (2.7°–7.5°), indistinguishably from real azoles**, so six decoys
outrank the best true active. We report two further findings: (i) transfer to a *C. auris*
Y132F resistance pocket preserves broad separation (AUC 0.805) while eliminating top-rank
enrichment via a sub-0.1 Å reshuffle; and (ii) at typical held-out sample sizes (7 actives,
355 molecules) EF@1% is arithmetically a **Bernoulli variable** taking only the values 0.00x
or 12.68x — it has no resolution and should not be used as a pass/fail criterion, including
by us.

The leading escape route — an ensemble of experimental receptor conformers — was
subsequently pre-registered and tested: it yields a genuine but sub-threshold lift (AUC 0.638
→ 0.682) that fails to clear the 0.70 bar and does not replicate on an independent 200-azole
set, so the top-rank ceiling stands on an ensemble receptor as well as a rigid one.

We conclude that mechanism-based rescoring buys *broad* enrichment — a trustworthy top-~5%
shortlist — and that a mechanism-matched decoy set structurally bounds what any
mechanism-based criterion can achieve at the extreme top. All pre-registrations,
including those that failed, are published with their hashes.

---

## 1. Background

Structure-based virtual screening rests on a scoring function whose correlation with real
affinity is, on most targets, weak. A standard response is to rescore poses against a
*mechanistic* criterion instead — a pharmacophore, a metal-coordination constraint — and rank
on that. The idea is decades old and is not the contribution here.

Fungal CYP51 (lanosterol 14α-demethylase) is an unusually clean test bed for it. The azole
antifungals inhibit by a single, well-defined structural event: an azole-ring nitrogen
coordinates the heme iron as the sixth axial ligand. The mechanism reduces to one measurable
geometric quantity, which removes the usual arbitrariness in choosing a pharmacophore. The
target also matters clinically: *Candida auris* is on the WHO fungal priority pathogens list
and defeats azoles largely through mutations in this enzyme.

The question this work set out to answer was not "does mechanism-based rescoring beat the
scoring function" — it does — but **how far does it carry, and where exactly does it stop.**

## 2. Methods

### 2.1 System

Receptor: PDB **5TZ1**, *C. albicans* CYP51, chain A plus heme `HEM A 601`, rigid. Docking:
AutoDock Vina, box 26 Å centred (70.61, 66.28, 4.18), `exhaustiveness 32`, `num_modes 20`,
seed 42. Every parameter is read at both docking and grading time from a single hash-frozen
`protocol.yaml` (sha256 `6eca16b1…`), whose integrity is verified by the grading script at run
time, so what was docked and what was graded cannot silently diverge.

### 2.2 Actives and decoys

Actives are known azole antifungals; the primary evaluation uses **7 held out** from all
criterion development (efinaconazole, luliconazole, fenticonazole, bifonazole, sertaconazole,
oxiconazole, butoconazole).

Decoys were generated from ChEMBL (`scripts/make_decoys.py`), 50 per active, under three
constraints chosen to defeat two specific ways this kind of benchmark cheats:

| constraint | value | what it prevents |
|---|---|---|
| property-matched | MW ±25, cLogP ±1.5, rotatable bonds ±2, HBD ±1, HBA ±2 | Docking score correlates with heavy-atom count at **r = −0.91** on this system; unmatched decoys let the pipeline "enrich" by preferring heavy molecules. |
| ≥ 1 aromatic nitrogen | required | Without it, a nitrogen-to-iron criterion would separate perfectly while only asking *"does this molecule contain a nitrogen?"* — circular. |
| Morgan(r=2) Tanimoto < 0.35 vs **every** active | required | Same physical profile, different skeleton — the basis for presuming non-binding. |

348 decoys survived to the graded pool (355 molecules total; 2 docking failures, both decoys).
**Decoys are presumed inactive, not experimentally confirmed** — see §6.

### 2.3 Ranking criteria

- **Mode A — Vina score** (the control).
- **Mode C — iron-approach distance:** per molecule, the minimum over poses of the distance
  from the nitrogen nearest the iron to the heme Fe. This is the primary criterion.
- **Mode D —** combined rank of A and C.
- **Mode E — axial (added later):** `S = |v| · (1 + sin θ)`, where `v = N − Fe` and `θ` is the
  acute angle between `v` and the heme normal `n = normalize((NC−NA) × (ND−NB))`, fixed once
  from the receptor as `(+0.6048, −0.7591, −0.2407)`. λ = 1, untuned. Rewards approaches that
  are both close *and* on-axis.

### 2.4 Metrics and the pre-registered bar

AUC, EF@1%, EF@5%, EF@10%, BEDROC(α=20), plus a 20,000-shuffle permutation test and a
bootstrap CI over actives. The pass bar, declared before any held-out result was seen and never
lowered, was **AUC ≥ 0.70 AND EF@1% ≥ 5.0x**, evaluated on the primary mode only.

### 2.5 Pre-registration procedure

Each run has a pre-registration file committed *before* the run, hashed, and the hash verified
unmodified by the grading script at grading time. Each pre-registration declares the criterion,
the aggregation, the metrics, the bar, **and the conclusion to be drawn on failure** — so a
negative cannot be reinterpreted after the fact. Where a run reuses existing poses, it says so
and grades without new docking.

## 3. Results

### 3.1 Geometry beats the scoring function (pre-registered, blind, held out)

On 7 held-out azoles among 348 property-matched decoys, sorting the *identical* poses by
iron-approach distance rather than by Vina score inverts the outcome:

| mode | AUC | EF@1% | EF@5% | BEDROC | role |
|---|---:|---:|---:|---:|---|
| **C — iron-approach distance** | **0.794** | 12.68x | 5.63x | 0.325 | primary — **PASS** |
| A — Vina score alone | 0.471 | 0.00x | 0.00x | 0.002 | control — **worse than random** |
| D — combined rank | 0.691 | 0.00x | 2.82x | 0.126 | secondary |

Permutation test (20,000 label shuffles): **p = 0.0028**. Bootstrap 95% CI over actives: AUC
0.590–0.946. Held-out active ranks (of 355): 3, 8, 29, 43, 49, 124, 274.

That n = 7 CI lower bound (0.590) sat *below* the 0.70 bar, so this run alone could not exclude
a true AUC beneath it. A pre-registered power run (look #10, [RESULTS_active_power.md](RESULTS_active_power.md),
pre-reg sha `9327563a…`) later re-ran the *identical* gate — same receptor, criterion, decoy
pool, protocol and seed — with only the active side widened, to a blind **N = 300** sample of
measured-active azoles: **AUC 0.750, bootstrap 95% CI 0.717–0.782**, permutation p = 5×10⁻⁵.
The point estimate is a little below the n = 7 value, but the interval is ~4× tighter and its
lower bound clears the bar, so a sub-0.70 true AUC is now excluded (see Limitation 1).

This is the durable result. On this target, the docking program's own affinity estimate is
anti-informative, and a single mechanistic distance recovers real drugs it has never seen.

A prior run (`RESULTS_validation_gate.md`) *failed* this gate — a 3.0 Å hard metal cutoff
discarded 4 of 8 actives, and enlarging the box from 26 Å to 30 Å without scaling search effort
degraded sampling. That failure is published alongside the pass.

### 3.2 The top-1% enrichment does not survive replication

The run in §3.1 docked each molecule once. Re-docking the same 355 molecules under five seeds
{42, 1, 2, 3, 4} and aggregating gives:

| aggregator | AUC | **EF@1%** | EF@5% | EF@10% | BEDROC | actives in top 15 |
|---|---:|---:|---:|---:|---:|---:|
| single search (seed 42) | 0.794 | **12.68x** | 5.63x | 4.23x | 0.325 | 2/7 |
| consensus — min over 5 seeds | 0.853 | **0.00x** | 8.45x | 5.63x | 0.389 | 3/7 |
| reliability — median over 5 seeds | 0.830 | **0.00x** | 5.63x | 4.23x | 0.310 | 2/7 |

Every *broad* measure improves or holds under more thorough sampling. The top-1% figure goes to
zero and stays there. The two aggregators bracket the plausible options — `min` takes the best
pose found, `median` is pool-size-neutral and robust — and pre-registration bounded the search
to exactly these two, with the "this is a ceiling, not noise" conclusion pre-committed to the
failure branch. Both failed.

The `median` result is what distinguishes the two explanations. If the collapse were extreme-value
bias — five searches giving decoys five chances to get lucky — a pool-size-neutral robust
estimator would recover the signal. It does not. The signal is not there to recover.

### 3.3 It is not the distance feature either: decoys coordinate axially

The one remaining explanation was that N–Fe distance is the wrong number. Azole inhibition is
*axial*; a decoy might reach the iron equatorially and be a geometric impostor. We pre-registered
exactly one orthogonal feature — approach direction — and graded it on the same 5-seed poses with
no new docking:

| ranking key | AUC | **EF@1%** | EF@5% | EF@10% | BEDROC |
|---|---:|---:|---:|---:|---:|
| N–Fe distance, median over 5 seeds | 0.830 | 0.00x | 5.63x | 4.23x | 0.310 |
| **axial `S`, median over 5 seeds** | **0.862** | **0.00x** | **8.45x** | **5.63x** | 0.340 |
| axial angle θ **alone** | 0.861 | 0.00x | 5.63x | — | 0.240 |
| seed-42 single-search `S` | 0.819 | 0.00x | 8.45x | — | 0.320 |

AUC 0.862 is the best of any look taken. EF@1% is zero. Critically, adding the direction term to
the *original seed-42 poses* — the ones that produced 12.68x — also gives 0.00x: the direction
term pushes decoys up, not down.

The reason is visible directly in the top of the ranking:

| rank | molecule | median `S` | median θ | median N–Fe (Å) | |
|---:|---|---:|---:|---:|---|
| 1 | decoy_CHEMBL6328_for_efinaconazole | 2.881 | 3.4° | 2.702 | decoy |
| 2 | decoy_CHEMBL275096_for_oxiconazole | 2.928 | 2.7° | 2.790 | decoy |
| 3 | decoy_CHEMBL10329_for_oxiconazole | 3.017 | 3.9° | 2.793 | decoy |
| 4 | decoy_CHEMBL447532_for_butoconazole | 3.042 | 7.5° | 2.691 | decoy |
| 5 | decoy_CHEMBL10330_for_sertaconazole | 3.079 | 4.3° | 2.842 | decoy |
| 6 | decoy_CHEMBL10547_for_oxiconazole | 3.094 | 4.3° | 2.827 | decoy |
| **7** | **luliconazole** | **3.133** | **6.2°** | **2.781** | **ACTIVE** |

Six decoys outrank the best real drug, and they are **not off-axis impostors**: their approach
angles (2.7°–7.5°) are as small as or smaller than luliconazole's (6.2°). They are doing
azole-like axial iron coordination. No weighting of distance and direction can demote them,
because on both axes they are doing the thing the mechanism says to do.

**This is the central finding of the paper.** The top-1% ceiling is not an aggregation artifact
and not a missing feature. It is a statement about the pose geometry available to a rigid-receptor
docking run: *closest-nitrogen distance and its approach angle saturate against a mechanism-matched
decoy set.*

### 3.4 Transfer to a resistant pocket: broad separation survives, the tip does not

The mission-relevant question is whether a wild-type-validated criterion still works once the
pocket carries the dominant *C. auris* resistance substitution **Y132F**. Controlled single
variable: same 7 actives, same decoys, same frozen protocol; only the residue-132 side chain
changed (Tyr→Phe by deterministic deletion of the para-hydroxyl).

| mode | AUC | EF@1% | EF@5% | BEDROC | role |
|---|---:|---:|---:|---:|---|
| **C — iron-approach distance** | **0.805** | **0.00x** | 5.63x | 0.235 | primary — **FAIL** |
| A — Vina score | 0.411 | 0.00x | 0.00x | 0.002 | control |
| D — combined rank | 0.677 | 0.00x | 1.41x | 0.040 | secondary |

AUC is essentially unchanged from wild type (0.805 vs 0.794). Actives rank 7, 16, 42, 46, 50, 79,
264 of 355. The failure is confined to the tip: **the top six slots are all decoys, reaching the
iron at 2.62–2.72 Å against the best real drug's 2.722 Å.** A sub-0.1 Å reshuffle among near-tied
molecules is enough to zero the metric — the same phenomenon as §3.3, now triggered by a single
atom's deletion.

### 3.5 External corroboration of the broad signal on an unlabeled library

Pointed at 2,776 docked compounds from the Broad Drug Repurposing Hub — screened blind, never
told which are antifungals — the criterion placed seven known azoles in the top 3.2%:

flutrimazole #1 (2.371 Å), ravuconazole #9, fosfluconazole #22, fluconazole #31, voriconazole #57,
luliconazole #73, croconazole #89. The heme-coordinating human CYP inhibitor **letrozole** ranked #23.

This is an embedded positive control on a genuinely novel library and it behaves consistently with
§3.1–3.3: strong broad enrichment. It also exposes the per-molecule noise floor — **clotrimazole,
a real azole antifungal, ranked #2774 of 2776** (best pose 17.97 Å from the iron), because docking
simply failed to find a coordinating pose for it in 20 modes.

Letrozole's rank is worth stating plainly as a limitation rather than a success: the criterion
measures heme coordination, and human CYPs have heme too. **Nothing in this method addresses
fungal-vs-human selectivity**, which is the hard problem in antifungals specifically.

### 3.6 EF@1% has no resolution at this sample size

With 7 actives among 355 molecules, the top 1% is the top 4 slots. Enrichment factor is
`(a/k)/(A/N)`, so:

| actives in top 4 | EF@1% |
|---:|---:|
| 0 | **0.00x** |
| 1 | **12.68x** |
| 2 | 25.36x |

The metric is, in practice, a **Bernoulli variable**: every EF@1% we ever reported was either 0.00x
or 12.68x, and the entire "12.68x → 0.00x collapse" is one molecule crossing rank 4. EF@5% is
similarly quantized in steps of 2.82x (2 actives in the top 18 → 5.63x; 3 → 8.45x).

Verified directly against the implementation this project grades with
(`rdkit.ML.Scoring.CalcEnrichment` via `openafr/scoring.py`), not derived analytically:

```python
from openafr.scoring import metrics
N, A = 355, 7
for k in range(3):                                  # k actives in the top 4
    flags = [1]*k + [0]*(4-k) + [1]*(A-k) + [0]*(N-4-(A-k))
    print(k, metrics(flags)[1][0])                  # -> 0.0, 12.68, 25.36
```

We flag this against our own design: **the pre-registered gate should not have been conditioned on
EF@1%**, because at n = 7 the metric cannot take a value between "fail" and "3.5x over the bar."
This does not rescue the negative — §3.3's angle analysis establishes the ceiling mechanistically
and independently of any enrichment metric — but it does mean that the *pass* in §3.1 was
over-read at the time, and the subsequent "failures" were partly an instrument with two settings.
The robust, resolved statements are AUC, EF@5%, and BEDROC.

## 4. Discussion

### 4.1 The decoy-construction bind

Benchmark decoy sets are usually criticised for being **too easy** — property-matched-but-trivially-separable
decoys inflate reported performance, the well-known concern with DUD-E-style sets in
machine-learning docking evaluations. Our result is the mirror image and, we think, the more
instructive one.

To avoid circularity, we *required* every decoy to be a plausible iron-coordinator (≥1 aromatic
nitrogen) and to match the actives on size and lipophilicity. That is the methodologically correct
choice. But it means the decoy set was constructed to be indistinguishable from the actives *on
precisely the axis the criterion measures*. §3.3 shows the consequence: at the tip, they are.

This generalizes into a constraint on mechanism-based rescoring:

> A mechanism-anchored criterion can separate actives from decoys only to the extent the decoys do
> not share the mechanism. Match decoys on the mechanistic feature, and the criterion's measured
> ceiling *is* the benchmark you built.

The bind is real and has no computational escape. Non-matched decoys give an inflated, circular
result; mechanism-matched decoys give a criterion that saturates at the top. The only clean exit is
**experimentally confirmed inactives** — compounds measured not to inhibit — rather than
computationally presumed ones. We do not have those; obtaining them for CYP51 is, in our view, the
most useful next contribution anyone could make to this benchmark.

### 4.2 Practical recommendations

1. **Do not report EF@1% from a single docking run.** Ours moved from 12.68x to 0.00x on
   re-sampling alone, with no change to the criterion. Report multi-seed aggregates.
2. **Do not gate on EF@1% at small active counts.** Compute the achievable values first; if the
   metric is Bernoulli, it is not a gate (§3.6).
3. **Report AUC, EF@5%, and BEDROC together**, all under replication. On this system those were
   stable across four independent looks (AUC 0.79–0.86) while the tip was not.
4. **Test transfer explicitly.** A criterion validated on a wild-type structure should not be
   assumed to hold on the resistant variant that motivates the work (§3.4). Ours did not, at the
   tip.
5. **State what the mechanism does not cover.** Heme coordination does not distinguish fungal from
   human CYP (§3.5).
6. **Check assay precedent before proposing a shortlist for the bench.** Novel-to-a-fingerprint is
   not never-tested: a molecule may already have a deposited ChEMBL/PubChem record showing it
   inactive against fungal CYP51, and because failing measurements are under-published, a geometry
   rank can quietly re-run a lost experiment. We wired an exact-match plus near-neighbour precedent
   check as the last stage of the pipeline (`scripts/check_precedent.py`) and ran it over the top 100
   of the repurposing screen (§3.5) as a worked example. It validated in-run — the only enzyme-level
   flags were the known azole controls (fluconazole, voriconazole), never a novel hit — and demoted
   seven top-ranked repurposing candidates that a whole-cell assay had already found inactive, while
   the near-neighbour probe against the verified-inactive azole analogues returned nothing at
   Tanimoto ≥ 0.85, consistent with the hits being genuinely off the azole scaffold. Absence of a
   record is reported as `no-precedent`, never "untested." The full read-out is in the repository's
   candidate-shortlist result.

### 4.3 What could break the ceiling — and the one route we tested

Signals a rigid-receptor, single-conformer Vina pose does not contain, in rough order of promise:
explicit iron–nitrogen coordination chemistry at a QM or semi-empirical level; induced-fit or
ensemble receptor treatment; desolvation.

We tested the second, because it was the most promising and the most tractable. An exhaustive
ensemble of the three experimental *C. albicans* CYP51 crystal structures (5TZ1 + 5FSA + 5V5Z,
hash-pinned, mode C unchanged, aggregated by the pre-committed ensemble-minimum) was pre-registered
and graded on the within-*class* question — separating measured-active azoles from measured-*inactive*
azole analogues, the pool where a single rigid conformer scored only AUC 0.650. On the 7 held-out
azoles the ensemble **passed** at 0.750 (look #8, `RESULTS_ensemble.md`), but with an n=7 CI of
0.683–0.825 and a pre-committed demand for confirmation on an independent active set. That
confirmation ran on an independent N=200 sample of ChEMBL measured-active azoles and **failed at AUC
0.682 < 0.70** (look #9, `RESULTS_ensemble_confirm.md`, tight CI 0.647–0.715). The ensemble delivers
a real but sub-threshold lift — 0.638 (5TZ1 alone) → 0.682, permutation p<0.0001 — and the conformer
carrying it is itself sample-dependent (5FSA in look #8, 5V5Z in look #9), the signature of a
small-sample interaction rather than a robust mechanism. **The within-azole tip therefore does not
survive an ensemble receptor any more than it survives a rigid one**, and §5's stopping rule now
closes the line on both. The remaining signals (QM coordination chemistry, desolvation) are
undecidable from these data, and we make no claim any of them would work. All of this concerns
*within-class* ranking only; the novel-chemotype guard is unaffected by the ensemble (AUC 0.823 vs
the 0.810 single-conformer comparator), consistent with §4.4.

### 4.4 What the method is good for

Unchanged from what the data supports: **a trustworthy top-~5% shortlist for experimental triage**
(AUC ~0.86, EF@5% ~8.5x), not razor-top ranking. That is a real and useful thing — it makes a
wet-lab yes/no more likely than random and cheaper to reach — and it is a smaller claim than the
one the first run appeared to license.

One boundary on that claim is worth stating precisely, because §3.4's "broad separation survives
Y132F" invites over-reading: **this method cannot, in principle, tell you a compound *beats*
resistance.** We tried to validate resistance-awareness directly — assemble compounds with measured
activity against azole-*resistant* *C. auris* and check whether the criterion ranks them above the
azoles the mutant defeats (a 2026-08-22 curation pass, `data/ligands/PROVENANCE.md`). It cannot be
done with this tool, for a structural rather than a data-gathering reason. The compounds that retain
activity against the resistant strain are all large, long-tailed molecules (VT-1598, 43 heavy atoms;
opelconazole, 50; the oteseconazole/posaconazole class) — over the rigid receptor's declared ≤45
heavy-atom applicability domain (§6.4), because they need the induced fit a single rigid conformer
cannot provide. The small molecules the tool docks reliably are exactly the ones resistance defeats.
So **resistance-retention is confounded with over-domain molecular size**, and re-discovering
Y132F-defeated azoles on the mutant pocket (§3.4) demonstrates robustness of *broad enrichment*, not
resistance-breaker discovery. The honest ceiling is broad azole-like enrichment on a wild-type-like
pocket. Lifting it would need a flexible-receptor / induced-fit protocol that widens the domain —
listed in §4.3, not claimed here.

### 4.5 Which nitrogen coordinates — a criterion-identity caveat

The criterion ranks on `min_nitrogen_iron_distance`: the closest approach to the heme iron of the
nearest ligand nitrogen *of any type*. On the azole actives this is unambiguous — the mechanism is an
aromatic triazole/imidazole ring N donating to the iron, and a purely topological classifier (a
nitrile is the only nitrogen with a single sub-1.25 Å bond; a ring N carries two ~1.34 Å bonds)
confirms **6/6 controls coordinate through that ring N** (`RESULTS_coordinator.md`, pre-reg sha
`b595de48…`). On a novel library the identity of the coordinating nitrogen is *not* guaranteed.
Re-examining the strongest hits of a repurposing focus set by coordinator identity, two of the four
strongest cleanly-selective candidates — taranabant and lersivirine — reach the iron through a
terminal **nitrile** nitrogen, with their nearest aromatic ring N 6.87 and 4.05 Å away. Their tight
ranks are coordination through chemistry the method was never validated on.

The obvious fix is to restrict the criterion to the aromatic ring nitrogen. We pre-registered and
tested it (`RESULTS_aromatic_criterion.md`, sha `9a164558…`), and it does not survive contact with
the validation set: on the re-docked run-2 held-out gate the aromatic-only criterion scores **AUC
0.729 versus the any-N 0.794** — a 0.065 drop — because a genuine held-out active, **luliconazole,
coordinates the iron through its own nitrile (2.76 Å) while its imidazole sits 5.28 Å away**.
Enforcing the mechanism demotes a real drug. The lesson is not that non-aromatic coordination is
spurious but that it is *unvalidated, not disproven*: at least one working azole uses it. We
therefore leave the frozen criterion unchanged and report coordinator identity as a per-candidate
confidence flag — a molecule ranked on a nitrile is a weaker hypothesis than one ranked on a ring
nitrogen, not a rejected one. The general point stands for anyone deploying a distance-to-metal
criterion on a novel library: **the distance does not tell you which atom achieved it, and on a
novel scaffold that atom may not be the one the mechanism validated on.**

## 5. Pre-registration and multiple-comparison ledger

Every look ever taken at the held-out data, in order, with its committed outcome:

| # | run | question | pre-reg sha256 | outcome |
|---:|---|---|---|---|
| 0 | `RESULTS_validation_gate.md` | initial gate, 3.0 Å hard cutoff | — | **FAIL** (published; blocked the screen) |
| 1 | `RESULTS_run2_holdout.md` | iron-approach distance, blind holdout | `7c89d78c…` | **PASS** |
| 2 | `RESULTS_auris_Y132F.md` | does it transfer to the resistant pocket? | `689de2fd…` | **FAIL** |
| 3 | `RESULTS_consensus.md` | 5-seed consensus (min) | `1c6ea042…` | **FAIL** |
| 4 | `RESULTS_reliability.md` | pool-size-neutral median — noise or ceiling? | `f8aa1ff9…` | **FAIL** |
| 5 | `RESULTS_axial.md` | one orthogonal feature: approach direction | `3ff59c29…` | **FAIL** |
| 6 | `RESULTS_verified_inactives.md` | **new dataset**: same criterion vs 279 *measured* inactives | `d2926f77…` | **PASS** (AUC 0.716) |
| 7 | `RESULTS_channel_engagement.md` | one orthogonal feature (mode F) on the azole-analogue pool | `62592a10…` | **FAIL** (AUC 0.590) |
| 8 | `RESULTS_ensemble.md` | 3-conformer ensemble receptor, within-class, n=7 held-out | `8041b38e…` | **PASS** (AUC 0.750) — *superseded by 9* |
| 9 | `RESULTS_ensemble_confirm.md` | ensemble receptor, **independent N=200** azole actives | `7754050e…` | **FAIL** (AUC 0.682) |
| 10 | `RESULTS_active_power.md` | run-1 gate, actives widened to **blind N=300** (power run) | `9327563a…` | **PASS** (AUC 0.750, CI 0.717–0.782) |

Four looks were taken at the original wild-type poses (1, 3, 4, 5). Each pre-registration bounded the
number of subsequent attempts and pre-committed the failure conclusion; the reliability
pre-registration foreclosed further *aggregators*, and the axial pre-registration foreclosed further
*features* on this dataset. Once that stopping rule was in force, no further criterion could be tried
on those poses: any further attempt required an independent dataset — new actives, a new decoy
construction, or a new receptor — not another feature on the same one.

Looks 6–9 are exactly those independent tests, each pre-registered before any of its molecules were
docked. Look 6 rebuilt the decoys from *measurement* rather than assumption (279 measured inactives)
and passed broadly (AUC 0.716) while confirming the tip is worse, not better. Look 7 spent the one
permitted orthogonal feature (channel engagement) and failed (0.590). Looks 8–9 changed the
*receptor*, not the criterion: an experimental 3-conformer ensemble, which passed on 7 held-out
azoles (0.750) but failed to replicate on an independent 200-azole sample (0.682 < 0.70). With the
within-class line now closed on a rigid receptor, an orthogonal feature, *and* an ensemble receptor,
the settled conclusion is that top-rank separation of working from failed azoles is a property of the
question, not of any feature, aggregator, or conformer set. What remains open is only what §4.3
lists as undecidable from these data.

The `RESULTS_axial.md` run also included a harness check: the reliability-distance slice reproduced
the previously published reliability numbers (0.830 / 0.00x / 5.63x / 4.23x / 0.310) bit-for-bit.

## 6. Limitations

1. **The primary gate ran on n = 7 held-out actives** (bootstrap 95% CI 0.590–0.946, lower bound
   below the 0.70 bar). This was retired by a pre-registered power run (look #10,
   `RESULTS_active_power.md`): the *identical* gate with the active side widened to a blind
   **N = 300** sample of measured-active azoles gives **AUC 0.750, 95% CI 0.717–0.782** —
   the lower bound now clears the bar, so a true AUC beneath it is excluded. The N = 300 point
   estimate (0.750) is marginally below the n = 7 value (0.794); the honest headline is the tight
   interval, not the small-sample point estimate. (The *within-azole* n = 7 → N = 200 widening is a
   separate question and went the other way — see Limitation 3 and look #9.)
2. **Decoys are presumed inactive, not verified.** Any true binder among the 348 depresses measured
   performance; any systematic selection bias inflates it. This is the single largest caveat and the
   subject of §4.1.
3. **One decoy set for the rigid-pose looks.** The rigid-receptor ceiling is demonstrated *on that
   set*; look #6 re-tested it against an independently constructed measured-inactive set (AUC 0.716)
   and looks #8–#9 against an independent 200-azole active set, so the conclusion no longer rests on
   a single decoy construction.
4. **Primarily one rigid receptor**, *C. albicans* 5TZ1, solved with the short ligand VT-1161. No
   induced fit for the main analyses; the posaconazole/itraconazole long-tail class is a known blind
   spot and the declared applicability domain is ≤ 45 heavy atoms. An experimental 3-conformer
   ensemble (5TZ1+5FSA+5V5Z) was tested for the within-class question (§4.3) and did not lift the tip
   above the bar; full flexible/induced-fit treatment remains untested.
5. **The Y132F model is a rigid single-atom deletion.** Clinical Y132F also involves conformational
   change, H-bond network rearrangement, and co-occurring efflux and expression changes. K143R
   (clade I) and F126L (clade III) were not tested.
6. **One docking program** (Vina), one conformer per molecule, SMILES protonation as supplied.
7. **The heme normal is an idealization** estimated from four atoms of a rigid porphyrin.
8. **Nothing here has been tested against a living fungus.** Every ranked molecule is a hypothesis.
9. **Resistance-retention is confounded with molecular size, so no resistance-breaker claim is
   possible** (2026-08-22 curation, `data/ligands/PROVENANCE.md`; see §4.4). The verifiable pool of
   compounds active against azole-resistant *C. auris* is dominated by molecules above the ≤45
   heavy-atom applicability domain (§6.4), where the rigid receptor fails; the compounds it docks
   reliably are the ones resistance defeats. A validation of resistance-*breaking* is therefore not
   achievable with this rigid tool, and the method's ceiling is broad azole-like enrichment only.

## 7. Data and code availability

All code, all nine pre-registrations with their hashes, the frozen `protocol.yaml`, the actives and
decoy SMILES, the full 2,776-compound ranking, and every `RESULTS_*.md` are in this repository.
Docking intermediates (~180 MB) are regenerated rather than versioned; each results file carries an
exact reproduce recipe. Licensed PolyForm Noncommercial 1.0.0.

---

### Note on framing

The ingredients here are not new. Metal-coordination and pharmacophore rescoring of docking poses is
a decades-old idea, and we are not the first to distrust a docking score. What this repository
contributes is a pre-registered, hash-frozen, blind-holdout execution of that idea on a
clinically important target — including seven published failures and a mechanistic explanation for
them — in a field where ranked lists are routinely published with no evidence the ranking means
anything. The honest negative is the result.
