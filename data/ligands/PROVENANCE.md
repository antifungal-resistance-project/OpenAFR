# Ligand set provenance + active-set expansion spec (T1)

This file documents how the active sets were built and specifies the disciplined
expansion called for by the 2026-08-21 eng review + outside voice
(`docs/designs/refocus-resistant-cyp51-triage.md`, task T1).

## Current sets (as of 2026-08-21)

### Training actives — `actives.smi` (n=8)
Classic azole antifungals, used to define the geometric criterion:
fluconazole, voriconazole, itraconazole, posaconazole, ketoconazole, clotrimazole,
miconazole, isavuconazole.

### Held-out actives — `actives_holdout_final.smi` (n=7)
Blind test set, fixed in `work/PREREGISTRATION_run2.md` and reused unchanged by the
Y132F transfer test (`work/PREREGISTRATION_auris_Y132F.md`). Ten azoles were retrieved;
three (econazole 0.79, tioconazole 0.76, terconazole) were removed as near-duplicates of
training actives (Morgan r=2 Tanimoto ≥ 0.70). Remaining 7:
efinaconazole, luliconazole, sertaconazole, oxiconazole, bifonazole, fenticonazole,
butoconazole. (`actives_holdout.smi` is the pre-filter list of 10.)

### Decoys — `decoys.smi`, `decoys_holdout.smi`
Property-matched presumed-inactives (see `scripts/make_decoys.py`).

### Verified inactives — `verified_inactives.smi` (n=279, added 2026-08-23)
Compounds **measured** not to inhibit *C. albicans*, not presumed inactive — the independent
decoy construction the preprint's stopping rule requires (§5) and the exit §4.1 names for the
decoy-construction bind. Built by `scripts/make_verified_inactives.py` from an 88,695-record
ChEMBL bioactivity pull (87,879 whole-cell *C. albicans* + 816 on the CYP51 enzyme target
CHEMBL1780); the evidence rule lives in `openafr/inactives.py`.

A compound qualifies only under **consistent inactivity**: at least one inactive-grade record
(MIC `>` X ≥ 32 µg/mL, MIC `=` X ≥ 64 µg/mL, or an explicit "Not Active" comment) and **zero**
active or ambiguous records anywhere, including on the enzyme target. Of 43,004 compounds ever
tested against *C. albicans*, 8,908 qualify; 279 survive the filters.

Those filters are byte-identical to `make_decoys.py` — property matching (MW ± 25, cLogP ± 1.5,
rotatable bonds ± 2, HBD ± 1, HBA ± 2), ≥ 1 aromatic nitrogen, Tanimoto < 0.35 against every
active, ≤ 50 per active — plus the ≤ 45 heavy-atom applicability domain. The swap is
single-variable: only the provenance of "inactive" changes. Achieved match vs the 7 held-out
actives: MW 379.9 / 392.5, cLogP 4.8 / 5.7, heavy atoms 26.3 / 25.4, rotatable bonds 4.6 / 5.4.

**Read this before comparing it to `decoys_holdout.smi`:** 58.4 % of the verified inactives
carry an azole ring, against 27.4 % of the presumed decoys, because the only compounds anyone
measures against *C. albicans* are antifungal medicinal-chemistry programs. More than half of
this set is **failed azole analogues** — same warhead as the drugs, measured not to work. That
makes it a harder benchmark, not a like-for-like one. Frozen with
`work/PREREGISTRATION_verified_inactives.md`; sha256 in `work/PREREG_verified_inactives.sha256`.

### Continuous-potency sets — `potency_enzyme.{tsv,smi}`, `potency_wholecell.{tsv,smi}` (added 2026-08-29)

For the continuous-potency validation (issue #81, `work/PREREGISTRATION_potency.md`): do the
N–Fe geometry scores track *how potent* a molecule is, not just active/inactive? Built by
`scripts/make_potency_set.py` from a ChEMBL pull of every bioactivity record **carrying a
`pchembl_value`** (−log10 molar potency) for two targets, each compound collapsed to one median
pChEMBL by `openafr.potency.median_pchembl` (exact-relation IC50/Ki/Kd/EC50/Potency only; MIC
carries no pchembl and never enters). Only compounds inside the validated applicability envelope
(≥1 aromatic N, ≤45 heavy atoms — same triage as the front door) are written to `.smi` for
docking, because the method makes no claim outside it.

- `potency_enzyme` — target **CHEMBL1780** (*C. albicans* CYP51 enzyme), the direct on-target
  binding potency: **32 in-domain** compounds (28 novel, 4 near-known azoles), pChEMBL 5.94–8.00.
  The PRIMARY test — no whole-cell confound, but a narrow ~2-log range at small n.
- `potency_wholecell` — target_organism *Candida albicans*, whole-cell growth-inhibition
  potency: **441 in-domain** compounds, pChEMBL 4.00–10.97. The SECONDARY — well-powered but
  confounded by uptake/efflux/target abundance, so context rather than proof.

The `.tsv` carries every compound's median pChEMBL, record count, standard_types, heavy-atom
count, triage verdict and nearest active; the `.smi` is the in-domain dockable subset. Frozen
with `work/PREREGISTRATION_potency.md`; sha256 in `work/PREREG_potency.sha256`.

### External validation holdout — `external_actives.smi` (n=50), `external_decoys.smi` (n=335) (added 2026-08-30)

For the fresh external-holdout validation (issue #82, `work/PREREGISTRATION_external.md`):
does the frozen geometric criterion generalize to an active/decoy set **no gate has ever
scored**, closing the overfitting objection? Built by `scripts/make_external_holdout.py`,
same provenance discipline as every set — only the provenance of "held out" changes:

- `external_actives` — a deterministic seed-**82** random sample of the measured-active
  *C. albicans* azoles in `verified_actives.smi` that are **not** in `verified_actives_sample.smi`.
  The ensemble-confirm look (#9) docked only the seed-42 n=200 sample, so these ~3,400 actives
  were never docked or scored; the sample of 50 is genuinely never-touched. Domain re-asserted
  (≤ 45 heavy, ≥ 1 aromatic N).
- `external_decoys` — property-matched presumed-inactives built by the frozen rule
  (`openafr.inactives.select_inactives`: MW ± 25, cLogP ± 1.5, rotB ± 2, HBD ± 1, HBA ± 2;
  ≥ 1 aromatic N; Morgan r=2 Tanimoto < 0.35 vs every external active; ≤ 45 heavy), fetched
  fresh from ChEMBL under seed 82 and **disjoint by ChEMBL id** from every set the project has
  used (all actives, both decoy sets, verified inactives, both potency sets — 4,903 ids
  excluded), so no decoy can secretly be a known active or a reused decoy.

Result (`work/RESULTS_external.md`): mode-C **AUC 0.715** (permutation p < 0.0001, bootstrap
95% CI 0.630–0.795), EF@5% 4.62× — SUPPORTED-BUT-WIDE per the pre-registered rule (point
estimate clears the 0.70 bar; CI dips below at n=50). Frozen with
`work/PREREGISTRATION_external.md`; sha256s in `work/PREREG_external.sha256`.

### Temporal (prospective) holdout — `temporal_actives.smi` / `temporal_decoys.smi` (scaffold, release-gated, added 2026-08-31)

The **stronger follow-up** the external holdout's own comment named (issue #82,
`work/PREREGISTRATION_temporal.md`). #108 above answered the overfitting objection with a
never-*scored* subset of the *same* curated population; its one honest residual was "a clean
never-touched holdout, but **not a later time slice**." This closes that: the split is
**temporal**, not just never-scored — every molecule's earliest antifungal / CYP51 activity is
reported in a ChEMBL document dated strictly after the tuning corpus (document year > `2026`),
so it could not have informed the criterion or its 3.0 Å cutoff even in principle. Built by
`scripts/make_temporal_holdout.py` (`is_temporal` predicate), same frozen decoy discipline as
every set, disjoint by ChEMBL id from every used set **including the #108 external sets**.

**Release-gated (scaffold only today).** The set materializes only once a ChEMBL release newer
than the 2026-08 tuning release accrues post-cutoff (> 2026) antifungal / CYP51 documents — the
"needs a newer ChEMBL release" wall #82's #108 comment named. The current CYP51/azole universe
is exhausted, so the builder writes nothing (NOT-YET-AVAILABLE) and `validate_gate_temporal.py`
reports `not-yet-run (release-gated)`, mirroring the metal-aware baseline harness (#83). The
rule, cutoff, disjointness, domain and three-outcome interpretation (reused verbatim from #108)
are all **frozen now** (`work/PREREG_temporal.sha256`) and unit-tested offline
(`tests/test_temporal_gate.py`); only the fetch + dock + grade wait on the newer release.

## The problem T1 fixes

Two independent findings converged on the active set:

1. **Every active is a defeated azole.** Re-discovering drugs that resistant *C. auris*
   already beats does not demonstrate the tool can find resistance-*breakers* — the actual
   goal of any wet-lab shortlist. Mutant AUC 0.805 ≈ WT 0.794: the criterion is nearly
   *insensitive* to the resistance mutation, so "broad enrichment transfers to Y132F" is
   close to a tautology.
2. **n=7 is too thin.** EF@5% on 7 actives / 350 decoys counts ~18 molecules; one active
   moves it ~1.4x. AUC and the permutation test have large variance at 7 positives.

## Expansion design (differential-activity, not "more azoles")

The discriminating test: add CYP51 inhibitors that **retain measured activity against
azole-resistant *C. auris*** (the newer tetrazoles), so a genuinely resistance-aware
criterion should rank the still-active compounds **above** the defeated azoles *in the
mutant pocket*. That converts the validation from "rank azoles vs decoys" (tautological on
the mutant) to "rank still-active vs defeated, in the resistant pocket" (a real claim).

**Candidate additions — SMILES VERIFIED from PubChem 2026-08-22 (isomeric/stereo form,
copied from the authoritative record, not AI-generated):**

> **CORRECTION 2026-08-22 (identity error in the 2026-08-21 curation — caught before any
> pre-registration froze).** The 2026-08-21 pass recorded "VT-1598 (quilseconazole),
> CID 91886002" — but CID 91886002 is **quilseconazole = VT-1129** (a *Cryptococcus*-focused
> tetrazole, `…OC(F)(F)F` tail, C22H14F7N5O2, 36 heavy atoms), **not** VT-1598. VT-1598 and
> VT-1129 are different Viamet tetrazoles. The C. auris evidence cited below (PMID 30530603)
> is genuinely VT-1598, so the *drug choice* was right; the *structure recorded* was wrong —
> exactly the "a wrong SMILES silently poisons the run" failure discipline constraint #3
> guards against. Corrected identities:

- **VT-1598** — PubChem **CID 126715974**, C31H20F4N6O2, MW **584.5**, **43 heavy atoms**.
  `C1=CC(=CC=C1COC2=CC=C(C=C2)C#CC3=CN=C(C=C3)C([C@](CN4C=NN=N4)(C5=C(C=C(C=C5)F)F)O)(F)F)C#N`
  IUPAC: 4-[[4-[2-[6-[(2R)-2-(2,4-difluorophenyl)-1,1-difluoro-2-hydroxy-3-(tetrazol-1-yl)propyl]-3-pyridinyl]ethynyl]phenoxy]methyl]benzonitrile.
  A large *extended* tetrazole (difluorophenyl-difluoro-tetrazolylpropanol-pyridine core +
  ethynyl-phenoxy-benzonitrile arm) — **not** a close analogue of quilseconazole.
- **Quilseconazole (VT-1129)** — PubChem **CID 91886002**, C22H14F7N5O2, MW 513.4, 36 heavy
  atoms. `C1=CC(=CC=C1C2=CN=C(C=C2)C([C@](CN3C=NN=N3)(C4=C(C=C(C=C4)F)F)O)(F)F)OC(F)(F)F`.
  On record only; **NOT added** — its published activity is vs *Cryptococcus*, no
  C. auris-specific MIC surfaced (same rule that drops oteseconazole below).
- **Oteseconazole (VT-1161)** — PubChem **CID 77050711**, C23H16F7N5O2, MW 527.4.
  `C1=CC(=CC=C1C2=CN=C(C=C2)C([C@](CN3C=NN=N3)(C4=C(C=C(C=C4)F)F)O)(F)F)OCC(F)(F)F`.
  On record only; **NOT added** (no *C. auris*-specific MIC — see below).

All three carry a **tetrazole** metal-binding group (`…N=N…`) rather than the
imidazole/triazole of every current active — genuinely distinct chemistry that clears the
Morgan-Tanimoto novelty filter.

**Applicability-domain flag (VT-1598).** At **43 heavy atoms / MW 584** VT-1598 is the
*largest* molecule in the entire active set (every existing active is ≤ 36 heavy atoms) and
sits right at the run-2 declared domain edge (valid ≤ 45 heavy atoms; posaconazole/itraconazole
class ≥ 52 failed in run 1 from induced fit a rigid receptor cannot model). A rigid-receptor
false negative (VT-1598 fails to place a coordinating pose because of induced fit, not
inactivity) would rank it LAST and sabotage the very test it is added for. This risk MUST be
pre-declared in the run-3 pre-registration.

**Widening pass 2026-08-22 (before concluding n=1) — two more named candidates evaluated:**

- **Ravuconazole** (CID 467825, C22H17F2N5OS, MW 437.5, 31 heavy; active vs *C. auris*,
  PMID 33609301) — **REJECTED on novelty:** Morgan r=2 Tanimoto **0.879** vs isavuconazole
  (a training active). It is essentially isavuconazole; fails the < 0.70 filter.
- **Opelconazole / PC945** (CID 121383526, C38H37F3N6O3, MW 682.7, **50 heavy atoms**; potent
  vs a global *C. auris* collection, PMID 31325309) — novel (max Tanimoto 0.553, nearest
  posaconazole) but **over the applicability domain:** 50 heavy atoms is posaconazole-class,
  and posaconazole/itraconazole (≥ 52) failed reproducibly in run 1 from induced fit a rigid
  receptor cannot model. Docking it would almost certainly reproduce that false negative.
- The recent (2024–25) resistance-active CYP51 hits are unnamed research compounds
  ("compound 7", "V23", "compound 22") with no PubChem registration → fail discipline
  constraint #3 (verified SMILES) → not includable.

### The structural confound this widening pass exposes (the real T1 finding)

Every compound that **retains activity against azole-resistant *C. auris*** is a large,
long-tailed molecule — VT-1598 (43 heavy), opelconazole (50), the oteseconazole/posaconazole
class — precisely the size range the rigid-receptor geometry tool is **NOT validated for**
(domain ≤ 45; posaconazole-class ≥ 52 failed in run 1). The small compounds the tool docks
reliably (fluconazole, voriconazole ≤ ~33 heavy) are exactly the ones Y132F defeats.
**Resistance-retention is confounded with over-domain molecular size.** This is not a
data-gathering shortfall that a harder search fixes; it is a property of the chemistry
(selectivity/potency against the mutant pocket is bought with extended side chains) meeting a
property of the tool (rigid receptor, small-ligand domain). It answers the refocus plan's
unresolved decision #2: **a resistance-*breaker* validation is not achievable with this rigid
tool** — the molecules that would prove it are the ones it cannot honestly dock.

**Activity evidence vs *C. auris* (literature search 2026-08-21 — verify numbers against the
primary papers before freezing):**

- **VT-1598 — SUPPORTED. Use this one.** Wiederhold et al., "The Fungal Cyp51-Specific
  Inhibitor VT-1598 Demonstrates In Vitro and In Vivo Activity against *Candida auris*,"
  *Antimicrob Agents Chemother* 2018, PMID 30530603 (doi 10.1128/aac.02233-18). Reported
  modal MIC **0.25 µg/mL**, range **0.03–8 µg/mL** across 100 clinical *C. auris* isolates,
  activity against azole-resistant isolates + in vivo burden reduction. This is the
  differential-activity evidence the design needs. (Numbers are from the abstract/summary;
  confirm exact values in the paper before the pre-registration freezes.)
- **Oteseconazole (VT-1161) — NOT SUPPORTED for this use; DROP.** Its strong MIC data is vs
  *C. albicans*/*C. glabrata* (incl. fluconazole-resistant) from vulvovaginal-candidiasis
  studies; no *C. auris*-specific MIC surfaced. Including it as a "resistance-retained" active
  would be the exact unsupported assumption this design guards against. Structure stays on
  record (CID 77050711) but it is NOT added as an active without a real *C. auris* source.

**Net:** add **VT-1598 only** (CID 126715974) as the resistance-retained active. This is the
*only* verifiable resistance-active CYP51 compound the literature pass surfaced — the pool is
n=1, not "a set of newer tetrazoles." (The earlier "0.823 self-similarity" caveat was an
artifact of the identity error: it was oteseconazole vs *quilseconazole* — the two smaller
tetrazoles differing only at the –OCH2CF3 / –OCF3 tail — not oteseconazole vs the real
VT-1598.)

**Consequence for the design (resolves the plan's unresolved decision #2).** The refocus plan
asked whether a credible resistance-*breaker* validation is even possible without a set of
known-active-vs-resistant compounds. The honest answer, after the widening pass: **no — and the
reason is structural, not a data shortfall** (see the confound above). With n=1 verifiable
in-domain resistance-active (VT-1598, itself at the domain edge) there is no
differential-activity permutation test, and the compounds that would populate one are
systematically over-domain. This is the more valuable result than any single-probe run would
have been: it says *what the tool can and cannot honestly claim*, which is exactly what the
T3 wet-lab package must lead with.

**Decision (2026-08-22):** do NOT spend a cloud run on the VT-1598 single-probe
(`work/PREREGISTRATION_run3.md` stays a DRAFT, frozen-ready but de-prioritized — a domain-edge
n=1 anecdote is low-yield). Carry the confound finding into T3 instead: OpenAFR credibly
triages **broad azole-like enrichment** on both the WT and Y132F pockets (AUC ~0.79–0.81), and
the honest pitch stops there — it must **not** be sold as finding resistance-breakers, both
because n=1 and because resistance-retention is confounded with over-domain size. Revisit run 3
only if a flexible-receptor / induced-fit protocol is ever added (would lift the domain ceiling
that makes the probe uninformative today).

## Discipline constraints (do NOT shortcut)

1. **New pre-registration, not an edit of the frozen holdout.** The existing holdout is
   sha-pinned and already graded. Any expanded set defines a NEW pre-registered run
   (`work/PREREGISTRATION_run3.md` or similar) with its own frozen list + hash. Never add
   compounds to `actives_holdout_final.smi` and re-grade an existing run — that is the
   p-hacking move D5 already rejected.
2. **Novelty filter vs BOTH sets.** Each new active must pass the same Morgan r=2 Tanimoto
   < 0.70 filter against training AND the existing holdout, so it is genuinely new chemistry.
3. **Verified inputs only.** Each new compound needs: a PubChem CID, the canonical isomeric
   SMILES copied from that authoritative record (NOT AI-generated), and a primary-literature
   citation for its activity against a resistant strain. Record CID + citation inline here.
   A wrong SMILES silently poisons the run — treat this like the receptor hash gate.
4. **Docking is gated on the cloud Linux/x86 VM** (vina/obabel absent on osx-arm64), same as
   every other screen. Curation + pre-registration here are local; the run is not.

## Status

- [x] Current sets documented (this file).
- [x] **CID/SMILES identity corrected 2026-08-22** (see CORRECTION box above): VT-1598 =
      **CID 126715974** (C31H20F4N6O2, MW 584.5, 43 heavy atoms), NOT CID 91886002
      (= quilseconazole/VT-1129). Oteseconazole = CID 77050711. All re-verified from PubChem.
- [x] Literature check for MIC vs resistant *C. auris*: VT-1598 SUPPORTED (PMID 30530603,
      modal MIC 0.25 µg/mL); oteseconazole + quilseconazole NOT supported for *C. auris* →
      dropped. The verifiable resistance-active CYP51 pool is **n=1 (VT-1598)**.
- [x] Novelty-filtered on the CORRECTED VT-1598 vs the 15 existing actives (Morgan r=2,
      2048-bit, `conda run -n openafr`, rdkit 2025.09.5, 2026-08-22): VT-1598 max Tanimoto
      **0.312** (nearest fluconazole) → PASS (< 0.70). (The earlier "0.358 / 0.823 caveat"
      was computed on quilseconazole, not VT-1598; discard it.)
- [x] **Widening pass done (2026-08-22):** ravuconazole rejected (novelty 0.879), opelconazole
      over-domain (50 heavy). Conclusion: verifiable in-domain resistance-active pool = n=1,
      and resistance-retention is confounded with over-domain size (see finding above).
- [x] **Decision:** run 3 de-prioritized — no cloud run for an n=1 domain-edge probe. Draft
      pre-registration kept frozen-ready (`work/PREREGISTRATION_run3.md`, DRAFT). The T1
      deliverable is the confound finding, which feeds T3.
- [ ] **NEXT (T3):** fold the confound finding into the wet-lab package — pitch broad
      azole-enrichment triage (AUC ~0.79–0.81, WT + Y132F), explicitly NOT resistance-breaker
      discovery. See `docs/designs/refocus-resistant-cyp51-triage.md` T3.
- [ ] (Deferred) Dock the VT-1598 probe only if a flexible-receptor protocol lifts the domain
      ceiling; grade on AUC/EF@5% (never EF@1%).
