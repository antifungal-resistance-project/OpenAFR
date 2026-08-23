# Design pass: lift the T1 rigid-receptor domain ceiling

Generated 2026-08-22 (design pass, not implementation)
Branch: jajjer/what-next-v2
Status: Stage 0 RAN 2026-08-23 and FIRED ITS STOP CONDITION — see below.

> **UPDATE 2026-08-23 — Stage 0 done, and it gates the rest of this doc.** The
> validation-set audit (`scripts/domain_lifted_pool_audit.py`,
> `work/RESULTS_domain_lifted_pool.md`) ran locally. **Lifting the domain ceiling grows the
> verified resistance-active pool only from n=1 to n=2** (adds opelconazole; everything else
> stays excluded for no-C.auris-evidence or non-novelty, which a bigger pocket cannot fix).
> n=2 < the n≥4 bar, so **the binding constraint is the compound universe, not (only) the
> rigid receptor**, and Stages 1–3 below (any receptor build) are **NOT justified on the
> resistance-breaker-validation rationale**. This is exactly the cheap kill Stage 0 was
> placed first to deliver. The Options/Stages below stand as the design that *would* apply if
> the pool were ever large enough; they are not the current recommendation.

## Why we're here

The last cycle (T1) closed as a **negative**: the rigid CYP51 receptor structurally
*cannot* evaluate resistance-*breakers*, because "resistance-retention" is confounded with
"over the receptor's ~45-heavy-atom domain." Every compound that keeps activity against
resistant *C. auris* (VT-1598 = 43 heavy atoms, opelconazole = 50, oteseconazole class) sits
at or past the edge of what a frozen 5TZ1 pocket can accommodate, while the compounds the
tool docks reliably are exactly the ones Y132F defeats. See
`docs/designs/refocus-resistant-cyp51-triage.md` (T1) and `data/ligands/PROVENANCE.md`.

The user's call (2026-08-22): **don't do wet-lab outreach with no candidate — build the tool
until it can actually find resistance-breakers.** The single wall between here and that is the
domain ceiling. This doc scopes how to lift it, honestly.

## The one insight that shapes everything: there are TWO ceilings, not one

They have different causes and need different fixes. Conflating them is how we'd fool
ourselves.

1. **Domain ceiling (T1).** A *rigid* pocket physically can't seat a >45-atom ligand in a
   productive, iron-coordinating pose. This is a **receptor-geometry** problem. A flexible or
   ensemble receptor is the right lever.

2. **Decoy ceiling (preprint §4.4).** Property-matched decoys coordinate the heme iron both
   closely (2.69–2.84 Å) and on-axis (2.7°–7.5°), *indistinguishably from real azoles*, so
   six decoys outrank the best true active. This is a **discriminating-criterion** problem,
   and it is untouched by receptor flexibility. **Worse:** adding receptor degrees of freedom
   gives *decoys more ways to reach the iron too*, so a flexible protocol may make the decoy
   ceiling lower, not higher, even as it lifts the domain ceiling.

**Consequence for this design:** lifting the domain ceiling is *necessary but not sufficient*.
Any plan that opens the pocket must carry a decoy control that measures whether it wrecked the
enrichment it was supposed to preserve.

## The actual prize (why lifting the domain ceiling is worth it)

The causal chain the whole build rests on:

- The rigid domain ceiling is *why* the in-domain resistance-active pool is **n=1** (T1). The
  other resistance-retaining actives were excluded for being **over-domain**, not for being
  irrelevant.
- Lift the domain ceiling → VT-1598 / opelconazole / oteseconazole class become dockable →
  the resistance-active validation pool grows from **n=1 toward n≥4**.
- Only *then* is a pre-registered resistance-breaker gate gradeable (n=1 can't support any
  AUC/EF/permutation claim).
- A validated resistance-breaker gate → a credible novel-hit shortlist → a candidate worth a
  lab's one-shot first contact.

So the receptor work is not the goal; **rebuilding the validation set is.** The receptor is
the means. This matters for sequencing (below): we can test whether the prize even exists
before building the receptor.

## Options considered

### Option A — Vina flexible side chains (`--flex`)
Mark the pocket residues that gate large ligands (e.g. Y132, and the residues whose van der
Waals walls define the ≤45-atom cap) as flexible; Vina samples their torsions during docking.
Directly models induced fit — the physically honest picture of how a big tetrazole actually
binds.

- **Attacks:** domain ceiling, directly and physically.
- **Determinism:** preserved *if* the flex-residue set is frozen (hash-pinned, like
  `protocol.yaml`) and the seed is fixed. Vina flexible docking is still seeded/reproducible.
- **Feasibility gap:** building the rigid/flex PDBQT split needs `prepare_flexreceptor`
  (Meeko / ADFR suite). **Not in `environment.yml`** (we have obabel + vina 1.2.7 only). Adds
  a pinned dependency to prove out — and Meeko's determinism must itself be verified
  (same-input → same-sha256), or it violates the credibility contract the way a guessed
  rotamer would.
- **Risk:** most exposed to the decoy ceiling (most added DOF for decoys to exploit) and to
  search-space blowup degrading the already-fragile top-rank enrichment.

### Option B — ensemble cross-docking into deposited CYP51 structures
Dock *rigidly* into each of several real, deposited fungal CYP51 pockets (an ensemble of
induced-fit states captured by crystallography with different ligands bound), then aggregate
the geometric score across members with the **median-over-members** aggregator we already
built and tested (`validate_gate_reliability.py` — same math as the 5-seed run).

- **Attacks:** domain ceiling *only insofar as some deposited pocket happens to be open
  enough* for the big ligand. Weaker guarantee than A, but real, experimental, and honest.
- **Determinism:** cleanest fit to the project's DNA. Each ensemble member is a real structure
  extracted deterministically and frozen with its own sha256 — exactly like
  `work/receptor_A.pdb`. No new docking mode; per-member dock is the existing rigid path; the
  aggregator is already tested. Minimal new code, maximal reuse.
- **Feasibility gap:** we have only 5TZ1 in `data/structures/`. Need to pull additional fungal
  CYP51 depositions (C. albicans 5V5Z/5FSA family, other fungal CYP51s). **No experimental
  C. auris structure exists** (per prior work), so a *resistant-pocket* ensemble would mean
  grafting Y132F onto each member — compounding the graft assumption across the ensemble.
- **Risk:** may simply not open the pocket enough (the deposited pockets could all be
  small-azole-shaped). Then it buys robustness, not domain.

### Option C — hybrid (ensemble of flexible-sidechain docks)
A on top of B. Maximum induced-fit coverage, maximum cost, maximum decoy exposure. Not a
starting point; note it only as where C could go if A and B each half-work.

## Recommendation: sequence it, and let cheap steps kill the expensive ones

The expensive work (any docking) is gated on the cloud Linux/x86 VM anyway (no vina/obabel on
osx-arm64). So put the cheap, locally-runnable, potentially-fatal checks first.

**Stage 0 — validation-set audit (local, rdkit-only, no docking, no VM). DO THIS FIRST.**
Before touching the receptor, answer: *if the domain ceiling were lifted, does the
resistance-active validation pool actually reach n≥4?* Pure cheminformatics + literature: for
each candidate resistance-retaining compound, record heavy-atom count, novelty vs the azole
panel (the existing 0.70 Tanimoto filter), and quality of the *resistant-strain* activity
evidence. This mirrors "size the run before you run it" (cf. the recaller). **If the pool is
still ~n=1 even with the domain cap removed, the receptor was never the blocker — the compound
universe is — and NO flexibility work is justified.** This step can end the whole line of work
for near-zero cost, which is exactly why it goes first.

**Stage 1 — receptor approach (local build + offline logic tests). Primary: B, fallback: A.**
Start with ensemble cross-docking (B): it reuses the frozen-anchor discipline and the tested
median aggregator, so most of it is receptor-prep plumbing whose logic is offline-testable
(like the recaller orchestration was de-risked before its real run). Reach for flexible side
chains (A) only if the deposited ensemble demonstrably fails to open the pocket for the Stage-0
compounds — and if so, pin + determinism-verify the flex-prep tool first.

**Stage 2 — decoy control (VM, pre-registered). Gate before any screen.** Re-run the
property-matched decoy check under the new protocol. Pre-register the bar: does it still
separate actives from mechanism-matched decoys *at all*? If flexibility lifted the domain
ceiling but collapsed enrichment (the two-ceilings risk), we stop and say so — another honest
negative, cheaply bought.

**Stage 3 — resistance-breaker gate (VM, pre-registered).** Only if Stage 0 (pool ≥ n≥4) and
Stage 2 (enrichment survives) both clear: pre-register and run the blind-holdout +
permutation gate on the in-domain resistance-actives. This is the first run that could
produce a candidate worth an email.

## Credibility requirements (non-negotiable, carried from the existing contract)
- Any flex-residue set or ensemble membership is **frozen and hash-pinned** exactly like
  `protocol.yaml` / `work/receptor_A.pdb`; changing it is a conscious re-freeze with a
  recorded reason.
- Any new prep tool (Meeko/ADFR) must be **deterministic by construction**: same input →
  same-sha256 output, or it's refused — the same bar that made `mutate_receptor.py` refuse
  atom-adding mutations.
- The heme-iron position stays load-bearing and checked per member (every geometry distance is
  measured to it).
- Pre-register Stage 2 and Stage 3 *before* looking at results; keep any FAIL as a FAIL.

## What can be built locally now vs. gated on the VM
- **Local, now:** Stage 0 entirely; Stage 1 receptor-prep code + its offline logic tests;
  extending `mutate_receptor.py`/`prep_receptor.py` for the ensemble members; the aggregation
  wiring (reuses `validate_gate_reliability.py`).
- **VM only:** every actual `vina` dock — Stages 2 and 3, plus any redock sanity check.

## Open questions
- **Stage 0 outcome** — is the domain-lifted resistance-active pool actually ≥4? (If not, stop.)
- **Which deposited CYP51 structures** form an honest ensemble, and do any have a pocket open
  enough for a 43–50-atom ligand? (Cross-dock VT-1598 into each as the probe.)
- **Meeko/ADFR determinism** — can a flex-receptor prep be made same-input→same-sha256 in a
  pinned env? If not, Option A is off the table on credibility grounds, same as the K143R
  repacker.
- **Grafting onto ensemble members** — is a Y132F graft onto each deposited pocket defensible,
  or does compounding the graft across members multiply the assumption past honesty?
- **Pass bar for the resistance-breaker gate** — AUC/EF@5% (the durable metrics), never EF@1%
  (the preprint's discredited Bernoulli metric).
