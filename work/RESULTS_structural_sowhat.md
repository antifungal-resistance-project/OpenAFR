# Structural so-what — the azole-fit consequence of a flagged mutation (issue #24)

**What:** the actionable tip of the C. auris azole early-warning loop. The emergence
detector (#22) says a mutation is *new* or *rising*; the mapper (#23) says *where* it sits
in the modeled azole pocket and whether we can build it. Neither answers the question a
clinician actually asks: **if this mutation spreads, does fluconazole still fit the
pocket — and how much worse?** This stage answers it, or honestly declines to.

**Code:** `openafr/structural.py` (library) · `scripts/structural_sowhat.py` (CLI:
`estimate <tokens>`, `from-emergence`) · `tests/test_structural.py` (14 offline tests).
Reproduce: `python scripts/structural_sowhat.py estimate Y132F K143R F126L V125A TR34`.

It builds directly on the #23 mapping verdict and reuses the same modeled receptor the
whole project docks against — the 5TZ1 chain-A scaffold (`work/receptor_A.pdb`), same heme,
same iron, same frame as the frozen protocol and the #13 *C. auris* Y132F port.

## The core problem this stage is careful about

A wrong **"the drug still works"** is the single most dangerous output the whole track
could emit — it would tell a clinician to keep prescribing an azole that is quietly
failing. So the estimator is deliberately narrow about what it will claim, and sorts every
flagged mutation into one of three **evidence tiers**, each with an explicit confidence.

## The three evidence tiers

| tier | confidence | when | what it claims |
|---|---|---|---|
| **empirical-dock** | `measured` | a real, pre-registered #13 dock exists | the graded gate outcome + numbers |
| **geometric-estimate** | `estimate` | deletion-buildable, not yet docked | a *direction* only, never a number |
| **no-call** | `none` | unbuildable, or not a residue at all | "no confident call" (plausible / low-prior) |

## The verdicts on the canonical C. auris ERG11 mutations

```
Y132F  [MEASURED  weaker-fit]              gate FAIL; docked in #13
V125A  [ESTIMATE  minimal-direct-effect]   buildable, removed atoms are second-shell
F126L  [NO-CALL   plausible]               pocket-contact but add-atom -> needs repacker
K143R  [NO-CALL   low prior]               second-shell AND add-atom -> unbuildable
TR34   [NO-CALL   —]                        not a point substitution -> no residue to place
```

### Tier 1 — Y132F: measured, reusing the #13 dock (not re-run)

For the one mutation we have actually run through the port pipeline, the estimator reports
the **measured, pre-registered** outcome from
[RESULTS_auris_Y132F.md](RESULTS_auris_Y132F.md) — it does *not* re-run a multi-hour docking
exercise. The recorded result:

> Gate **FAIL**. Broad iron-reach separation **survives** (AUC 0.805 vs wild-type 0.794) —
> azoles as a class still coordinate the catalytic iron — but the **top-rank fit advantage
> collapses**: property-matched decoys tie and edge past the best real drugs at the extreme
> tip, dropping EF@1% from **12.68× to 0.00×**.

**Plain verdict:** fluconazole is **not structurally *excluded*** — azoles still reach the
iron — but the best-fitting azoles lose their top-rank edge. That is a real, if subtle,
degradation, and it ships with every caveat the #13 result declared, always including the
one that matters most here: **#13 docked a held-out azole *class*, not fluconazole
specifically**, so the fluconazole-level statement is an inference, not a measurement.
The collapse also sits at the tool's sub-0.1 Å per-pose noise floor, and a rigid single-atom
deletion under-models clinical Y132F.

### Tier 2 — V125A: a deterministic direction, never a number

V125A is deletion-buildable, so the estimator computes exactly which atoms the deletion
removes (the two γ-methyls CG1/CG2) and measures their distance to the co-crystal azole:
both are **≥ 8 Å** — well outside the 4.5 Å contact shell. No drug contact is lost, so the
predicted **direct** fit change is **minimal**; any real effect would be indirect (pocket
dynamics) and a static model cannot call it. The verdict prints the `mutate_receptor.py`
command to dock it if a number is wanted — the estimate never fabricates one.

### Tier 3 — F126L / K143R / TR34: honest no-calls

- **F126L** lines the drug-contact surface (3.68 Å), so a fit effect is *plausible* — but
  F→L **adds** side-chain atoms and the pinned env has no repacker, so the mutant cannot be
  built or docked. Reported as **plausible but unresolved**, never as a guessed side chain.
- **K143R** is **second-shell** (7.39 Å, no contact) **and** add-atom unbuildable → **low
  prior** for a direct effect and no dockable path. Reported as unresolved, *not* "no effect".
- **TR34** is a promoter / tandem-repeat token, not a point substitution → no residue to
  place in the pocket → **no structural fit call at all**.

## The independent corroboration

The Tier-1 (measured) and Tier-2 (geometric) machinery agree where they overlap. Y132F's
deterministic geometry — the removed para-hydroxyl **OH sits 4.01 Å from the co-crystal
azole**, inside the contact shell — predicts **weaker fit**, the *same direction* as the
measured dock. The lost atom is the para-OH, exactly the hydrogen-bond donor the literature
names as Y132F's azole-resistance mechanism. A test cross-checks that the geometric
direction and the registry direction stay consistent, so a later edit that flipped one
without the other would fail.

## The honesty bar (same as the main repo)

- **Direction without a dock is labeled an estimate, never a measurement** — Tier 2 carries
  no ΔΔG and prints the command to get a real one.
- **"The drug still fits" is the claim we are most reluctant to make** — even the measured
  Y132F tier says broad iron-reach *survives* while the top-rank *advantage* collapses, and
  flags the class-vs-fluconazole gap every time.
- **"No confident call" is a first-class, expected outcome** — 3 of the 5 canonical tokens
  land there, and each says *why* (unbuildable / not a residue) rather than defaulting to a
  reassuring "no effect."

## The loop closes

`structural_sowhat.py from-emergence` runs the #22 demo detector and estimates its own
output end-to-end: the four flagged tokens come back as one measured verdict (Y132F), and
three honest no-calls — a complete **detect → map → estimate** loop with no manual hand-off.
`openafr.structural.estimate_signals` accepts emergence signal dicts or bare tokens alike,
so #25 (alert composition) can consume one object per flagged mutation.

## What is deferred, honestly

- **Add-atom mutations (F126L, K143R) get no fit call.** Their mutant side chains need a
  pinned repacker — the same documented follow-up as the auris port (Path A, step 2). The
  pocket verdict is real; the fit verdict waits on repacker tooling.
- **Only Y132F carries a measured number.** Every other buildable mutation is a
  direction-only estimate until its dock is run and added to `EMPIRICAL_DOCKS` with its
  pre-registration provenance.
- **Fluconazole-specific dock is not yet run.** The strongest evidence to date is the
  held-out azole *class* result; a fluconazole-specific dock against the Y132F pocket would
  turn the Tier-1 inference into a direct measurement.
