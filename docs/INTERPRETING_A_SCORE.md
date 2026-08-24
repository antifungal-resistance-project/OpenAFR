# Interpreting a score

You ran `scripts/score_molecule.py` (or the triage-only `scripts/applicability_triage.py`) on
your own molecule. This page explains what the output means — and, just as importantly, what it
does **not** mean. Read it once before you act on a number.

The one-line version: **a score here is a hypothesis for a wet lab, never a hit.**

---

## The two stages you see

`score_molecule.py` runs in two stages, and the first one can stop the second.

### Stage 1 — applicability triage (runs anywhere, no docking)

Before spending any docking compute, the tool asks whether the validated method is even
*claimed* to cover your molecule. Every molecule gets exactly one verdict, from a fixed
priority order:

| verdict | meaning | proceeds to docking? |
|---|---|---|
| `unparseable` | RDKit could not read the SMILES | no |
| `out-of-scope` | no aromatic nitrogen — it cannot coordinate the heme iron, so the criterion is undefined for it | no |
| `out-of-domain` | more than 45 heavy atoms — outside the pre-registered applicability domain | no |
| `near-known` | inside the envelope, but a near-duplicate (Morgan r=2 Tanimoto ≥ 0.70) of a validated active | yes |
| `in-envelope-novel` | inside the envelope **and** novel chemistry — the method's evidence covers it | yes |

Two honest facts about this gate:

- **The 45-heavy-atom domain is honest to a fault.** It is a *pre-registered* blind spot: in
  early runs, larger ligands failed reproducibly because a rigid receptor cannot model their
  induced fit. The cap is not tuned to flatter the method — it is strict enough that
  **posaconazole and itraconazole, two of the defining azole drugs, are themselves
  `out-of-domain`.** If your molecule is large, an `out-of-domain` verdict is a statement about
  the tool's validated envelope, not a judgment of your molecule.
- **`out-of-scope` is a definitional stop, not a low score.** The criterion measures a
  nitrogen-to-iron distance. A molecule with no aromatic nitrogen has no such distance, so the
  method makes no claim about it. The tool reports this and refuses to dock it, rather than
  inventing a meaningless number.

Nothing is silently dropped. A molecule the method cannot cover still gets a row — it just
carries a "not docked" note explaining why.

### Stage 2 — the N–Fe geometry read (needs the `openafr` env)

Only `near-known` and `in-envelope-novel` molecules are docked. For those, the tool docks the
molecule under the hash-frozen protocol and reports the single validated criterion:

> **the closest distance any nitrogen in your molecule reaches to the heme iron (smaller = the
> molecule coordinates the iron the way the real drugs do).**

The readout contextualizes that distance against the **validated actives' iron-bound band,
2.47–2.88 Å** — the range spanned by every azole drug that produced an iron-bound pose in the
validation run. Landing inside that band means your molecule reaches the iron like the known
drugs; it does **not** mean it is active (see below) — most molecules measured *in* that band are
compounds measured **not** to inhibit.

If docking cannot run on your machine (missing `vina`/`obabel` or an unprepared receptor), the
geometry line reads `PENDING` and the tool prints the exact commands to finish the run. It never
fakes the number.

---

## What a tight N–Fe does NOT mean

This is the part that matters most. Every one of these caveats is load-bearing and is carried
in the tool's own printed header for a reason.

1. **It is a hypothesis, not a hit.** Docking produces hypotheses; published hit rates are low.
   A tight distance is a reason to *consider* an assay, not evidence of activity.

2. **Coordinating the iron is necessary but NOT sufficient — and this is now measured, not
   assumed.** Look-alike molecules reach the iron too. That was first shown with *presumed*-inactive
   decoys (the **decoy ceiling**, `work/RESULTS_axial.md`), and it has since been confirmed against
   279 compounds **measured** not to inhibit *C. albicans*
   (`work/RESULTS_verified_inactives.md`): **seven of them coordinate the iron more closely than
   the best real drug**, and of the 39 molecules landing inside the 2.47–2.88 Å band, **33 are
   compounds measured not to work and only 2 are true actives.** A single molecule landing in the
   band therefore tells you very little on its own. The method's evidence is *statistical
   enrichment across a ranked library*, never a verdict on one molecule.

   **The sharpest form of this caveat:** on that measured benchmark the criterion scores AUC 0.810
   against non-azole chemotypes but only **0.650 against azole analogues that were measured not to
   work** — below its own 0.70 pass bar. Most of what it detects is *"is this azole-like?"*, not
   *"does this azole work?"* If your molecule is itself an azole, the criterion is close to
   uninformative about whether it will be active.

3. **The durable metric is a shortlist, not a razor-top rank.** The validated performance is
   **AUC ≈ 0.79** — a trustworthy top-~5% shortlist. An earlier eye-catching "12.68× enrichment
   at the top 1%" figure was a **single-run artifact**: it collapsed to 0.00× once each ligand
   was docked over 5 seeds (`work/RESULTS_reliability.md`). Do not read a single molecule's rank
   as a razor-sharp prediction.

4. **It is meaningful only if the gate has passed.** The whole criterion is only trustworthy on
   a protocol where `scripts/validate_gate2.py` has PASSED under the frozen `protocol.yaml`. If
   you change the receptor, box, or protocol, re-run the gate before trusting any geometry read.

5. **It says nothing about resistance-breaking.** The method was validated on wild-type CYP51.
   A separate pre-registered test found the criterion **failed** its gate in the *C. auris*
   Y132F resistance pocket (`work/RESULTS_auris_Y132F.md`). Do not read an in-band score as
   evidence a molecule defeats azole resistance.

---

## So what is a score good for?

Used honestly, this tool does one useful thing well: **it triages a library so you don't waste
an assay — or a claim — on molecules the validated method has no basis for.** An
`in-envelope-novel` molecule that lands in the actives' iron-bound band is a defensible
*candidate to prioritize for the next, more expensive step* — not a discovery. That prioritization
is strongest when the molecule is a **new chemotype**; for another azole analogue, the measured
benchmark says the criterion is close to uninformative (caveat 2).

If you have a compound you want to check, that is exactly the intended use. Bring the SMILES,
run the front door, and read the verdict with this page open.

---

## Where these claims come from

- Validation run and pre-registration: `work/RESULTS_run2_holdout.md`, `work/PREREGISTRATION_run2.md`
- Reliability / the top-1% artifact: `work/RESULTS_reliability.md`
- Decoy ceiling: `work/RESULTS_axial.md`
- Resistance-pocket negative: `work/RESULTS_auris_Y132F.md`
- The full method write-up: `work/PREPRINT_geometry_ceiling.md` (+ `.pdf`)
