# Pre-registration — FKS1 novel-variant POSITION CONTEXT (moonshot #4, Rung 1)

> **STATUS: FROZEN 2026-10-12, before the derivation was wired into the live weather
> page.** This document fixes exactly what position context claims, the four things it
> refuses to claim, and how it is verified. It crosses the FKS1 track's declared
> "detection-only, no structural verdict" boundary, so under this project's
> proof-before-claim discipline it must declare its scope and non-claims in advance
> (`docs/DIAGNOSTIC_GO_NO_GO.md`). The repro gate (`scripts/repro.py`) re-checks this
> document's SHA-256 via `work/PREREG_fks1_position_context.sha256` and fails on any drift.

## Why this exists — surveillance is blind on first sighting

The certified FKS1 caller (`openafr/fks1_caller.py`) tags only the literature-pinned
resistance panel (HS1 F635C/Y, S639F/P/Y, D642Y; HS2 R1354S; HS3 M690I, W691L). A
substitution it has never seen inside a hot-spot flows through the verdict
(`openafr/verdict.py`) and the detection-only alert (`openafr/alert.py`) as a bare
`UNCHARACTERIZED_VARIANT` with `structural = None` and **no context at all** — a reader
cannot tell whether the novel change lands on the S639 anchor (where F/P/Y are canonical
resistance markers) or several residues away in benign background. That is exactly the
variant genomic surveillance exists to catch early, surfaced with the least possible signal.

## What position context IS — a deterministic derivation, not a predictor

`openafr/fks1_position.position_context(token)` is a pure function over the **already-pinned
caller panel and reference numbering** (`FKS1_WINDOWS` / `RESISTANCE_PANEL` / `LOW_PPV_PANEL`,
reference `XM_085597048.1`, CDS SHA-256 `d50ececa…`). For an uncharacterised hot-spot token it
reports geography only:

- `window` — which hot-spot (HS1/HS2/HS3) contains the residue;
- `anchor` + `anchor_distance` — the window anchor (e.g. S639) and signed residue offset;
- `at_panel_residue` + `panel_letters_at_position` — is this residue one the panel already
  knows is mutable, and which letters;
- `nearest_core_marker` + `distance_to_nearest` — closest **core** (high-PPV) marker position
  (the low-PPV positions D642/M690 are excluded, so distance measures distance to a marker we
  would actually call);
- `tier` — one of, strongest first: `known_low_ppv_marker` (the token itself is a known-but-
  low-PPV letter the verdict abstains on) > `novel_substitution_at_panel_residue` (a novel
  letter at a panel residue, e.g. S639L, W691C) > `in_hotspot_adjacent_to_panel` (within 2 aa
  of a core marker) > `in_hotspot_background`.

Because it reuses the caller's panel, the two can never disagree; a token the caller tags as a
**core** marker is already `RESISTANCE_MARKER_DETECTED` and is given **no** context (returns
`None`), as is any token outside every hot-spot window or any unparseable string.

## The four things it DOES NOT claim (the frozen non-claims)

1. **No drug-FIT change.** It never asserts the novel residue weakens or strengthens
   echinocandin binding. The azole pocket engine (`openafr/structural.py`) rests on
   fluconazole↔heme-iron geometry in CYP51; echinocandins coordinate **no metal** and FKS1 is
   a glucan synthase, so that geometry does **not** transfer (a real structural verdict is
   Rung 2: a glucan-synthase structure + its own prereg). 
2. **No phenotype.** It never predicts resistant/susceptible for the novel letter.
3. **No probability or score.** There is no numeric likelihood field, consistent with the
   concordance-only verdict contract (`docs/DIAGNOSTIC_VERDICT_CONTRACT.md`).
4. **`structural` stays null for FKS1.** Context is a *separate* field (`position_context` /
   per-token `position`); it never populates `structural`, so the detection-only boundary and
   its standing caveat (`alert.FKS1_DETECTION_ONLY_CAVEAT`) remain intact. An FKS1 alert still
   never reaches ACT-NOW (that tier requires a weaker-fit verdict this does not supply).

## Verification (fixed now) — correctness, not a phenotype bar

Nothing predictive is claimed, so there is **no sensitivity/specificity/AUC to measure and no
GCP compute to spend** — the honest verification of a deterministic derivation is that it
computes the pinned geography correctly and preserves the non-claims. That is locked by unit
tests (`tests/test_fks1_position.py`, plus the verdict/alert integration in
`tests/test_verdict.py` and `tests/test_fks1_alert.py`), which assert, among others:

- a core panel marker (S639F, W691L) → `None`;
- a low-PPV letter (D642Y, M690I) → `known_low_ppv_marker`, never "novel";
- a novel letter at a panel residue (S639L; the anti-circularity W691C) →
  `novel_substitution_at_panel_residue` at the correct anchor/window;
- adjacent-vs-background split by residue distance only;
- outside-window / unparseable → `None`;
- the derivation stays consistent with `fks1_caller.is_core_panel_token` (no drift);
- the FKS1 verdict carries `position_context` while `structural` is **still** `None`, and the
  ERG11 object never gains the key.

## Interpretation, committed now

- **This is an escalation of surveillance signal, not a resistance claim.** A
  `novel_substitution_at_panel_residue` alert means "a never-seen change sits on a residue
  where other substitutions are known resistance markers — worth characterising", and is
  reported with the detection-only caveat still attached. It raises *attention*, never a
  verdict.
- **It does not change any count or threshold.** Prevalence and emergence are computed exactly
  as before (`is_panel_token` unchanged); only the *presentation* of the uncharacterised tail
  gains context.

## What is NOT in scope (declared, not silently dropped)

- **A structural fit/proximity verdict (Rung 2).** Mapping the novel residue onto a
  glucan-synthase (cryo-EM / AlphaFold) structure and scoring echinocandin-contact proximity is
  a separate, later increment with its own frozen prereg, anti-circularity truth set, and
  (likely) a priced run. This Rung 1 deliberately stops at sequence/position.
- **Phenotype or MIC prediction for novel letters.** Out of scope, as above.
- **ERG11.** The ERG11/azole path already carries a true structural verdict where mappable
  (#135); position context is an FKS1-only field and the ERG11 object is unchanged.
