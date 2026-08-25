# Pre-registration — metal-aware baseline (GNINA) comparison

**Written 2026-08-24, BEFORE running GNINA on any pose.** Everything below is fixed in
advance. Inherits the run-2 protocol, receptor and held-out set unchanged; adds ONE new
scorer applied to the SAME poses. Read [PREREGISTRATION_run2.md](PREREGISTRATION_run2.md)
and [RESULTS_run2_holdout.md](RESULTS_run2_holdout.md) first.

## Background — the specific objection this closes

Run 2 established, on held-out data, that the iron-approach geometric criterion (AUC 0.794)
ranks known azoles far better than **AutoDock Vina's own affinity score** (AUC 0.471,
worse than random) on the *identical* poses. That contrast is real, but it invites one
predictable rebuttal:

> Vina is *known* to mishandle metal coordination — it ignores partial charges and cannot
> model the azole-N -> heme-Fe bond. Of course it fails here. Beating vanilla Vina on a
> metalloenzyme is beating the one baseline everyone already knows is broken. Does your
> criterion beat a *competent, metal-aware* scorer?

This pre-registration answers exactly that, and only that. It does not re-open the run-2
result; it stress-tests the *headline framing* against the strongest readily-available
opponent before that framing is defended in public.

## The opponent (fixed; one, chosen in advance)

**GNINA** (CNN-rescoring docking, McNutt et al.), using its `CNNaffinity` head as the
ranking score. Chosen because it is the single most-cited learned rescorer that is (a)
documented to fix precisely Vina's discrimination weakness and (b) a drop-in re-scorer of
existing poses. It is the strongest-news opponent: beating it settles the objection;
weaker metal-aware baselines (smina/Vinardo, AutoDock4Zn) would be less informative and are
**not** part of this pre-registration.

- `CNNaffinity` (predicted binding affinity, higher = better) is the ranking score.
- `CNNscore` (pose-quality probability) is reported alongside for context, never gated.
- **Same poses, same everything else.** GNINA runs in `--score_only` mode on the exact
  run-2 held-out poses (regenerated bit-for-bit from the frozen deterministic
  `scripts/screen.sh`, same box/exhaustiveness/modes/seed/receptor). No re-docking, no
  pose regeneration under a different protocol, no per-molecule tuning. GNINA's own default
  CNN model/weights are used as shipped and recorded in the results.

## The three rankings compared (all on the same poses)

    C  (ours)     iron filter (N-Fe < 3.0 A), then best Vina affinity among survivors
                  — the run-2 PRIMARY criterion, reproduced on these poses
    G  (opponent) best CNNaffinity over ALL poses — GNINA with NO mechanistic prior
    GC (probe)    best CNNaffinity among the iron-bound (N-Fe < 3.0 A) poses only

Metrics for every mode: AUC, EF@1/5/10%, BEDROC(alpha=20), all from `rdkit.ML.Scoring`
(the same canonical implementations the gate uses). Any active with no rankable pose is
ranked LAST, never dropped — identical to every prior gate.

## Primary comparison and endpoints (fixed)

**Headline comparison: G vs C** — the practitioner's question, *"should I just use GNINA
instead of your whole pipeline?"* **Endpoints: AUC and EF@5%.** (EF@1% is deliberately
excluded: [RESULTS_reliability.md](RESULTS_reliability.md) already established the top-1%
figure is not robust under seed sampling, so it is not a defensible endpoint for anything.)

**Secondary comparison: GC vs C** — the mechanistic question, *"given the same iron-bound
poses, does a learned scorer beat plain Vina affinity?"* Same endpoints. Reported, not
headlined.

Held-out set: the run-2 set unchanged — 7 held-out azoles (`actives_holdout_final.smi`) +
348 property-matched decoys, rigid 5TZ1 chain A + heme (`work/receptor_A.pdb`).

## Pre-committed interpretation (three outcomes, fixed before seeing numbers)

Read on **G vs C** at the AUC / EF@5% endpoints. "Beats" = strictly higher point estimate;
ties within overlapping bootstrap CIs are called a TIE, not a win, either way.

- **C > G (geometry wins).** The mechanistic criterion beats a trained CNN rescorer that
  *was built to* handle metal sites — a genuinely strong result. The run-2 framing stands
  and strengthens: "geometry beats the scoring function" survives contact with a competent,
  not just a broken, scorer.

- **C ~ G (tie).** The simple, free, interpretable geometric criterion *matches* an opaque
  learned model at zero training cost and full mechanistic transparency. Reported as a
  **feature, not a loss**: the geometry recovers what the CNN finds. The run-2 headline is
  softened from "beats the scoring function" to "matches the best learned scorer and beats
  the default one, interpretably."

- **G > C (GNINA wins).** A learned metal-aware scorer beats the geometric criterion on
  this set. Then **the "geometry beats the scoring function" line is retired as the
  headline** — honestly, in the results and the preprint — and the project leads entirely
  with its actual contribution: the pre-registered, self-validating gate and the honest
  reliability accounting. The geometric criterion remains a cheap, interpretable
  *first-pass filter* (its role in mode GC), not a claimed state-of-the-art ranker. This
  outcome does not invalidate run 2 (Vina *did* rank worse than random on those poses); it
  reframes what that contrast is allowed to claim.

Whichever occurs is reported as the outcome. **This is the ONLY baseline comparison run on
these poses under this pre-registration.** No second opponent is swapped in after seeing
GNINA's numbers to manufacture a preferred result; adding smina/Vinardo/AutoDock4Zn later
requires its own separate pre-registration and is disclosed as an additional comparison.

## Limitations (fixed in advance)

1. **n = 7 actives.** Same small-sample caveat as run 2: bootstrap CIs are wide, so a
   "tie" is the honest default unless one method clears the other's CI. A single held-out
   set cannot rank two good scorers with confidence — this establishes *order of magnitude*
   ("is our criterion in GNINA's league or beaten by it"), not a precise ranking.
2. **GNINA default model, as shipped.** No model selection, no ensemble, no CNN
   fine-tuning — using a stronger GNINA configuration could only *raise* the opponent's
   score, so a C >= G result here is conservative and a G > C result is not attributable to
   a hobbled opponent.
3. **CNN scoring determinism.** GNINA CPU scoring is deterministic given fixed weights; the
   exact GNINA version and model hash are recorded in the results so the comparison is
   reproducible. GPU vs CPU numerical drift, if any, is noted.
4. **Same rigid receptor, same presumed-inactive decoys, same <= 45 heavy-atom
   applicability domain** as run 2 — every run-2 limitation carries over unchanged. This
   compares two *scorers on one benchmark*, nothing more.

## Analysis plan / reproduce

    # on a GNINA-capable (Linux/x86) host — regenerate the exact held-out poses, then rescore:
    ./scripts/run_screen2.sh                              # deterministic; reproduces work/screen2 poses
    ./scripts/run_gnina_rescore.sh work/screen2 work/gnina_scores work/receptor_A.pdb

    # grade (runs anywhere with rdkit; only parses saved GNINA output + the Vina poses):
    python scripts/rescore_gnina.py work/gnina_scores work/screen2 work/receptor_A.pdb

`scripts/rescore_gnina.py` re-checks this document's SHA-256 and the frozen protocol hash
before reporting, and refuses to headline a comparison it cannot grade cleanly.
