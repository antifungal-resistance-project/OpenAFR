# T5 — the go/no-go validation run

Date: 2026-07-30
Test set: 8 known azoles hidden among 396 property-matched decoys (1:50), one identical
protocol for every molecule (box 30 A @ 67.4/68.8/3.4, exhaustiveness 16, seed 42).
Receptor: 5TZ1 chain A + heme, rigid.

## RESULT: GATE FAILED (exit code 1)

Pre-registered bar, declared before any results were seen: **AUC >= 0.70 AND EF@1% >= 5.0x**.

| mode | AUC | EF@1% | BEDROC | verdict |
|---|---|---|---|---|
| A: docking score alone | 0.561 | 10.10x | 0.162 | AUC fails |
| B: hard 3.0 A metal cutoff, then score | 0.338 | 0.00x | 0.005 | worse than random |

The gate blocked the screen, which is the designed behaviour. No new molecules were
ranked on top of an unvalidated pipeline.

## Root cause analysis

**1. The 3.0 A hard cutoff is miscalibrated.** It discarded 4 of 8 actives. Docking has
roughly 1-2 A of positional error, so demanding a 3.0 A coordination geometry rejects true
binders that land at 3.0-3.5 A. Re-running the four failures at 4x search depth: isavuconazole
recovered (2.74 A), itraconazole 3.22 A, miconazole 3.32 A — two of those sit just outside an
arbitrary line.

**2. Enlarging the box without scaling exhaustiveness degraded sampling.** The earlier
per-ligand work used a 26 A box and got 6/9 iron-bound. This screen used 30 A (1.5x the
volume) at the same exhaustiveness 16 and got 4/8. Search effort must scale with search
volume; it did not.

**3. Posaconazole remains unreachable** (5.13 A even at exhaustiveness 64). Consistent with
induced fit: 5TZ1 was solved with the short ligand VT-1161, and a rigid receptor cannot open
the channel its long tail needs.

## The signal is real, and it is geometric

Across the full screen, closest nitrogen-to-iron distance achieved:

    actives (n=8)    median 3.03 A
    decoys  (n=396)  median 4.69 A

Actives approach the iron markedly closer than property-matched decoys. The information is
present; the binary cutoff destroyed it by turning a graded signal into a yes/no.

## EXPLORATORY — post-hoc, NOT a re-graded pass

These ranking methods were evaluated AFTER seeing the results. They cannot be claimed as a
gate pass. They are a hypothesis for the next pre-registered run.

| ranking method | AUC | EF@1% | BEDROC |
|---|---|---|---|
| docking score alone (pre-registered) | 0.561 | 10.10x | 0.162 |
| **iron-approach distance alone** | **0.818** | **30.30x** | **0.547** |
| combined rank (distance + score) | 0.763 | 0.00x | 0.217 |
| hard cutoff 3.0 A then score | 0.338 | 0.00x | 0.005 |
| hard cutoff 3.5 A then score | 0.496 | 0.00x | 0.049 |
| hard cutoff 4.0 A then score | 0.444 | 0.00x | 0.022 |

**Hypothesis:** for azoles against CYP51, how closely a molecule can bring a nitrogen to the
heme iron is a substantially better discriminator than Vina's affinity score (AUC 0.818 vs
0.561; EF@1% 30x vs 10x). Geometry beats the scoring function on this target.

## What must happen before this can be claimed

1. **Pre-register** "rank by iron-approach distance" as mode C, with its threshold declared
   in advance, and re-run. Tuning on the data you then report is p-hacking.
2. **Validate on held-out actives.** The 8 azoles used here shaped the hypothesis; they cannot
   also test it. Candidate held-out set: oteseconazole (VT-1161, already have the crystal pose),
   efinaconazole, luliconazole, tioconazole, econazole, sertaconazole.
3. **Scale exhaustiveness with box volume** and fix one protocol before re-running.
4. **Report posaconazole-class ligands as a known blind spot** rather than silently scoring them.

## Honest summary

The pipeline as configured cannot yet re-discover the drugs we already know work, and it
said so out loud instead of producing a confident ranking. That is the system functioning
correctly. The run also produced a strong, testable hypothesis about why — and the fix is
cheap. Screening new molecules remains unjustified until a pre-registered re-run passes.

## Reproduce

    ./scripts/run_screen.sh
    conda run -n openafr --no-capture-output python scripts/validate_gate.py \
      work/screen work/receptor_A.pdb        # exits 1 on failure
