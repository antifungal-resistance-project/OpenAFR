# Channel engagement — the within-azole question (mode F, look #7)

Date: 2026-08-25
Pre-registration: `work/PREREGISTRATION_channel_engagement.md` (sha256 `62592a101ad893f1…`,
verified unmodified by the grading script at run time)
Test set: 7 held-out azoles + ~162 azole-bearing verified inactives (169 ranked)
Protocol: unchanged and hash-verified — box 26 Å @ 70.61/66.28/4.18, exhaustiveness 32,
num_modes 20, seed 42, rigid 5TZ1 chain A + heme. Feature: mode F (count of distinct
channel-shell residues contacted in the coordinating pose; channel-shell defined label-blind
from VT1 co-crystal, 22 residues within 4.5 Å).
Host: Apple Silicon, `openafr` env, ~2.5 h wall for 294 dockings on 12 cores.

## RESULT: GATE FAILED — the within-azole ceiling is not iron-distance-specific

**H1 fails: AUC 0.590 < 0.70.** An orthogonal, mechanistically-motivated tail-fit feature
also cannot rank working azoles above measured-inactive azole analogues on a single rigid
conformer. Mode F is *worse* than mode C on this pool (0.590 vs 0.650).

| ranking | pool | AUC | EF@5% | BEDROC | actives in top 15 |
|---|---|---|---|---|---|
| **mode F — channel engagement** | 7 actives + 162 azole-bearing inactives | **0.590** | 2.68x | 0.179 | 1/7 |
| mode F — powered (n=15, secondary) | 15 actives + 162 inactives | 0.623 | 3.93x | 0.295 | 3/15 |
| mode C comparator (look #6, same pool) | 7 + 162 | 0.650 | 2.68x | 0.108 | 1/7 |

Mode F is genuinely orthogonal to mode C (r = −0.301), so this is not a failure of
independence — it is a different signal that is also insufficient. The size-residualized AUC
(0.612) is similarly below the bar, so the failure is not explained by mode F merely
re-measuring molecular weight.

## Orthogonality and size diagnostics

    r(mode F, mode C = N–Fe distance) : −0.301   (genuinely different signal)
    r(mode F, heavy-atom count)        : +0.428   (moderate size correlation)
    size-residualized AUC              :  0.612   (size does not explain the failure)

The mechanistic claim (tail occupancy is the within-class discriminator) is not refuted — but
it cannot be evaluated on a single rigid conformer. A tail that would thread the channel under
induced fit may be modelled folded back by rigid docking, which would blunt mode F through no
fault of the hypothesis (pre-registration Limitation 2).

## Statistical support (descriptive; the gate is the AUC point estimate)

    permutation test (20,000 shuffles): p = 0.2192   (not significant)
    bootstrap 95% CI over actives     : AUC 0.359 – 0.785
    best active                       : fenticonazole, rank 2/169 (1.18th percentile)
    inactives above it                : 1

The p-value is non-significant and the bootstrap CI lower bound is well below the bar,
consistent with the point-estimate failure.

## Pre-committed interpretation (from pre-registration, §Pre-committed interpretation)

> FAIL (AUC < 0.70). The within-azole ceiling is **not specific to iron distance**: an
> orthogonal, mechanistically-motivated tail-fit feature also cannot separate working from
> failed azoles on a single rigid conformer. This strengthens the preprint's §4.3 induced-fit /
> single-conformer limitation and is the decisive evidence for **Option B** — re-scope the
> product to novel-chemotype triage, where the measured AUC is 0.810. The next attack on the
> within-class question would then have to be a *flexible / ensemble receptor* (new poses = a
> new dataset), separately pre-registered.

**Stopping rule is binding.** Mode F will not be re-defined with a different contact cutoff,
a different channel-shell definition, atom-counts instead of residue-counts, or a fitted
weight. One orthogonal feature was spent on the verified-inactive poses; no third feature will
be attempted on these rigid conformers.

## What this means for the product

The defensible product is **novel-chemotype triage** (non-azole chemotype pool, AUC 0.810,
EF@5% 4.98x). A shortlist this tool produces will be enriched for novel-scaffold compounds
that coordinate the heme iron — which is exactly what a wet-lab collaborator would want (not
re-ranked azoles that have already been tried). Within-class azole ranking requires
ensemble/flexible-receptor poses and a separately pre-registered attack; that is the correct
next direction if within-class triage is the goal.

## Reproduce

```
conda activate openafr
python scripts/prep_receptor.py
python scripts/prep_ligands.py data/ligands/actives_holdout_final.smi work/verified_ligands
python scripts/prep_ligands.py data/ligands/verified_inactives.smi   work/verified_ligands
python scripts/prep_ligands.py data/ligands/actives.smi              work/verified_ligands
scripts/screen.sh work/verified_ligands work/screen_verified
python scripts/validate_gate_channel.py work/screen_verified   # exits 1 (FAIL)
```
