# Active-power gate — the n=7 tip is retired: AUC lower bound clears 0.70 at N=300

Date: 2026-09-06
Pre-registration: `work/PREREGISTRATION_active_power.md` (sha256 `9327563ac69ff087…`,
verified unmodified by the grader at run time). Test set:
`data/ligands/active_power_sample.smi` (sha256 `ba8d1f3cb0571701…`, N=300, seed 7),
verified by the grader. Protocol `protocol.yaml` (sha256 `6eca16b1be2e5e33…`) unchanged.
Compute: GCP e2-standard-16 (project `openafr-recaller-20260818`, us-central1-a, VM
deleted after), conda env from `environment.yml` (vina 1.2.7, openbabel 3.1.0, rdkit).
Grader: `scripts/validate_gate_active_power.py`.

Single-variable widening of the Run-2 primary gate (`work/RESULTS_run2_holdout.md`): rigid
5TZ1, mode C (minimum nitrogen-to-heme-iron distance), the same 350 property-matched decoys,
the frozen protocol/seed — **only the actives changed**, from 7 held-out azoles to a blind
N=300 sample of measured-active azoles. Docking: 292/300 actives and 348/350 decoys posed
(8 actives + 2 decoys were obabel/Vina failures, ranked LAST per pre-registration — never
dropped, which is conservative for the gate).

## RESULT: GATE PASSED — Limitation #1 RETIRED

**H1 passes on both counts: mode C AUC 0.750, AND the bootstrap 95% CI lower bound 0.717
≥ 0.70.** The n=7 tip is discharged: with the sampling noise of seven molecules removed,
the geometry criterion's held-out AUC clears the project's 0.70 usefulness bar with its
*lower bound*, not merely its point estimate.

| contrast | pool | AUC | bootstrap 95% CI | perm p |
|---|---|---|---|---|
| **PRIMARY — actives vs property-matched decoys** | 300 act + 350 decoys | **0.750** | **0.717 – 0.782** | 5e-05 |
| Run 2 (n=7, same contrast) — for reference | 7 act + 348 decoys | 0.794 | 0.590 – 0.946 | 0.0028 |
| S1 product — actives vs NON-azole verified inactives | 300 act + 116 inact | 0.774 | 0.741 – 0.807 | 5e-05 |

The Run-2 point estimate (0.794) always rested on 7 molecules with a CI spanning 0.590–0.946.
The widened sample lands at 0.750 with a CI of **0.717–0.782** — the point estimate is a
little lower, but the interval is ~4× tighter and its lower bound sits **above** the bar. The
honest headline is the well-powered interval, not the n=7 point estimate.

## Secondaries (all pre-specified, none gating)

- **S1 — the product claim.** Novel-chemotype triage (measured-active azoles vs NON-azole
  verified inactives) holds at **AUC 0.774, CI 0.741–0.807**, well-powered and consistent
  with the single-conformer comparator 0.810 and the look #9 guard 0.823. The defensible
  product AUC is now backed by a tight interval, not a 7-active estimate.
- **S4 — size-artifact check.** Splitting the actives at the median heavy-atom count (32):
  small half AUC 0.766, large half AUC 0.732. Both halves clear the bar; the signal is not
  carried by molecular size.
- **EF is not gating and is weak here by construction** (near-balanced-to-modest pools):
  primary EF@1% 2.17×, S1 EF@1% 1.39×. Consistent with the preprint §3.2 finding that EF is
  seed-fragile and prevalence-sensitive; AUC + its CI are the honest, prevalence-neutral read.

## Pre-committed interpretation (from the pre-registration, verbatim branch)

> **RETIRED — H1 AUC ≥ 0.70 AND CI lower bound ≥ 0.70.** Limitation #1 is discharged: with
> the sampling noise of n = 7 removed, the geometry criterion's held-out AUC clears the
> usefulness bar with its lower bound, not just its point estimate. The preprint §6.1 is
> updated from "a true AUC under the bar cannot be excluded" to the measured tight interval.
> The Run-2 point estimate (0.794) stands; this run reports whether it is robustly above 0.70.

## What this earns

Preprint Limitation #1 ("n = 7 actives … a true AUC under the bar cannot be excluded")
updates to a measured, well-powered statement: on a blind N=300 sample of measured-active
azoles, the geometry criterion separates them from property-matched decoys at AUC 0.750
(95% CI 0.717–0.782), and from non-azole chemotypes at 0.774 (0.741–0.807). The single
largest small-n caveat on the headline result is retired. This does NOT re-open the
within-azole line (closed by look #9) and makes no wet-lab claim; it hardens the existing
novel-chemotype triage product.

## Reproduce

Docking runs on a Linux/x86 VM with Vina; see `.context/RUNBOOK_active_power.md`.

```
conda activate openafr
python -m openafr.protocol --verify
python scripts/prep_receptor.py
python scripts/prep_ligands.py data/ligands/active_power_sample.smi work/ap_actives
python scripts/prep_ligands.py data/ligands/decoys_holdout.smi      work/ap_decoys
python scripts/prep_ligands.py data/ligands/verified_inactives.smi  work/ap_inactives
RECEPTOR=work/receptor.pdbqt scripts/screen.sh work/ap_actives   work/screen_active_power
RECEPTOR=work/receptor.pdbqt scripts/screen.sh work/ap_decoys    work/screen_active_power
RECEPTOR=work/receptor.pdbqt scripts/screen.sh work/ap_inactives work/screen_active_power
python scripts/validate_gate_active_power.py work/screen_active_power work/receptor_A.pdb   # exits 0
```
