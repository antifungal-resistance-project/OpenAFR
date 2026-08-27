# Results — seed-stability of the candidate shortlist (issue #68)

Date: 2026-08-26. Pre-registered in
[PREREGISTRATION_shortlist_consensus.md](PREREGISTRATION_shortlist_consensus.md)
(sha256 `bab4c682…`, frozen before any candidate distance was read). Protocol hash `6eca16b1…`,
consensus criterion from [PREREGISTRATION_consensus.md](PREREGISTRATION_consensus.md). Fresh
5-seed docking of the shortlist against the rigid wild-type *C. albicans* CYP51
(`work/receptor.pdbqt`, 5TZ1 chain A + heme).

> **Outcome: the shortlist survives re-docking, but two hits do not.** 8 of the 10 novel
> survivors coordinate the heme iron reproducibly across all 5 seeds; **LDN-27219 and olprinone
> are demoted as single-dock artifacts** — exactly the failure mode this run exists to catch.
> The single-seed ranking that produced the shortlist mis-ranked at least one hit by pure
> search luck.

## What was run

16 molecules — the 10 novel shortlist survivors, 5 reference azoles from the same top-100
(flutrimazole, ravuconazole, letrozole, fluconazole, voriconazole), and ketoconazole — each
docked under 5 seeds {42, 1, 2, 3, 4}. Only the RNG seed varied; box, exhaustiveness (32),
num_modes (20) unchanged from the frozen protocol. Seed 42 reproduces the single-search number
the shortlist was built on. Consensus = min N-Fe over seeds; spread = max−min over the seeds
that coordinated (< 3.0 Å). Stable = coordinates in **all 5 seeds** AND spread ≤ 1.0 Å.

## The table (fungal receptor)

| candidate | role | fung min (Å) | spread | seeds coord. | **stable** |
|---|---|---:|---:|:---:|:---:|
| flutrimazole | azole ref | 2.371 | 0.017 | 5/5 | yes |
| verinurad | novel #2 | 2.417 | 0.027 | 5/5 | **yes** |
| vatalanib | novel #3 | 2.479 | 0.086 | 5/5 | **yes** |
| flucloxacillin | novel #4 | 2.517 | 0.041 | 5/5 | **yes** |
| taranabant | novel #5 | 2.564 | 0.032 | 5/5 | **yes** |
| lersivirine | novel #6 | 2.572 | 0.055 | 5/5 | **yes** |
| L-838417 | novel #10 | 2.580 | 0.087 | 5/5 | **yes** |
| LDN-27219 | novel #7 | 2.603 | 0.057 | **2/5** | **NO** |
| R-1479 | novel #8 | 2.621 | 0.225 | 5/5 | **yes** |
| ravuconazole | azole ref | 2.653 | 0.044 | 5/5 | yes |
| voriconazole | azole ref | 2.657 | 0.142 | 5/5 | yes |
| nolatrexed | novel #14 | 2.677 | 0.015 | 5/5 | **yes** |
| olprinone | novel #13 | 2.680 | 0.040 | **4/5** | **NO** |
| letrozole | azole ref | 2.705 | 0.026 | 5/5 | yes |
| fluconazole | azole ref | 2.711 | 0.047 | 5/5 | yes |
| ketoconazole | azole ref | 2.921 | 0.000 | **1/5** | **NO** |

## What this establishes

1. **8 of 10 novel survivors are seed-stable.** verinurad, vatalanib, flucloxacillin,
   taranabant, lersivirine, L-838417, R-1479, nolatrexed each reach the iron in all 5 seeds
   with small spreads (0.015–0.225 Å). Their shortlist coordination was not a single-dock
   fluke — this is the reinforcement the shortlist needed.

2. **LDN-27219 is a single-dock artifact — the headline catch.** Its shortlist rank (#7) came
   from seed 42's 2.603 Å, but its per-seed fungal distances are
   **{42: 2.603, 4: 2.66, 1: 5.574, 2: 5.509, 3: 5.566}** — it coordinates in only 2 of 5
   searches; in the other three the nitrogen never approaches the iron (~5.5 Å). One lucky seed
   put it on the shortlist. This is precisely the instability the reliability result warned of
   (a single dock is a noisy estimator; clotrimazole once ranked #2774/2776), now caught before
   any bench cost.

3. **olprinone is borderline (4/5).** It coordinates in four seeds (~2.68–2.72 Å) but misses in
   one (3.726 Å). Less alarming than LDN-27219, but it does not clear the pre-registered
   all-5-seed bar, so it is demoted per the frozen rule — honestly, as a borderline case, not a
   clean fail.

4. **Rigid azoles are rock-stable; the one flexible azole is not.** All five small reference
   azoles coordinate 5/5. Ketoconazole (36 heavy atoms, 8 rotatable bonds) coordinates in only
   1/5 — a docking-convergence limit for large flexible ligands from a single embedded
   conformer, which recurs on the human side (see the selectivity run) and matters for
   interpreting that run's positive control.

## The reinforced shortlist (seed-stable, this run)

verinurad, vatalanib, flucloxacillin, taranabant, lersivirine, L-838417, R-1479, nolatrexed.
**Demoted:** LDN-27219 (single-dock artifact, 2/5), olprinone (borderline, 4/5).

Selectivity against human CYP51 is assessed separately in
[RESULTS_human_selectivity.md](RESULTS_human_selectivity.md); the combined confidence view
lives in `work/candidates/shortlist_confidence.json`.

## Limitations

1. R = 5 reduces but does not eliminate search variance.
2. Rigid single-conformer *C. albicans* receptor, ≤45-heavy-atom domain — this run tests
   seed-stability only, not the conformer (#70) or *C. auris* receptor (#69) caveats.
3. Stable coordination is **necessary, not sufficient**: property-matched decoys also
   coordinate (the decoy ceiling). This run removes single-dock luck, not the ceiling.

## Reproduce

```
conda activate openafr
python scripts/prep_receptor.py                       # rebuild work/receptor.pdbqt
python scripts/prep_ligands.py work/candidates/shortlist_focus.smi work/candidates/focus_ligands
RECEPTOR=work/receptor.pdbqt scripts/screen_consensus.sh \
    work/candidates/focus_ligands work/candidates/dock_fungal "42 1 2 3 4"
python scripts/shortlist_confidence.py \
    --fungal work/candidates/dock_fungal --human work/candidates/dock_human
```

Docked poses are gitignored; the SMILES set (`work/candidates/shortlist_focus.smi`) and the
summary (`work/candidates/shortlist_confidence.json`) persist.
