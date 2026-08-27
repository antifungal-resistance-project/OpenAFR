# Results — candidate resistance-pocket retention (*C. auris* Y132F)

Date: 2026-08-27. Pre-registered in
[PREREGISTRATION_candidate_resistance.md](PREREGISTRATION_candidate_resistance.md)
(sha256 2241b85a…, frozen BEFORE any Y132F distance was read;
`work/PREREG_candidate_resistance.sha256`). Protocol hash unchanged from run 2.

> **This is a per-candidate confidence axis, not a gate.** It demotes/annotates the
> scorecard; it never promotes. It answers the first `not-yet-checked` line on every card
> (#69, #77): do the shortlist candidates keep their iron-coordinating geometry in the
> resistant pocket the mission actually targets, or only in the *C. albicans* stand-in they
> were ranked against?

## What was tested

The 16 scorecard molecules (10 novel candidates + 6 known-azole controls) were docked into
the *C. auris* **Y132F** resistance pocket (`work/receptor_auris_Y132F.pdbqt`, a pure para-
hydroxyl deletion on the validated 5TZ1 chain-A scaffold — same heme, same iron, same box)
under the frozen protocol, 5 seeds {42,1,2,3,4}. The **only** variable versus the wild-type
candidate dock is the residue-132 side chain. The wild-type baseline is the already-frozen
`fungal_consensus` in `work/candidates/shortlist_confidence.json` (not re-docked).

Verdicts (frozen): **RETAINED** = coordinates Y132F iron ≥3/5 seeds and |shift| ≤ 0.5 Å;
**SHIFTED** = ≥3/5 but |shift| > 0.5 Å; **LOST** = <3/5 seeds reach the iron.

## The result

| candidate | WT N–Fe | Y132F N–Fe | coord | shift (Å) | verdict |
|---|---:|---:|:--:|---:|---|
| flutrimazole* | 2.371 | 2.358 | 5/5 | −0.013 | RETAINED |
| verinurad | 2.417 | 2.369 | 5/5 | −0.048 | RETAINED |
| **vatalanib** | 2.479 | 2.510 | **2/5** | — | **LOST** |
| flucloxacillin | 2.517 | 2.538 | 4/5 | +0.021 | RETAINED |
| taranabant | 2.564 | 2.534 | 5/5 | −0.030 | RETAINED |
| lersivirine | 2.572 | 2.636 | 5/5 | +0.064 | RETAINED |
| L-838417 | 2.580 | 2.555 | 5/5 | −0.025 | RETAINED |
| LDN-27219 | 2.603 | 2.586 | **1/5** | — | **LOST** |
| R-1479 | 2.621 | 2.677 | 5/5 | +0.056 | RETAINED |
| ravuconazole* | 2.653 | 2.644 | 5/5 | −0.009 | RETAINED |
| voriconazole* | 2.657 | 2.656 | 5/5 | −0.001 | RETAINED |
| nolatrexed | 2.677 | 2.644 | 5/5 | −0.033 | RETAINED |
| olprinone | 2.680 | 2.667 | 5/5 | −0.013 | RETAINED |
| letrozole* | 2.705 | 2.674 | 5/5 | −0.031 | RETAINED |
| fluconazole* | 2.711 | 2.697 | 5/5 | −0.014 | RETAINED |
| ketoconazole* | 2.921 | 2.703 | 3/5 | −0.218 | RETAINED |

*controls (known azoles). Tally: **14 RETAINED, 2 LOST**.

## What this shows that no earlier axis did

- **vatalanib is a genuinely new demotion.** In the wild-type pocket it was one of the
  tightest, fully-stable coordinators (2.479 Å, **5/5** seeds, widest human-selectivity margin
  +1.82 Å) — it carried only one open concern (a hERG flag). In the Y132F pocket it coordinates
  in only **2/5** seeds (seeds 3,4 at ~2.51 Å; seeds 1,2,42 fall to 3.3–4.2 Å): the resistance
  substitution turns its reproducible coordination into search luck. **No prior axis — geometry,
  seed stability, selectivity, ADMET — surfaced this.** It is exactly the failure the mission
  cares about: a molecule that looks solid on the *C. albicans* stand-in but whose geometric
  basis is fragile in the resistant target.
- **The other LOST is unsurprising and consistent.** LDN-27219 was already seed-unstable in
  wild-type (2/5); its 1/5 on Y132F confirms it was never a reliable coordinator, in either
  pocket — not a new resistance signal, just the same instability.
- **The two scorecard front-runners survive.** lersivirine and L-838417 — the only novel hits
  clean on every previously-measured axis — both **RETAINED** (5/5, shifts +0.064 / −0.025 Å).
  The resistance axis does not open a new concern on either, so the scorecard's two zero-concern
  hits stay zero-concern.
- **Retained molecules barely move.** Every RETAINED candidate shifts <0.1 Å (except ketoconazole,
  whose WT baseline was a single poorly-converged seed — see below). This reproduces the held-out
  finding that the geometric criterion is **nearly insensitive** to Y132F.

## The honesty caveat (carried from the held-out Y132F run — do NOT overclaim)

**RETAINED ≠ resistance-breaking.** Because the criterion is nearly insensitive to the Y132F
change ([RESULTS_auris_Y132F.md](RESULTS_auris_Y132F.md): AUC 0.794 → 0.805 on the held-out
azoles), a RETAINED verdict means only that the coordination geometry which put a molecule on
the shortlist is **not destroyed** by the dominant resistance substitution — a robustness
check, not evidence the molecule *defeats* resistance. The known-azole controls RETAINING is
the clearest proof of this: those drugs are clinically defeated by Y132F, yet they still rank
well on the mutant pocket. Genuine resistance-breaker discovery needs actives *measured against
the resistant strain* (the T1 active-set expansion), which this is not. **The informative
direction of this axis is the LOST verdict**, which flags a shortlist molecule whose geometry
collapses in the pocket the mission targets — and vatalanib is one.

- **ketoconazole's improvement is not signal.** It reads 1/5 → 3/5 (shift −0.218 Å), but its
  wild-type 1/5 was a known floppy-ligand global-search artifact
  ([RESULTS_human_control.md](RESULTS_human_control.md)), not a receptor effect. The apparent
  "gain" is the WT baseline being a single poorly-converged seed, not Y132F helping. Reported,
  not leaned on.

## Limitations (all from the pre-registration, all still bind)

1. Rigid single-atom deletion under-models real Y132F resistance (no conformational repacking,
   no efflux/expression changes).
2. **Y132F only.** K143R (clade I) and F126L (clade III) need a side-chain repacker not in the
   pinned environment — still `not-yet-checked`.
3. Single rigid conformer; coordination necessary-not-sufficient.

## Reproduce

```
conda activate openafr
python scripts/prep_receptor.py --mutant work/receptor_auris_Y132F.pdb work/receptor_auris_Y132F.pdbqt
python scripts/prep_ligands.py work/candidates/shortlist_focus.smi work/candidates/cand_ligands
RECEPTOR=work/receptor_auris_Y132F.pdbqt scripts/screen_consensus.sh \
    work/candidates/cand_ligands work/candidates/dock_auris_Y132F
python scripts/candidate_resistance.py --json > work/candidates/shortlist_resistance.json
python scripts/candidate_resistance.py    # human-readable table
```
