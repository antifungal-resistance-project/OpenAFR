# Pre-registration — ensemble within-azole CONFIRMATION on an independent active set

**Status: FROZEN 2026-08-26, before any active/inactive label was scored under the ensemble.**
The independent active set was built and sampled first (deterministic, label-free — its sha256
is pinned below), then this document is frozen by its own SHA-256. The grader
(`validate_gate_ensemble_confirm.py`) re-checks this document's hash, `protocol.yaml`, the three
receptor-anchor hashes (unchanged from look #8), and refuses to certify a pass if any moved. No
AUC of the independent actives against the inactive pool had been read at freeze time.

## Why this run is required (it is the price look #8's PASS named)

Look #8 (`work/RESULTS_ensemble.md`) passed the within-azole gate at **ensemble-min mode C
AUC 0.750** — but on only the **7 held-out azoles**, with a bootstrap 95% CI of **0.683–0.825**
that straddles the 0.70 bar. Its pre-committed PASS branch said, in advance:

> "The within-class triage product re-opens, still at n=7 with a wide CI; **a separately
> pre-registered confirmation on an *independent* active set is required before any within-class
> wet-lab claim.**"

This is that confirmation. It answers exactly one question: **was the 0.750 a property of those
particular 7 molecules, or does the ensemble lift generalise to a broad, independent population
of measured-active azoles?**

## Single-variable design — everything is held fixed from look #8 except the actives

| | |
|---|---|
| ensemble | the EXHAUSTIVE 3-conformer *C. albicans* CYP51 set — 5TZ1, 5FSA, 5V5Z — receptor anchors **hash-pinned, identical to look #8** |
| criterion | mode C (min ligand-N-to-heme-Fe distance), **unchanged** |
| aggregation | **ensemble-minimum**, unchanged |
| protocol | `protocol.yaml` unchanged — box 26 Å @ 70.61/66.28/4.18, exhaustiveness 32, 20 modes, seed 42 |
| inactives | the **same frozen ~162 azole-bearing verified inactives** as looks #6/#7/#8 (`data/ligands/verified_inactives.smi`, `azole_bearing()`) |
| **actives** | **CHANGED** — an independent measured-active azole set (below), replacing the 7 held-out azoles |

Only the actives change. So a difference from look #8 is attributable to the active sample, not
to method drift.

## The independent active set — built by a frozen rule, not hand-picked

Built by `scripts/make_verified_actives.py` as the **single-variable mirror** of the
verified-inactive construction (`openafr/inactives.py`): from the SAME ChEMBL *C. albicans*
bioactivity cache (`calbicans_cell` 87,879 records + `calbicans_cyp51` 816 records, 43,004
distinct compounds; manifest sha in `work/chembl/MANIFEST.json`), keeping a compound iff **all**:

    * verdict == has-active-evidence   (>= 1 active record; openafr.inactives.compound_verdict.
                                        Its asymmetry — one active outranks any inactives — is
                                        CONSERVATIVE here too: a false active ranks low and only
                                        DEPRESSES measured separation)
    * azole-bearing                    (the same AZOLE_SMARTS the gate subgroup uses — a
                                        within-class test: measured-active vs measured-inactive azoles)
    * <= 45 heavy atoms                (the frozen run-2 applicability domain)
    * >= 1 aromatic nitrogen           (the mechanism requirement)
    * max Morgan-r2 Tanimoto < 0.70    (vs all 8 training + 7 held-out actives AND oteseconazole/
                                        VT-1161, the 5TZ1 co-crystal ligand — genuinely distinct
                                        chemistry; exact-matches of the 15 and the reference
                                        ligand are excluded, so no self-docking active can leak in)

Funnel (`make_verified_actives.py build`, 2026-08-26): 9,270 has-active-evidence compounds →
reject 5,310 non-azole, 261 over-domain, 46 too-similar-to-existing, 1 unparseable → **3,609
independent measured-active azoles** (heavy atoms 9–45, mean 31.8). Population written to
`data/ligands/verified_actives.smi` (sha256 `1ee6e45e0c550583…`).

**The frozen test set (what this run grades).** 3,609 is far more than can be docked; the test
set is a **deterministic random sample of N = 200, seed 42** (`make_verified_actives.py sample
--n 200 --seed 42`) — unbiased (a seeded shuffle over the sorted population, not a head-of-list
slice), reproducible, and ~30× the original n = 7. Written to
`data/ligands/verified_actives_sample.smi`, **sha256
`078c247f742a773744e968d7cbef7d330e60837c76eb039653dd017ec08ae46a`** — pinned here and verified
by the grader. Composition (label-blind, from structure only): 200/200 azole-bearing, heavy
atoms 16–44 (mean 32.0), all in-domain.

## Primary hypothesis (H1) — the only thing that gates

Ranking the pooled molecules (the 200-azole independent sample + the ~162 azole-bearing measured
inactives) by **ensemble-min mode C over {5TZ1, 5FSA, 5V5Z}** separates measured-active azoles
from measured-inactive azole analogues at **AUC ≥ 0.70** (`protocol.yaml` `min_auc`, unchanged).
AUC is the gate because it is threshold-free and prevalence-independent (the active/inactive
balance differs from look #8; AUC is unaffected, which is why it, not EF, gates).

## Secondary — all pre-specified, none gating

1. **Per-conformer AUC** (5TZ1 / 5FSA / 5V5Z alone) and **ensemble-mean** AUC — does the lift
   still ride on 5FSA, as in look #8?
2. **Marginal-conformer increment** {5TZ1} → {+5FSA} → {+5V5Z}.
3. **Novel-chemotype guard** — the 200 independent actives vs the **non-azole** inactives
   (comparator 0.810 single-conformer / 0.949 look #8). Must stay ≥ ~0.78 to call a real win.
4. **Statistical support**: permutation p (20,000 shuffles), bootstrap 95% AUC CI, the tip.
   With n = 200 the CI will be *tight* — that tightness (not just the point estimate) is the
   whole reason this run exists.
5. **EF@1/5/10%, BEDROC(20)** reported for continuity, but with ~200 actives vs ~162 inactives
   the pool is near-balanced, so EF caps near ~1.8× and is a weak statistic in this direction
   (the opposite problem from look #8's n = 7); it does not gate.

## Pre-committed interpretation

- **PASS — AUC ≥ 0.70 AND guard ≥ ~0.78.** The within-azole lift **replicates** on a broad,
  independent, measured-active azole population. Look #8's 0.750 was not a fluke of the 7 held-out
  azoles; the ensemble genuinely ranks measured-active azoles above measured-inactive analogues.
  Within-class triage is confirmed **at the ensemble level** — with the standing caveat that this
  is a *triage* AUC, not top-of-list enrichment (look #8's EF@1/5% were 0.00×; read this run's
  EF/BEDROC and the tip before any wet-lab prioritisation). The candidate shortlist
  (`RESULTS_candidate_shortlist.md`) can now be re-scored under the confirmed ensemble.
- **PASS on within-azole but the guard collapses.** The ensemble discriminates within-class only
  by making everything coordinate more easily. Reported honestly as not a clean win, with both
  numbers; the single rigid conformer remains the product receptor.
- **FAIL — AUC < 0.70.** The within-azole lift does **not** replicate. Look #8's 0.750 was
  consistent with a 7-molecule sampling fluctuation (its CI straddled the bar). **No within-class
  claim is earned**; the defensible product remains novel-chemotype triage (AUC 0.810), and the
  within-azole line closes for good on the rigid *and* the exhaustive-ensemble receptor.

## Stopping rule (binding)

**One independent sample, one look.** The test set is the frozen N = 200, seed 42 sample — it
will not be re-sampled with a different seed, N, or novelty cutoff and re-graded if it
disappoints. The criterion (mode C), aggregation (ensemble-min), ensemble (the 3 pinned
conformers), and inactive pool are all unchanged from look #8. Any further within-class work is a
new, separately pre-registered dataset.

## Limitations (fixed in advance)

1. **Whole-cell "active" ≠ CYP51 "active".** Most ChEMBL *C. albicans* actives are whole-cell
   MICs; a whole-cell active could act off-target. This is the mirror of the inactive-set
   confound and is conservative in the same way (a mislabelled active ranks low, depressing AUC).
2. **Single active record suffices for `has-active-evidence`.** Some of the 3,609 may be weakly
   or singly-reported active; the random sample inherits that. Conservative (depresses AUC), and
   representative of the measured-active azole population this run is meant to generalise to.
3. **All three conformers are holo** (azole-bound), a thumb on the scale toward actives —
   unchanged from look #8, carried here.
4. **Near-balanced pool** makes EF/BEDROC weak in this run; AUC is the honest, prevalence-neutral
   read, which is why it alone gates.
5. **A rigid-docking false negative** (an in-domain azole whose tail is modelled folded back,
   losing its coordinating pose) ranks LAST, never dropped — conservative for the gate.

## Reproduce

```
conda activate openafr
python scripts/make_verified_inactives.py fetch                       # ChEMBL cache (network)
python scripts/make_verified_actives.py build                         # -> verified_actives.smi (3609)
python scripts/make_verified_actives.py sample --n 200 --seed 42      # -> verified_actives_sample.smi
python scripts/prep_ligands.py data/ligands/verified_actives_sample.smi work/verified_ligands
RECEPTOR=work/receptor.pdbqt      scripts/screen.sh work/verified_ligands work/screen_verified
RECEPTOR=work/receptor_5FSA.pdbqt scripts/screen.sh work/verified_ligands work/screen_5FSA
RECEPTOR=work/receptor_5V5Z.pdbqt scripts/screen.sh work/verified_ligands work/screen_5V5Z
python scripts/validate_gate_ensemble_confirm.py                      # exits 0 pass / 1 fail
```
