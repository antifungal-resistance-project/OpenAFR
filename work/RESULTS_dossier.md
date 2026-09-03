# Candidate dossiers — the fused, machine-readable artifact

> **UPDATE 2026-09-02 — regenerated under the coordinator-identity gate; tally is now
> A=0, B=14, C=2 (was A=1, B=13, C=2).** After this dossier first ran, hardening its sole
> Priority-A candidate (verinurad) exposed that the geometry criterion measures the nearest
> nitrogen of *any* type, not the aromatic ring N the azole mechanism was validated on (see
> [RESULTS_coordinator.md](RESULTS_coordinator.md)). The dossier now gates automatic-A on the
> coordinator's identity — the sixth A-condition, `in-mechanism coordination`, requires the
> nitrogen reaching the iron to be an aromatic ring N. Both a **dual-mode** rank (inflated by
> a non-aromatic N, aromatic N still in band) and an **out-of-mechanism** rank (set by a
> nitrile/amine/azide N, aromatic N not in band) therefore **cannot earn automatic-A**, and
> both are flagged prominently — **but neither is hard-demoted to C.** The pre-registered gate
> re-validation ([RESULTS_aromatic_criterion.md](RESULTS_aromatic_criterion.md)) found that a
> real active, **luliconazole**, itself coordinates via its nitrile in-pose (aromatic N 5.3 Å),
> so a non-aromatic coordinator is *unvalidated, not disproven* — forcing it to C would risk
> rejecting a genuine active. Consequence: **verinurad → B (dual-mode)**; **taranabant,
> lersivirine, LDN-27219, olprinone → B, flagged out-of-mechanism**; **no molecule clears
> automatic-A.** Ordering is unchanged (still N–Fe, the frozen criterion — the coordinator gate
> annotates priority, it never re-ranks). All 6 azole controls remain in-mechanism.

Date: 2026-09-02

The first time the dossier generator (`openafr/dossier.py` + `scripts/build_dossier.py`,
PR #118) is actually **run**. It fuses every axis the pipeline computes — geometry rank,
assay precedent, drug-likeness (ADMET), human-CYP off-target liability, and
synthesizability — into **one machine-readable dossier per molecule**, each carrying an
A/B/C experimental-priority verdict and a one-sentence, sourced selection hypothesis. This
is the artifact a wet-lab PI actually reads: the scattered `RESULTS_*` axes collapsed onto
one row per candidate.

Run over the **focus shortlist** (`work/candidates/shortlist_focus.smi`, 16 molecules: the
10 novel repurposing survivors + the azole reference controls), freshly **re-docked** under
the frozen protocol for this run (the 2026-08-05 screen poses are gitignored and do not
persist). Docking: `scripts/screen.sh`, receptor rebuilt deterministically from
`data/structures/5TZ1.pdb` (`prep_receptor.py`, anchor-matched to `work/receptor_A.pdb`),
box/exhaustiveness/seed from the hash-frozen `protocol.yaml`. Provenance stamp on every
dossier records `protocol_unmodified: true` — the anti-tuning guarantee held.

> **Every dossier is a hypothesis for a wet lab, never a hit.** Nothing here has been tested
> against any fungus by us. A/B/C is an *orthogonal annotation on top of the validated
> geometry ranking, never a re-rank* — priority never moves a molecule up or down the N-Fe
> order. The tool's defensible claim is novel-chemotype triage (AUC ~0.81), not within-azole
> potency ranking.

## The tally: A=0, B=14, C=2

`coord` = which nitrogen sets the reported N–Fe (RESULTS_coordinator.md): **in** = aromatic
ring N (validated mechanism); **dual** = non-aromatic N sets the rank but the aromatic N still
reaches the band; **out** = non-aromatic N, aromatic N not in band. Both dual and out **block
automatic-A** and are flagged; neither forces C (a non-aromatic coordinator is unvalidated, not
disproven — the active luliconazole also uses one, RESULTS_aromatic_criterion.md).

| rank | molecule | N–Fe (Å) | coord | priority | why (short) |
|---:|---|---:|:--:|:--:|---|
| 1  | flutrimazole   | 2.371 | in   | B | azole control; 1 Ro5 + 1 structural alert, imidazole pan-CYP warhead |
| 2  | **verinurad**  | 2.417 | dual | **B** | rank set by its **nitrile** N; aromatic pyridine N only 2.80 Å (band edge) — no longer automatic-A |
| 3  | vatalanib      | 2.479 | in   | B | in-mechanism; below top-decile + 1 Ro5 violation (not automatic-A) |
| 4  | flucloxacillin | 2.558 | in   | B | in-mechanism; 1 structural alert |
| 5  | taranabant     | 2.586 | out  | **B** | rank set by a **nitrile** N; aromatic N 6.9 Å from iron — flagged out-of-mechanism |
| 6  | lersivirine    | 2.599 | out  | **B** | rank set by a **nitrile** N (2 nitriles); aromatic N 4.1 Å — flagged out-of-mechanism |
| 7  | LDN-27219      | 2.603 | out  | **B** | rank set by an **amine** N; also single-dock unstable (#68) — flagged out-of-mechanism |
| 8  | R-1479         | 2.621 | out  | **C** | PAINS (azide/diazo) forces C; rank also set by the **azide** terminus |
| 9  | ravuconazole   | 2.653 | in   | B | azole reference (record-inconclusive precedent) |
| 10 | L-838417       | 2.667 | in   | B | in-mechanism, but below the top-decile geometry cutoff |
| 11 | olprinone      | 2.680 | out  | **B** | rank set by a **nitrile** N; also single-dock unstable (#68) — flagged out-of-mechanism |
| 12 | nolatrexed     | 2.692 | in   | B | in-mechanism; below top-decile (also hits human CYP51, #73) |
| 13 | letrozole      | 2.719 | in   | B | azole reference (record-inconclusive) |
| 14 | fluconazole    | 2.740 | in   | B | **reference: measured enzyme-active** (positive control) |
| 15 | voriconazole   | 2.799 | in   | B | **reference: measured enzyme-active** (positive control) |
| 16 | ketoconazole   | 2.921 | in   | **C** | PAINS + reference-active — a control, correctly deprioritised |

### No molecule clears automatic-A once the mechanism is enforced

The sole Priority-A (verinurad) was demoted to B: its rank-#2 2.417 Å is set by its nitrile
nitrogen, and automatic-A now requires the coordinating nitrogen to be the validated aromatic
ring N. Four more candidates whose rank is set by a nitrile/amine/azide nitrogen (taranabant,
lersivirine, LDN-27219, olprinone) are also blocked from A and carry a prominent
out-of-mechanism flag — but they stay **B, not C**: the gate re-validation showed a real active
(luliconazole) coordinating via a non-aromatic N, so hard-demoting these would risk rejecting a
genuine binder. The honest read: the fused pipeline, with the mechanism enforced, yields **no
automatic-A candidate** — the strongest hypotheses are the in-mechanism B tier (vatalanib,
flucloxacillin, L-838417, nolatrexed), each carrying its own named caveat (a Ro5 violation, a
structural alert, a below-decile rank, or a human-CYP hit). Only R-1479 and ketoconazole are C,
both on PAINS.

### The priority table validated itself on the controls
- **fluconazole / voriconazole** carry `precedent: tested-active` and are tagged
  **reference (positive control, not a novel hit)** — the table refuses to promote a known
  inhibitor to a discovery slot, exactly as designed.
- **ketoconazole** lands in C: it recovers its known PAINS/liability profile — the control
  proving the deprioritisation logic bites.
- **R-1479** → C on two independent rules now: the PAINS/ADMET axis (#75, azide/diazo
  chemically implausible) *and* the coordinator gate — its reported rank is set by the
  **azide** terminus, not an aromatic ring N. The fused artifact re-derives both from
  primitives, not from a hand-copied note.

## Why verinurad was demoted from Priority-A
> Ranks #2 by N-Fe geometry — but the coordinating nitrogen is its **nitrile**, not an
> aromatic ring N. Its pyridine N (the azole-mimetic coordinator) reaches only 2.80 Å.

In the first run verinurad was the sole automatic-A: it cleared top-decile geometry rank,
convergent poses, novelty, clean drug-likeness, and obtainability. The
"coordinator-fragile" flag (#87) worried its single aromatic N was a thin coordinator.
Hardening that flag showed the sharper truth ([RESULTS_coordinator.md](RESULTS_coordinator.md)):
its rank-#2 2.417 Å is set by its **nitrile** nitrogen in all 5 seeds — an out-of-mechanism
interaction the method was never validated on. Its aromatic pyridine N reaches only **2.790 Å**
(band edge), which on its own would rank verinurad ≈14th. Automatic-A now requires the
coordinating nitrogen to be the validated aromatic ring N (the sixth A-condition,
`in-mechanism coordination`), so verinurad is **dual-mode → B**: an azole-like binding mode
exists, but the geometry that made it the top pick was a nitrile artifact. (Note also: the
human-selectivity axis is no longer *provisional* — the #73 native-pose control validated the
receptor; verinurad is simply marginally selective there.)

**There is no Priority-A candidate in this regeneration.** The strongest hypotheses are the
in-mechanism B tier — **vatalanib** foremost (in-mechanism, seed-stable #68, cleanly
human-selective #73; held out of A only by a single Ro5 violation).

## An honest note on the run: the `--catalogued` provenance flag
On the first invocation every molecule showed `availability: not-checked`, which the
priority table reads as *not* physically obtainable, forcing **A=0** — a misleadingly
pessimistic result, since the entire focus library is Drug Repurposing Hub compounds already
documented as buyable (`RESULTS_synth.md` #86). The CLI hardcoded `catalogued=False` with no
override. This run adds a **`--catalogued`** flag: an explicit operator claim about library
provenance (not a per-molecule network lookup, so the dossier stays deterministic and
offline). Run **with** `--catalogued`, availability resolves to `catalogued-drug` and
verinurad reaches A. The flag defaults off, so existing behaviour is unchanged.

## What this does and does not establish
- **Does:** produce the fused, per-molecule deliverable the generator was built for; prove
  the A/B/C table on the fluconazole/voriconazole/ketoconazole controls; re-derive the R-1479
  demotion from primitives; and stamp `protocol_unmodified: true` provenance on every card.
- **Does not:** claim any molecule is a hit, re-rank anything (order stays N-Fe), or re-run
  the gate. All limitations of the underlying screen still bind (single rigid *C. albicans*
  5TZ1 conformer, coordination necessary not sufficient, per-molecule dock variance).

## Reproduce
    # 1. receptor + ligands (both gitignored outputs)
    conda run -n openafr python scripts/prep_receptor.py
    conda run -n openafr python scripts/prep_ligands.py \
        work/candidates/shortlist_focus.smi work/screen_ligands
    # 2. dock the 16 focus ligands under the frozen protocol -> work/screen/*.pdbqt
    conda run -n openafr bash scripts/screen.sh work/screen_ligands work/screen
    # 3. fuse into dossiers (precedent hits ChEMBL + PubChem live)
    conda run -n openafr python scripts/build_dossier.py \
        work/screen work/candidates/shortlist_focus.smi work/receptor_A.pdb \
        -o work/candidates/dossiers.json --split work/candidates/dossiers \
        --near data/ligands/verified_inactives.smi --catalogued

Outputs: `work/candidates/dossiers.json` (the array) and `work/candidates/dossiers/<name>.json`
(the per-molecule registry). Precedent calls hit ChEMBL (www.ebi.ac.uk) and PubChem live on
2026-09-02.
