# Candidate dossiers — the fused, machine-readable artifact

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

## The tally: A=1, B=13, C=2

| rank | molecule | N–Fe (Å) | priority | why (short) |
|---:|---|---:|:--:|---|
| 1  | flutrimazole   | 2.371 | B | azole control; 1 Ro5 + 1 structural alert, imidazole pan-CYP warhead |
| 2  | **verinurad**  | 2.417 | **A** | top-decile, convergent, novel, clean drug-likeness, obtainable |
| 3  | vatalanib      | 2.479 | B | 1 Ro5 violation (not automatic-A clean) |
| 4  | flucloxacillin | 2.558 | B | 1 structural alert |
| 5  | taranabant     | 2.586 | B | 2 Ro5 violations (insoluble, over size) |
| 6  | lersivirine    | 2.599 | B | clean, but below the top-decile geometry cutoff |
| 7  | LDN-27219      | 2.603 | B | 3 structural alerts (reactive hydrazide) |
| 8  | R-1479         | 2.621 | **C** | PAINS (azide/diazo) — assay-interference, unreliable readout |
| 9  | ravuconazole   | 2.653 | B | azole reference (record-inconclusive precedent) |
| 10 | L-838417       | 2.667 | B | clean, but below the top-decile geometry cutoff |
| 11 | olprinone      | 2.680 | B | below top-decile (also a demoted single-dock artifact, #68) |
| 12 | nolatrexed     | 2.692 | B | below top-decile (also hits human CYP51, #73) |
| 13 | letrozole      | 2.719 | B | azole reference (record-inconclusive) |
| 14 | fluconazole    | 2.740 | B | **reference: measured enzyme-active** (positive control) |
| 15 | voriconazole   | 2.799 | B | **reference: measured enzyme-active** (positive control) |
| 16 | ketoconazole   | 2.921 | **C** | PAINS + reference-active — a control, correctly deprioritised |

### The priority table validated itself on the controls
- **fluconazole / voriconazole** carry `precedent: tested-active` and are tagged
  **reference (positive control, not a novel hit)** — the table refuses to promote a known
  inhibitor to a discovery slot, exactly as designed.
- **ketoconazole** lands in C: it recovers its known PAINS/liability profile — the control
  proving the deprioritisation logic bites.
- **R-1479** → C on the same PAINS rule, consistent with the ADMET axis (#75) that already
  flagged its azide/diazo as chemically implausible. The fused artifact re-derives that
  demotion from primitives, not from a hand-copied note.

## The single priority-A: verinurad
> Ranks #2 by N-Fe geometry with convergent poses (2/20 contact the iron), no prior
> fungal-CYP51 measurement (novel, but absence is not proof of untested).

Automatic-A requires **all five**: top-decile geometry rank, convergent poses, not a known
active, clean drug-likeness (≤1 Ro5 violation and 0 structural alerts), and physically
obtainable. verinurad is the only novel molecule clearing every gate. Read the caveat the
other axes attach: verinurad is flagged **coordinator-fragile** (single aromatic N — a
stricter ≥2 rule would drop it, `RESULTS_domain_sensitivity.md` #87), and it is stable but
only *provisional* on the human-selectivity axis. Priority A means "test this first," not
"this is the answer."

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
