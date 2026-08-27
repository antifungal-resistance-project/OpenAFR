# Candidate shortlist — the repurposing hits, filtered by assay precedent

Date: 2026-08-25

The first time the project turns a *ranked* library into an actual **candidate
shortlist** — the deliverable the whole method exists to produce. It takes the
2026-08-05 Drug Repurposing Hub screen (`RESULTS_repurposing.md`, 2776 drugs docked
once under the frozen protocol) and runs the new **assay-precedent layer**
(`openafr/precedent.py`, `scripts/check_precedent.py`, PR #65) over the **top 100**,
asking the one question geometry cannot: *has a wet lab already measured this molecule,
and what did it say?* Two probes per candidate — an exact-match lookup against ChEMBL +
PubChem, and a near-neighbour check against the 279 verified-inactive azole analogues.

> **Every rank here is a hypothesis for a wet lab, never a hit.** Nothing below has been
> tested against any fungus *by us*. The list is meaningful only because the validation
> gate passed first, only within the declared applicability domain (≤45 heavy atoms,
> ≥1 aromatic N), and only as a triage over one rigid *C. albicans* 5TZ1 conformer.

## The precedent layer validated itself in-run

Of the top 100, the exact-match enzyme lookup flagged exactly **four** compounds — and
**all four are the known azole controls**, nothing else:

| rank | compound | enzyme verdict | ChEMBL |
|---:|---|---|---|
| 9  | ravuconazole | record-inconclusive | CHEMBL294029 |
| 23 | letrozole    | record-inconclusive | CHEMBL1444 |
| 31 | fluconazole  | **tested-active**   | CHEMBL106 |
| 57 | voriconazole | **tested-active**   | CHEMBL638 |

A precedent filter that did *not* flag fluconazole and voriconazole would be broken;
one that flagged a *novel* repurposing hit as a secret known inhibitor would be crying
wolf. It did neither. This is the live proof the lookup works — not just that it ran.

## What the filter caught among the novel hypotheses

**Enzyme-level: 96/100 `no-precedent`.** No deposited fungal-CYP51 enzyme assay for the
repurposing hits. Read this honestly, as the code insists: `no-precedent` is **not
"untested"** — a failing measurement is often never deposited (publication bias).
Absence is not evidence; it only means no lost experiment is *on record* at the enzyme.

**Whole-cell (`phenotypic`) catch — the file drawer that fired.** Seven novel
hypotheses were already measured **inactive in a cell-based assay**:

| rank | drug | usual class | N–Fe (Å) |
|---:|---|---|---:|
| 12  | **epacadostat** | IDO1 inhibitor (oncology) | 2.674 |
| 18  | **GDC-0879**    | B-Raf inhibitor | 2.705 |
| 47  | **TG-02**       | CDK/FLT3 inhibitor | 2.781 |
| 54  | **H-89**        | PKA inhibitor | 2.788 |
| 56  | **olmesartan**  | angiotensin-II receptor blocker | 2.797 |
| 86  | **A-769662**    | AMPK activator | 2.831 |
| 100 | **lesinurad**   | URAT1 inhibitor (gout) | 2.844 |

Caveat the code states loudly: a whole-cell MIC folds in permeability and efflux, so an
inactive cell assay does **not** prove the CYP51 enzyme was untouched — and it is not
organism-filtered. Treat these seven as "read the underlying assay before spending bench
time," a demotion, not a disproof. The `phenotypic: tested-active` rows (flucloxacillin,
tazobactam, cytarabine, …) are active in *some* cell assay somewhere — not a fungal
verdict, and for the antibiotics almost certainly antibacterial; neither confirmation
nor exclusion here.

**Near-neighbour check: 0/100 reached the Tanimoto ≥ 0.85 bar.** No top-100 hit is an
all-but-identical twin of any of the 279 verified-inactive azole analogues. This is an
expected, clean *negative*, not a null result: the verified inactives are azole-scaffold
analogues, and every hit already passed the <0.70-to-azole-panel novelty filter, so
≥0.85 similarity to them is near-impossible by construction. The check is orthogonally
clean here — it would bite on a library drawn closer to the azole chemotype, and is worth
keeping wired in for exactly that case.

## The surviving shortlist (novel, clean of precedent)

Controls removed, the seven cell-inactives demoted, top novel hypotheses in geometry
rank order — the molecules whose nitrogen comes as close to the heme iron as a real azole
does, with no deposited assay (exact or near) saying they've already failed:

| rank | drug | usual target / class | N–Fe (Å) | precedent |
|---:|---|---|---:|---|
| 2  | verinurad | URAT1 inhibitor (gout) | 2.417 | clean (no-precedent) |
| 3  | vatalanib | VEGFR tyrosine-kinase inhibitor | 2.479 | clean |
| 4  | flucloxacillin | β-lactam antibiotic | 2.558 | cell-active (antibacterial) |
| 5  | taranabant | CB1 inverse agonist | 2.586 | clean |
| 6  | lersivirine | HIV NNRTI | 2.599 | clean |
| 7  | LDN-27219 | transglutaminase-2 inhibitor | 2.603 | clean |
| 8  | R-1479 (balapiravir) | HCV polymerase inhibitor | 2.621 | clean |
| 10 | L-838417 | GABA-A partial agonist | 2.667 | clean |
| 13 | olprinone | PDE3 inhibitor (inotrope) | 2.680 | clean |
| 14 | nolatrexed | thymidylate synthase inhibitor | 2.692 | clean |

## What this run does and does not establish

- **Does:** produce a filtered, evidence-backed candidate shortlist — the method's
  intended output — remove seven molecules a wet lab has already found cell-inactive, and
  prove the precedent layer on the fluconazole/voriconazole controls while showing it does
  not over-flag novel chemotypes.
- **Does not:** claim any molecule is a hit, or that `no-precedent`/`near 0.85 = clean`
  means untested. The five limitations of the repurposing screen still bind (single rigid
  conformer, per-molecule variance — clotrimazole ranked #2774 — SMILES as-supplied,
  *C. albicans* stand-in receptor, coordination necessary not sufficient).

## Follow-up confidence layers (2026-08-26)

Two of the shortlist's biggest caveats have since been attacked directly, on the survivors above:

- **Seed stability** ([RESULTS_shortlist_consensus.md](RESULTS_shortlist_consensus.md), #68):
  re-docked under 5 seeds, 8 of the 10 novel survivors coordinate the iron reproducibly.
  **LDN-27219 (#7) and olprinone (#13) are demoted as single-dock artifacts** — LDN-27219
  coordinates in only 2 of 5 searches (three seeds put it at ~5.5 Å). Its rank was search luck.
- **Human-CYP51 selectivity** ([RESULTS_human_selectivity.md](RESULTS_human_selectivity.md),
  #73, **validated** 2026-08-26): vatalanib, flucloxacillin, taranabant and lersivirine reach the
  fungal iron but never the human iron (margins +1.8–2.7 Å) — a clean structural selectivity
  window; nolatrexed hits human CYP51 as readily as fungal (flagged). Originally provisional
  because the pre-registered ketoconazole control did not converge; a native-pose-recovery control
  ([RESULTS_human_control.md](RESULTS_human_control.md)) has since validated the human receptor
  (crystal ketoconazole re-coordinates at N-Fe 2.651 Å), isolating that failure as global-search
  power on a floppy ligand, not a receptor fault.

**Refined tiering:** strongest = vatalanib, flucloxacillin, taranabant, lersivirine (stable +
selective); stable/provisional = verinurad, L-838417, R-1479; demoted = LDN-27219, olprinone
(unstable), nolatrexed (hits human).

- **ADMET / structural-alert profile** ([RESULTS_admet.md](RESULTS_admet.md), #75): the
  chemistry-realism axis. Informational, never a gate. The reference azoles validate it
  (fluconazole/voriconazole clean; ketoconazole recovers its known Ro5/solubility/hERG
  liabilities). Among the novel hits it flags **R-1479** (azide/diazo/PAINS — chemically
  implausible), **taranabant** (2× Ro5, insoluble), and **LDN-27219** (reactive hydrazide);
  verinurad, lersivirine, ravuconazole, L-838417 profile clean.
- **Synthesizability / availability** ([RESULTS_synth.md](RESULTS_synth.md), #86): every
  candidate is buyable (`catalogued-drug`, from the Drug Repurposing Hub) and lands in the
  same easy/moderate SAscore band as the marketed azoles — none is a synthetic outlier
  (highest R-1479 at 4.46, still < the "hard" 6.0 line).
- **Domain-threshold sensitivity** ([RESULTS_domain_sensitivity.md](RESULTS_domain_sensitivity.md),
  #87): the ≤45-heavy-atom cap sits in a 12-atom gap (no active or candidate near it; largest
  candidate 36) and the ≥1-aromatic-N floor is a mechanistic scope gate, not a discriminator.
  No candidate is size-fragile; **verinurad, flucloxacillin and taranabant** are flagged
  *coordinator-fragile* (a single aromatic N — a stricter ≥2 rule would drop them).

**Consolidated scorecard (2026-08-27, #88).** All axes above are now rolled onto one row
per candidate — the artifact you actually reason over — in
[RESULTS_scorecard.md](RESULTS_scorecard.md) (`work/candidates/shortlist_scorecard.json`,
builder `work/candidates/build_scorecard.py`). Seen together, **lersivirine and L-838417**
are the only two novel hits with zero open concern on any *measured* axis; every card also
prints the axes still `not-yet-checked` (candidate-path ensemble #70, RMSD clustering #71,
auris/mutant retention #69/#77, protomer/stereo #84/#85, broad CYP panel #74) so no molecule
looks fully vetted.

## Reproduce

    conda run -n openafr python work/candidates/build_shortlist.py 100
    conda run -n openafr python scripts/check_precedent.py \
        work/candidates/top100.smi --near data/ligands/verified_inactives.smi \
        -o work/candidates/top100_precedent.tsv
    conda run -n openafr python work/candidates/merge_view.py \
        work/candidates/top100_precedent.tsv

Inputs persist from the 2026-08-05 screen (`work/repurposing/ranked.tsv`,
`library_in_domain.smi`); the docked poses do not (gitignored), so this shortlist is
built from the persisted ranking + SMILES, not re-docked. Precedent calls hit ChEMBL
(www.ebi.ac.uk) and PubChem live on 2026-08-25.
