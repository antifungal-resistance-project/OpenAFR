# Human-CYP off-target liability of the candidate shortlist (#74)

Date: 2026-08-29

The selectivity axis the shortlist was missing. `human_CYP51_selectivity` (#73) already asks
whether a candidate binds the *human sterol-14α-demethylase* it most resembles. This axis asks
the other, louder off-target question: does the candidate inhibit the human **drug-metabolizing**
cytochrome P450s — CYP3A4, 2C9, 2D6, 2C19, 1A2 — the enzymes whose inhibition causes the
drug–drug interactions and hepatotoxicity the marketed azoles are infamous for? A geometry tool
that only checks the on-target is silent on this, and it is silent on the *most predictable*
liability an azole-shaped molecule carries.

The uncomfortable part is that the liability is the **same feature we select for**. Azole
antifungals inhibit human CYPs by *type-II* binding: an sp2 ring-nitrogen lone pair ligating the
ferric heme iron as the sixth axial ligand — atom-for-atom the coordination our geometry
criterion rewards on CYP51. The heme does not know which P450 it sits in. So the warhead that
shortlists a molecule is, unchanged, the warhead that makes it a CYP3A4 inhibitor. This axis
makes that visible per candidate (`openafr/cyp_offtarget.py`, `scripts/cyp_offtarget_screen.py`),
a fully **local, deterministic RDKit computation** — no network, no docking.

> **A flag, not an IC50.** Docking these five large, famously flexible pockets rigidly would
> manufacture false precision (the issue says as much); a trained QSAR would need a model file
> and a provenance we do not have offline. So this axis stays at the honest tier — a mechanistic
> + substrate-pharmacophore *liability level* (HIGH / MODERATE / LOW per isoform), the same claim
> class as the hERG risk-factor hint in `openafr.admet`. Informational, never a gate. For a
> candidate that is already an approved drug the authoritative source is the **drug label**.

## Two signals, deliberately separated

1. **Type-II heme ligation (mechanistic, pan-CYP).** An accessible azole warhead —
   imidazole / 1,2,4- or 1,2,3-triazole / tetrazole carrying a pyridine-type (two-coordinate,
   unprotonated, no-H) ring nitrogen — is the strong, documented antifungal-class type-II ligand
   and is flagged **HIGH against all five isoforms at once**. A weaker ligator — a pyrazole
   (1,2-diazole) or a bare non-azole azine nitrogen (pyridine / pyrimidine / …) — is a real but
   softer type-II center and is flagged **MODERATE pan-CYP**.
2. **Per-isoform substrate pharmacophore (coarse, literature-grounded hints).** Each human CYP
   has a recurring substrate shape; a match flags **MODERATE** on that isoform only:
   3A4 large+lipophilic (MW≥400, cLogP≥3) · 2D6 basic N + aromatic · 2C9 acidic group +
   lipophilic · 2C19 moderately lipophilic aromatic amide/basic · 1A2 small planar poly-aromatic.

Per isoform the reported liability is the **stronger** of the two, and a `basis` list records
every reason that fired.

## The reference drugs validate the axis in-run

The shortlist carries the real azoles as controls, and they profile exactly as they must:

- **ketoconazole, fluconazole, voriconazole, letrozole, ravuconazole, flutrimazole — all HIGH
  type-II across the whole panel.** Ketoconazole is *the* textbook potent CYP3A4 inhibitor,
  black-boxed for drug interactions and hepatotoxicity; the axis recovers exactly that from the
  imidazole warhead alone, with no label lookup. A CYP off-target axis that did **not** flag
  ketoconazole would be broken.

Controls are **exempt from the demotion** on the scorecard: their HIGH type-II is the point of a
control (it validates the axis), not a strike against a discovery.

## What it caught among the novel hypotheses

Of the 10 novel candidates, **2 carry the strong azole warhead and inherit the full HIGH
pan-CYP liability**:

| candidate | class | warhead | verdict |
|---|---|---|---|
| **L-838417** | GABA-A modulator | 1,2,4-triazole | **HIGH type-II ×5** — the same CYP3A4 DDI liability as the marketed triazoles; undercuts any "cleaner than existing azoles" thesis |
| **olprinone** | PDE3 inhibitor | imidazole | **HIGH type-II ×5** — imidazole warhead; expect azole-class CYP inhibition |

Both gain a `human-CYP type-II inhibition liability (azole warhead; CYP3A4 DDI risk)` concern on
the consolidated scorecard. This is the axis doing its job: two candidates that looked clean on
most chemistry axes carry the one liability that is most predictable from their shape.

The remaining 8 novel candidates are **MODERATE only** — a bare azine nitrogen (verinurad,
vatalanib, taranabant, LDN-27219, R-1479, nolatrexed), a weaker pyrazole (lersivirine), or a
pure substrate-pharmacophore match (flucloxacillin: 2C9-acid + 2C19). MODERATE hits are reported
on the card but do **not** demote — nearly every heme-seeking scaffold carries one, so a strike
there would be uninformative. The honest reading is comparative: **most of the novel shortlist
carries at least a weak human-CYP liability, and two carry the full azole liability** — a direct
consequence of screening for a heme-ligating nitrogen.

## What this is not

- **Not a docking result and not a QSAR IC50.** It is a mechanistic + pharmacophore liability
  *flag*. The per-isoform pharmacophore rules are coarse 30-year-old-heuristic-tier hints, not
  predictions; the type-II flag fires on *any* accessible azole/azine nitrogen regardless of the
  steric context that would modulate a real Kᵢ (an over-inclusive hint by design).
- **Not the human-CYP51 selectivity axis (#73).** That is the *sterol-demethylase* off-target
  (the enzyme we most resemble, docked); this is the *drug-metabolizing* panel (a flag). They
  answer different questions and are separate columns.
- **No statistical pre-registration.** Like the ADMET axis (#75), this is a fixed rule set over
  the input SMILES — a pure deterministic function with no graded pass/fail bar and nothing to
  p-hack, so the freeze discipline that protects the docking gates does not apply. The rules are
  fixed in `openafr/cyp_offtarget.py` and pinned by `tests/test_cyp_offtarget.py`; the
  scorecard re-derivation is gated by `scripts/repro.py`.

This layer is one column in the per-candidate scorecard (#88), not a verdict on its own. With it
measured, the only axis still `not-yet-checked` for the shortlist is the remaining *C. auris*
mutant panel (K143R/F126L, #77).

## Reproduce

```
conda activate openafr
python scripts/cyp_offtarget_screen.py work/candidates/shortlist_focus.smi          # table
python scripts/cyp_offtarget_screen.py work/candidates/shortlist_focus.smi --tsv    # flat table
python scripts/cyp_offtarget_screen.py work/candidates/shortlist_focus.smi --json   # full profiles
python work/candidates/build_scorecard.py                                           # fold into scorecard
python -m pytest tests/test_cyp_offtarget.py tests/test_scorecard.py -q             # 11 + tests
```

Artifacts: `work/candidates/shortlist_cyp.json`, `work/candidates/shortlist_cyp.tsv`,
folded into `work/candidates/shortlist_scorecard.json`.
Module `openafr/cyp_offtarget.py`; CLI `scripts/cyp_offtarget_screen.py`;
tests `tests/test_cyp_offtarget.py`.
