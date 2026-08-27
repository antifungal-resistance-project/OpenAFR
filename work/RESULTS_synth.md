# Synthesizability / availability of the candidate shortlist (#86)

Date: 2026-08-26

The practicality axis of the hand-off: a ranked hypothesis is only useful if a wet lab can
put the molecule in a tube. `openafr/synth.py` + `scripts/synth_filter.py` add two signals
of very different strength — a **synthetic-accessibility score** (local, deterministic) and
an **availability** annotation.

> **Informational, never a gate.** SAscore is a heuristic, not a route; a low score is not
> a guarantee and a high one is not a veto. Availability here is honest about its own
> weakness (below).

## Synthesizability — SAscore (local, no network)

Ertl & Schuffenhauer's synthetic-accessibility score (RDKit-Contrib), 1 (trivial) → 10
(very hard), from fragment contributions plus a complexity penalty. For the shortlist:

| candidate | SAscore | band |
|---|---:|---|
| vatalanib | 2.12 | easy |
| LDN-27219 | 2.22 | easy |
| letrozole | 2.38 | easy |
| verinurad | 2.61 | easy |
| nolatrexed | 2.68 | easy |
| olprinone | 2.72 | easy |
| flutrimazole | 2.71 | easy |
| fluconazole *(ref)* | 2.77 | easy |
| lersivirine | 2.85 | easy |
| L-838417 | 2.90 | easy |
| taranabant | 3.43 | easy |
| ketoconazole *(ref)* | 3.44 | easy |
| ravuconazole | 3.48 | easy |
| voriconazole *(ref)* | 3.59 | moderate |
| flucloxacillin | 3.67 | moderate |
| **R-1479** | **4.46** | **moderate** |

Read alongside the reference azoles: every novel hit sits in the same easy/moderate band as
the marketed drugs (fluconazole 2.77, voriconazole 3.59, ketoconazole 3.44) — none is a
synthetic outlier. The highest is **R-1479** (4.46), the azido-nucleoside — consistent with
its ADMET flags (#75): a stereochemically dense sugar with a reactive azide is both harder
to make and chemically fragile. Still "moderate", not "hard" (SAscore < 6), so synthesis is
not the reason to deprioritise it — its reactive-alert profile is.

## Availability — honest about its limits

Every shortlist entry is annotated **`catalogued-drug`**: the shortlist is drawn from the
Broad **Drug Repurposing Hub**, a curated library of approved / clinical / catalogued
compounds that are by construction vendor-listed and buyable. So for *this* shortlist the
sourcing question is already answered — the check matters most for any future non-drug
library, and is wired now so it is ready then.

The module also ships an **optional** PubChem presence probe (`--online`): a found CID is a
*weak* proxy for "buyable" — presence in PubChem is not a vendor quote — so it is off by
default and clearly labelled. A true buyable check needs a vendor-catalogue API
(ZINC/eMolecules), which is out of scope here; the probe parsing is tested offline with
canned JSON, behind an injectable opener, exactly like the precedent layer.

## Reproduce

```
conda activate openafr
python scripts/synth_filter.py work/candidates/shortlist_focus.smi --catalogued        # table
python scripts/synth_filter.py work/candidates/shortlist_focus.smi --catalogued --tsv  # flat
python scripts/synth_filter.py lib.smi --online --json                                 # + PubChem probe
python -m pytest tests/test_synth.py -q                                                # 7 tests
```

Artifacts: `work/candidates/shortlist_synth.tsv`, `work/candidates/shortlist_synth.json`.
Module `openafr/synth.py`; CLI `scripts/synth_filter.py`; tests `tests/test_synth.py`. One
more column for the eventual per-candidate scorecard (#88).
