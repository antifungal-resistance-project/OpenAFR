# Applicability-domain threshold sensitivity (#87)

Date: 2026-08-26

The shortlist is only valid "within the declared applicability domain (≤45 heavy atoms, ≥1
aromatic N)". This checks whether those two lines are **principled** — set by where the
known actives live and by docking reliability — or arbitrary cutoffs that happen to include
our hits. Fully local, deterministic RDKit reanalysis (`scripts/domain_sensitivity.py`,
`tests/test_domain_sensitivity.py`), no docking, no network. The two thresholds are imported
from `scripts/filter_library.py` (single source of truth), never re-hardcoded.

Data: 15 validated actives (`data/ligands/actives.smi` + `actives_holdout_final.smi`), 399
property-matched decoys, and the 16-molecule focused shortlist (`shortlist_focus.smi`).

## The heavy-atom cap (45) sits in a gap — it is maximally robust

Sweeping the cap (aromatic-N held at ≥1), retained / total:

| cap | validated actives | decoys | shortlist |
|---:|---|---|---|
| 30 | 11/15 | 197/399 | 13/16 |
| 35 | 12/15 | 257/399 | 14/16 |
| **36** | 13/15 | 268/399 | **16/16** |
| 40 | 13/15 | 300/399 | 16/16 |
| **45** *(frozen)* | **13/15** | **301/399** | **16/16** |
| 48 | 13/15 | 336/399 | 16/16 |
| 50 | 14/15 | 370/399 | 16/16 |
| 55 | 15/15 | 399/399 | 16/16 |

The retained *actives* and the retained *shortlist* are **identical for every cap from 36
to 48** — the frozen 45 sits in the middle of a 12-atom dead zone. Concretely: the largest
in-domain validated active is **ketoconazole at 36 heavy atoms**, the smallest out-of-domain
active is **itraconazole at 49**, and the two excluded (itraconazole 49, posaconazole 51)
are exactly the long-tailed azoles the pre-registration
(`PREREGISTRATION_run2.md`) cites as reproducible induced-fit failures on the rigid
receptor. The cap is therefore set at the **docking-reliability boundary**, not tuned to the
shortlist: no active and no shortlist candidate lives anywhere near it.

**No shortlist candidate is size-fragile.** The largest is 36 heavy atoms
(taranabant, ketoconazole) — a **9-atom margin** below the cap. You would have to lower the
cap to 36 before losing any candidate at all, and to ~30 before losing three.

## The aromatic-N floor (≥1) is a scope gate, not a discriminator

Sweeping the aromatic-N floor (cap held at 45):

| aromatic-N ≥ | validated actives | decoys | shortlist |
|---:|---|---|---|
| 0 | 13/15 | 301/399 | 16/16 |
| **1** *(frozen)* | **13/15** | **301/399** | **16/16** |
| 2 | 13/15 | 191/399 | 13/16 |

Going 0 → 1 changes **nothing**: **all 399 decoys already carry ≥1 aromatic N** (they were
property-matched to the azole actives), and **every validated active has ≥2**. So the ≥1
rule does not separate actives from decoys — it is a *mechanistic-scope* gate (the criterion
is a nitrogen-to-iron distance; a molecule with no aromatic N cannot coordinate the heme
iron and is undefined under it), exactly as documented. It is correctly the loosest rule
consistent with the mechanism.

## Threshold-fragile candidates (flagged)

Tightening the floor to ≥2 aromatic N would drop **three** candidates — the ones whose
entire in-domain status rests on a *single* potential iron coordinator:

| candidate | heavy | aromatic N | flag |
|---|---:|---:|---|
| **verinurad** | 25 | 1 | coordinator-fragile |
| **flucloxacillin** | 30 | 1 | coordinator-fragile |
| **taranabant** | 36 | 1 | coordinator-fragile |

This is a genuinely useful flag, not a threshold artifact: for these three the geometry rank
depends on one aromatic N being the actual heme coordinator. For flucloxacillin (that N sits
in an isoxazole) and taranabant (a pyridine N on a distal ring) it is worth confirming the
docked pose coordinates through *that* nitrogen before trusting the rank — the same caution
the pose-convergence (#71) and heme-geometry-audit (#72) axes will sharpen. The other 13
candidates carry ≥2 aromatic N and are robust to this line.

## Verdict

Both cutoffs are principled, and the shortlist is insensitive to them:

- **≤45 heavy atoms** is set at the induced-fit reliability boundary (posaconazole/
  itraconazole fail), sits in a 12-atom gap containing no active and no candidate, and
  leaves every candidate a ≥9-atom margin. No candidate is size-fragile.
- **≥1 aromatic N** is the mechanistic floor, non-discriminating between actives and
  decoys; three candidates sit exactly on it and are flagged coordinator-fragile.

## Reproduce

```
conda activate openafr
python scripts/domain_sensitivity.py            # full report
python scripts/domain_sensitivity.py --json     # machine-readable
python -m pytest tests/test_domain_sensitivity.py -q   # 4 tests
```

Artifact: `work/candidates/domain_sensitivity.json`. Script `scripts/domain_sensitivity.py`;
thresholds from `scripts/filter_library.py`; tests `tests/test_domain_sensitivity.py`.
