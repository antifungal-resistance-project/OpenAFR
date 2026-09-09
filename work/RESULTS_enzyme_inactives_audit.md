# The §4.1 "most useful next contribution" is not buildable: enzyme-confirmed CYP51 inactives do not exist at scale in public data

Date: 2026-09-08
Script: [`scripts/audit_enzyme_inactives.py`](../scripts/audit_enzyme_inactives.py) (reuses the
frozen decoy filters in [`openafr/inactives.py`](../openafr/inactives.py) unchanged).
Source: ChEMBL target **CHEMBL1780** — *C. albicans* (SC5314) lanosterol 14α-demethylase, the
UniProt-**P10613** enzyme behind the docked 5TZ1 structure. Pulled live (816 records, the full
target); the cache is regenerated, not versioned, like every other ChEMBL pull in this repo.

This is an availability **audit**, not a pre-registered gate. It decides no criterion and docks
nothing. It answers one question before any run is designed: *do the data the preprint calls for
even exist?*

## Why this audit, and what it was going to feed

The preprint (`work/PREPRINT_geometry_ceiling.md`) names, twice, the single most useful next
contribution to this benchmark:

> §4.1 — "The only clean exit is **experimentally confirmed inactives** — compounds measured not
> to inhibit … obtaining them for CYP51 is, in our view, the most useful next contribution anyone
> could make to this benchmark."
> §6, Limitation 2 — "CYP51-**enzyme-confirmed** (not merely whole-cell) inactives remain the ideal
> test and the most useful next contribution (§4.1)."

The motivation is concrete. Limitation #2 was powered at N=300 against 279 compounds measured not
to inhibit *whole-cell* C. albicans (look #6/#11, `RESULTS_verified_power.md`): pooled **AUC 0.688
[0.656, 0.720]**, sub-bar, because (a) 58% of that set are *failed azole analogues* — a
chemotype-vs-potency confound (S2 within-warhead AUC 0.622) — and (b) a whole-cell MIC is not a
target assay: a compound can inhibit CYP51 yet fail the cell through permeability/efflux, a false
inactive. The clean exit from both is compounds measured not to inhibit the **enzyme itself**. The
plan was the exact single-variable swap look #6 used — hold every decoy filter byte-identical,
change only the *provenance* of "inactive" from whole-cell-measured to enzyme-measured — then
re-grade and, if it cleared the bar, retire Limitation #2 the way look #10 retired Limitation #1.

**That set cannot be built.** The audit below is why, and it is the deliverable.

## Result: the ideal inactive set is essentially empty

| evidence rule | enzyme-inactive compounds | after the frozen decoy filters |
|---|---:|---:|
| **CLEAN** (only unambiguously-classifiable rows) | 5 / 184 | **0** |
| **GENEROUS** (indefensible upper bound) | 24 / 184 | **6** |

A powered contrast needs ~O(100) property-matched non-inhibitors (cf. the 279 whole-cell set that
still only reached CI lower bound 0.656). The clean rule yields **zero**; even a deliberately
over-loose rule yields **six**, and those six rest on a threshold no benchmark should use.

## Why the data is not there — three findings, each reproducible

1. **The target is tiny and active-dominated.** All of CHEMBL1780 is **816 records across just 184
   distinct compounds**. Under the clean rule: 49 have active evidence, **5** are confirmed
   non-inhibitors, 130 have no usable verdict. The IC50/Ki/Kd values (n=44, all in nM) are *all
   potent* — max 1153 nM, no censored high non-binders — because binding databases publish hits,
   not misses (the same file-drawer effect documented in `METHODS_file_drawer.md`).

2. **The only non-inhibition signal lives in %-assays, and most of it is directionally ambiguous.**
   650 of 816 records are in `%` units. The `Inhibition`-type %-rows (n=101) are directly
   interpretable (low % = non-inhibitor), but the bulk — 496 `Activity`-type %-rows — are
   **per-sterol GC-MS composition** readouts ("… assessed as ergosterol/eburicol/lanosterol
   composition at 8 µg/ml"). Their direction is opposite per sterol: a CYP51 inhibitor *lowers*
   ergosterol but *raises* the eburicol/lanosterol/obtusifoliol substrate, so a low value is
   "inactive" or "active" depending only on which sterol the row measures. Classifying them by
   value alone would fabricate labels; the clean rule refuses them. The 5→0 result comes from the
   honest subset; the generous 24→6 comes from ignoring this and it still fails.

3. **Even this "enzyme" target is 81% cell-based.** 660/816 assay descriptions are cellular
   ("Inhibition of CYP51 **in** Candida albicans SC5314 …", "ergosterol biosynthesis in cells",
   "… in membrane"), not cell-free purified-enzyme inhibition. So CHEMBL1780 does not cleanly
   escape the permeability/efflux confound that motivated moving off whole-cell data in the first
   place — the truly cell-free enzyme measurements are a small minority.

**Orthologues do not rescue it.** ChEMBL has 17 sterol-14α-demethylase targets (*T. cruzi* 861,
human 155, *A. fumigatus* 5, *C. glabrata* 3, …), but a compound that fails to inhibit *T. cruzi*
or human CYP51 says nothing about *C. albicans* CYP51 and is untested against the docked enzyme —
counting it inactive would be a wrong cross-species label, exactly the move
`RESULTS_precedent_bindingdb.md` flags as a deliberate cross-species claim. There is no valid
same-target expansion.

## What this establishes

1. **§4.1's "most useful next contribution" is unperformable on public data, not merely
   unperformed.** The exhaustive, reproducible count is 0 (clean) / ≤6 (generous) property-matched
   enzyme-confirmed inactives for the docked target. This converts the preprint's hand-wave ("we do
   not have those") into a hard, audited number.

2. **The whole-cell set (279) is not a stopgap for a better set that exists — it is the best set
   that exists.** Limitation #2's powered 0.688 is therefore the *ceiling of what public data can
   say*, not an artifact of settling for whole-cell data. That strengthens, not weakens, the
   paper's honesty: the caveat cannot be discharged by anyone, absent new wet-lab measurement.

3. **The confound is structural to the literature.** Failing enzyme measurements are under-published
   (file-drawer), CYP51 assays are overwhelmingly run as cellular sterol-profile readouts, and the
   only compounds anyone tests are antifungal med-chem programs (hence the 58% azole enrichment of
   the whole-cell set). No curation move fixes a set that was never deposited.

## Consequence for the preprint (§4.1, §6 Limitation 2)

Update §4.1 and Limitation 2 from "the most useful next contribution *anyone could make*" to an
audited statement: *we searched the docked target (CHEMBL1780, 184 compounds) exhaustively;
property-matched enzyme-confirmed inactives number 0 under a clean rule and ≤6 under an
indefensibly loose one, 81% of the target's records are cell-based, and no valid same-target
expansion exists — so the enzyme-confirmed-inactive test that would settle Limitation 2 is not
achievable on public data and requires new experimental non-inhibitor measurement.* This closes the
open action item rather than leaving it as an aspiration.

## Reproduce

```bash
conda activate openafr
python scripts/audit_enzyme_inactives.py          # live ChEMBL pull of CHEMBL1780, prints the funnel
# --refetch to force a fresh pull; --cache DIR to relocate the (regenerated) cache
```
