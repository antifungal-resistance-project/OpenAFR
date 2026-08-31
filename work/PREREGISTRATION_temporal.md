# Pre-registration — temporal (prospective) holdout validation (issue #82)

**Written 2026-08-31, BEFORE building or docking a single molecule of this set.** Everything
below is fixed in advance. Read [PREREGISTRATION_external.md](PREREGISTRATION_external.md) and
[RESULTS_external.md](RESULTS_external.md) first — this is the **stronger follow-up** that
pre-registration named, and it reuses that protocol, receptor, grading math and three-outcome
rule unchanged. Only the *provenance* of the holdout changes: from "never-scored subset of the
same curated population" to "deposited strictly after the tuning corpus."

## The specific objection this closes

The disjoint-population holdout (#108) answered the overfitting objection with an active/decoy
set no gate had ever scored (mode-C AUC 0.715, p < 0.0001, CI 0.630–0.795, SUPPORTED-BUT-WIDE).
It left exactly one honest residual on the record:

> Your "never-scored" actives are still a subset of the **same** curated ChEMBL *C. albicans*
> azole population the criterion was developed against. That is a clean never-*touched* holdout,
> but not a **later time slice** — chemistry that post-dates the criterion. A genuinely
> prospective test needs molecules that did not exist when you set the 3.0 Å cutoff.

This pre-registration answers exactly that, and only that. It does **not** re-open run 2, the
external holdout, or the criterion; it grades the frozen protocol, once, on a temporally
disjoint holdout.

## The temporal split (fixed)

Every molecule in this holdout has its **earliest** measured antifungal / CYP51 activity
reported in a ChEMBL document dated **strictly after the tuning corpus** — so it could not have
informed the criterion, its cutoff, or any decoy-selection choice, even in principle.

- **Cutoff — `CUTOFF_YEAR = 2026`.** The tuning corpus (`verified_actives`, `verified_inactives`,
  every decoy set, both potency sets, and the #108 external holdout) was pulled from the ChEMBL
  release current in **2026-08**, whose documents can carry year 2026. A molecule is temporal
  iff its earliest antifungal / CYP51 document year is **> 2026** (i.e. ≥ 2027). Strict `>`
  (not `≥`) guarantees zero overlap with the corpus's latest possible year. The predicate is
  `scripts.make_temporal_holdout.is_temporal`; the per-compound year is `earliest_document_year`
  (min over its records that carry a year; an unknown year is **not** temporal — absence of
  proof is not proof).

## The holdout (fixed; frozen by sha256 before docking, exactly as #108)

Built by `scripts/make_temporal_holdout.py` under the same provenance discipline as every other
set (`data/ligands/PROVENANCE.md`).

- **Actives — `data/ligands/temporal_actives.smi`.** Measured-active azole CYP51 antifungals
  (`openafr.inactives.compound_verdict == has-active-evidence`; azole-bearing by the gate's
  `AZOLE_SMARTS`) whose earliest antifungal / CYP51 document year > 2026; in the frozen domain
  (≤ 45 heavy, ≥ 1 aromatic N); **novel** (max Morgan r=2 Tanimoto < 0.70 vs every prior active
  + the 5TZ1 co-crystal, the same `NOVELTY_MAX` the verified-active build uses); and disjoint by
  ChEMBL id from **every** set the project has used, the #108 external sets included.
  **sha256 `__ACTIVES_SHA__` (filled by the release-gated build; NOT-YET-AVAILABLE today).**
- **Decoys — `data/ligands/temporal_decoys.smi`.** Property-matched presumed-inactives, built by
  the frozen rule (`openafr.inactives.select_inactives`: MW ± 25, cLogP ± 1.5, rotB ± 2, HBD ± 1,
  HBA ± 2; ≥ 1 aromatic N; Morgan r=2 Tanimoto < 0.35 vs every temporal active; ≤ 45 heavy),
  themselves drawn **only from post-cutoff documents** and disjoint by id from every used set.
  **sha256 `__DECOYS_SHA__` (filled by the release-gated build).**

Any active with no rankable pose is ranked **LAST**, never dropped — the same anti-inflation
rule as every prior gate.

## Release gate (why the numeric run is deferred)

This is **release-gated**, the direct analogue of the metal-aware baselines' cloud gate (#83):
the temporal holdout materializes only once a ChEMBL release **newer than the 2026-08 tuning
release** accrues antifungal / CYP51 documents dated > 2026. **Today the current CYP51/azole
universe is exhausted** — the same wall #82's #108 comment named ("needs a newer ChEMBL
release") — so `make_temporal_holdout.py` writes nothing and reports NOT-YET-AVAILABLE, and
`validate_gate_temporal.py` reports `not-yet-run (release-gated)`. The **rule, cutoff,
disjointness, domain and interpretation are all frozen now** and unit-tested offline
(`tests/test_temporal_gate.py`); only the fetch + dock + grade wait on the newer release. This
PR ships the reproducible scaffold and the frozen pre-commitment, not a number.

## Protocol (fixed; identical to run 2 / #108)

- Receptor: rigid 5TZ1 chain A + heme, `work/receptor_A.pdb`
  (sha256 `a1a8410868d1b9ee2835117cf6327ba4721ca3830a4d5689b7e9ea0ee1eee22f`); docking into
  `work/receptor.pdbqt` prepared by `scripts/prep_receptor.py`.
- Box + search: from the hash-frozen `protocol.yaml`
  (sha256 `6eca16b1be2e5e33ae42d4b1b7a63a0e9b891a5dd6cbfb19642dd3e9a149aba4`) — box 26 Å @
  70.61/66.28/4.18, exhaustiveness 32, num_modes 20, seed 42. No per-molecule tuning.
- Ligand prep: `obabel … -p 7.4`, exactly as `scripts/screen.sh`.

## Primary criterion, endpoints and pass bar (fixed)

**Primary (gates): mode C** — rank by the minimum azole-nitrogen-to-heme-iron distance over a
molecule's poses (lower = better), the run-2 PRIMARY criterion, reproduced unchanged.
**Secondary (report only): mode A** (best Vina affinity) — the run-2 worse-than-random baseline.

**Endpoints: AUC and EF@5%** (issue #82). EF@1% is excluded as an endpoint
([RESULTS_reliability.md](RESULTS_reliability.md): not robust under seed sampling) though printed
for continuity. Reported alongside: a one-sided permutation test (20,000 shuffles, seed 42) and a
percentile bootstrap 95% CI over the actives (10,000 resamples, seed 42), from the same
`openafr.scoring` functions every prior gate uses.

**Pass bar: mode-C AUC ≥ 0.70** — the frozen protocol's `gate.min_auc`.

## Pre-committed interpretation (three outcomes, fixed before seeing numbers)

Reused verbatim from #108 (`scripts.validate_gate_external.classify_outcome`), read on mode-C
**AUC** at the ≥ 0.70 bar with the bootstrap CI qualifying it:

- **PASS — AUC ≥ 0.70 and bootstrap 95% CI lower bound ≥ 0.70.** The criterion generalizes to
  chemistry that did not exist when it was built, with the whole interval above the bar.
  Overfitting to the development *era* is excluded — the strongest available outcome.
- **SUPPORTED-BUT-WIDE — AUC ≥ 0.70 but CI lower bound < 0.70.** Point estimate passes and the
  permutation test rules out chance, but the interval still admits a true AUC below the bar.
  Prospective generalization is supported, not proven at the bar — reported exactly as such.
- **FAIL — AUC < 0.70.** The ranking result does **not** transfer to a post-dating set;
  temporal overfitting cannot be excluded. Reported as a genuine negative. Per pre-registration
  the criterion is **not** re-tuned to recover a pass, and no alternative ranking is substituted
  and reported as a win.

## Honest scope (fixed in advance)

1. **Temporal, but same organism/target world.** The holdout post-dates the corpus but is still
   *C. albicans* azole CYP51 chemistry; it tests generalization across **time**, not across a
   different pathogen or a non-azole class (those are separate, already-scoped questions).
2. **Decoys are presumed inactive, not measured** — any real binder among them depresses the
   measured AUC; any decoy-selection bias inflates it. Same limitation as run 2 / #108.
3. **One rigid *C. albicans* receptor**; rigid docking cannot model induced fit; domain ≤ 45
   heavy atoms. Same blind spots as run 2.
4. **This validates a ranking criterion, not a drug.**
5. **Sample size is release-determined.** The temporal `n` is whatever the newer release yields;
   if it is small the CI will be wide, and the SUPPORTED-BUT-WIDE outcome (not a forced PASS) is
   the honest read — identical to the n=50 external holdout's stance.

## Reproduce (on the newer release; release-gated today)

    python scripts/make_temporal_holdout.py status              # report the release gate
    python scripts/make_temporal_holdout.py actives             # fetches ChEMBL (release-gated)
    python scripts/make_temporal_holdout.py decoys              # fetches ChEMBL (release-gated)
    RECEPTOR=work/receptor.pdbqt scripts/screen.sh work/temporal_ligands work/screen_temporal
    python scripts/validate_gate_temporal.py work/screen_temporal work/receptor_A.pdb

Until a newer release exists the grader prints `not-yet-run (release-gated)` and exits
incomplete — the honest current state, mirroring the metal-aware baseline harness (#83).
