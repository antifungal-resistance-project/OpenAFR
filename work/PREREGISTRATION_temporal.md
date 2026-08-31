# Pre-registration — temporal (publication-date) holdout validation (issue #82)

**Written 2026-08-31, BEFORE docking a single molecule of this set.** The active/decoy sets are
built and frozen by sha256 below; no pose has been generated. Everything below is fixed in
advance. Read [PREREGISTRATION_external.md](PREREGISTRATION_external.md) and
[RESULTS_external.md](RESULTS_external.md) first — this is the **stronger follow-up** that
pre-registration named, and it reuses that protocol, receptor, grading math and three-outcome
rule unchanged. Only the *provenance* of the holdout changes: from "never-scored random subset
of the curated population" to "the newest-literature azoles the criterion never saw."

## The specific objection this closes

The disjoint-population holdout (#108) answered the overfitting objection with an active/decoy
set no gate had ever scored (mode-C AUC 0.715, p < 0.0001, CI 0.630–0.795, SUPPORTED-BUT-WIDE).
It left exactly one honest residual on the record:

> Your "never-scored" actives are still a random subset of the same curated population. That is
> a clean never-*touched* holdout, but not a **later time slice** — chemistry that post-dates
> the criterion. A stronger test holds out the *newest* molecules and asks whether a criterion
> built on old drugs generalizes forward in time.

This pre-registration answers exactly that, and only that. The criterion — the iron-approach
filter (N–Fe < 3.0 Å), its 3.0 Å cutoff, and the frozen run-2 protocol — was **defined and
validated on the classic approved azoles** (fluconazole, voriconazole, itraconazole,
posaconazole, ketoconazole, clotrimazole, miconazole, isavuconazole; held-out through
efinaconazole/luliconazole — antifungal drugs approved into the 2010s). This holdout is azoles
whose **earliest** antifungal / CYP51 literature appearance is recent, that no gate has docked.

## Scope, stated honestly (fixed)

This is a **publication-year split of the current ChEMBL release** (ChEMBL_37, 2026-05-01 — the
release this project's corpus was pulled from; confirmed the latest via the ChEMBL `/status`
API on 2026-08-31). It is **not** a *future*-prospective set post-dating the release: ChEMBL_37
is the latest release, so a set drawn from depositions that post-date it is **calendar-gated
for everyone** until ChEMBL_38 and cannot be built today by any means (a compute host does not
change this — the data does not exist yet). What this holdout *does* deliver now, executed
end-to-end, is a genuine **temporal generalization** test: does a criterion built on old
approved drugs rank the newest research azoles in the literature above property-matched decoys?
When ChEMBL_38 ships, re-running `scripts/make_temporal_holdout.py years` and re-freezing lifts
this to the fully future-prospective variant under the identical protocol.

## The holdout (fixed; frozen by sha256 before docking)

Built by `scripts/make_temporal_holdout.py` under the same provenance discipline as every set
(`data/ligands/PROVENANCE.md`). Per-molecule earliest antifungal year from the pinned map
`data/ligands/chembl_document_years.json`
(sha256 `6cda545ebb1ec9ee0e3f4ac8641712ec152738bda8467ee44d0f01651b64c3b1`), built by the `years`
command from the whole-cell *C. albicans* + CYP51 (CHEMBL1780) activity documents.

- **Actives — `data/ligands/temporal_actives.smi`, n = 50.** A deterministic seed-**82** sample
  of the `verified_actives` azoles whose earliest antifungal / CYP51 document year **> 2023**
  (i.e. ≥ 2024), that are **not** in the seed-42 docked sample (`verified_actives_sample.smi`)
  and **not** in the #108 external holdout (`external_actives.smi`) — so genuinely never docked
  by any gate. 320 such never-docked candidates exist; 50 are sampled (matching #108's n for a
  like-for-like contrast). Achieved: earliest year 2024–2025 (41× 2024, 9× 2025), heavy 15–41.
  Applicability domain re-asserted (≤ 45 heavy, ≥ 1 aromatic N).
  **sha256 `6ab47fcd39a044f2dce736cbcb813148cf5f68a80a1b872456a29d72a976a102`**
- **Decoys — `data/ligands/temporal_decoys.smi`, n = 262.** Property-matched presumed-inactives
  built by the frozen rule (`openafr.inactives.select_inactives`: MW ± 25, cLogP ± 1.5, rotB ± 2,
  HBD ± 1, HBA ± 2; ≥ 1 aromatic N; Morgan r=2 Tanimoto < 0.35 vs every temporal active; ≤ 45
  heavy), fetched fresh from ChEMBL and **disjoint by ChEMBL id** from every set the project has
  used (all actives, both decoy sets, verified inactives, both potency sets, AND the #108
  external sets) plus the temporal actives. Decoys are **not** date-constrained: the temporal
  variable is the *actives'* publication era — the single changed variable vs #108 is which
  actives are held out.
  **sha256 `c23ed6bc90fd5c543e07403bd08a89a0ba4cf06f7eeef705e22611c2f0486d95`**

312 molecules total. Any active with no rankable pose is ranked **LAST**, never dropped — the
same anti-inflation rule as every prior gate.

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
  the newest chemistry in the literature, with the whole interval above the bar. Overfitting to
  the development *era* is excluded — the strongest available outcome.
- **SUPPORTED-BUT-WIDE — AUC ≥ 0.70 but CI lower bound < 0.70.** Point estimate passes and the
  permutation test rules out chance, but the interval still admits a true AUC below the bar.
  Temporal generalization is supported, not proven at the bar — reported exactly as such.
- **FAIL — AUC < 0.70.** The ranking result does **not** transfer to newer-literature azoles;
  temporal overfitting cannot be excluded. Reported as a genuine negative. Per pre-registration
  the criterion is **not** re-tuned to recover a pass, and no alternative ranking is substituted
  and reported as a win.

## Honest scope (fixed in advance)

1. **Temporal within one release, not future-prospective.** See "Scope" above — the holdout is a
   later time slice of ChEMBL_37, not a set post-dating it. The fully future-prospective run is
   calendar-gated until ChEMBL_38 and re-runs this protocol unchanged.
2. **Same organism/target world.** Still *C. albicans* azole CYP51 chemistry; this tests
   generalization across **time**, not across a different pathogen or a non-azole class.
3. **Decoys are presumed inactive, not measured** — any real binder among them depresses the
   measured AUC; any decoy-selection bias inflates it. Same limitation as run 2 / #108.
4. **One rigid *C. albicans* receptor**; rigid docking cannot model induced fit; domain ≤ 45
   heavy atoms. Same blind spots as run 2.
5. **This validates a ranking criterion, not a drug.**

## Reproduce

    python scripts/make_temporal_holdout.py years            # (re)pin the year map (network)
    python scripts/make_temporal_holdout.py actives          # offline, deterministic
    python scripts/make_temporal_holdout.py decoys           # fetches ChEMBL
    RECEPTOR=work/receptor.pdbqt scripts/screen.sh work/temporal_ligands work/screen_temporal
    python scripts/validate_gate_temporal.py work/screen_temporal work/receptor_A.pdb
