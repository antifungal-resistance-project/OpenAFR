# Pre-registration — fresh external-holdout validation (issue #82)

**Written 2026-08-30, BEFORE docking a single molecule of this set.** Everything below is
fixed in advance. Read [PREREGISTRATION_run2.md](PREREGISTRATION_run2.md) and
[RESULTS_run2_holdout.md](RESULTS_run2_holdout.md) first — this reuses that protocol,
receptor and grading math unchanged and only swaps in a never-scored active/decoy set.

## The specific objection this closes

Every AUC/EF the project reports was measured on sets the geometric criterion — the
iron-approach filter (N–Fe < 3.0 Å), its 3.0 Å cutoff, and the frozen run-2 docking
protocol — was **developed and repeatedly re-examined against** (run 1 → run 2 holdout →
the verified-inactive and ensemble looks). A skeptic's rebuttal:

> Your criterion and its cutoff were chosen with these molecules in view. AUC 0.794 on the
> run-2 holdout could be the protocol fitting that particular set. Does it hold on data the
> criterion has genuinely never touched?

This pre-registration answers exactly that, and only that. It does **not** re-open run 2 or
change the criterion; it grades the frozen protocol, once, on a clean holdout.

## The holdout (fixed; frozen by sha256 before docking)

Built by `scripts/make_external_holdout.py` under the same provenance discipline as every
other set (`data/ligands/PROVENANCE.md`) — only the provenance of "held out" changes.

- **Actives — `data/ligands/external_actives.smi`, n = 50.** A deterministic seed-**82**
  random sample of the measured-active *C. albicans* azoles in `verified_actives.smi` that
  are **not** in `verified_actives_sample.smi`. The ensemble-confirm look (#9,
  `RESULTS_ensemble_confirm.md`) docked only the seed-42 n=200 sample; the remaining ~3,400
  actives have never been docked or scored by any gate. Every molecule here is genuinely
  never-scored. Applicability domain re-asserted at build time: ≤ 45 heavy atoms, ≥ 1
  aromatic nitrogen (min/mean/max heavy = 18 / 31.0 / 45).
  **sha256 `fe5a696c64914db1b1070dd0a823a9181c7d591368c753b67ee0c33c6e786cd9`**
- **Decoys — `data/ligands/external_decoys.smi`, n = 335.** Property-matched presumed-
  inactives, built by the frozen rule (`openafr.inactives.select_inactives`: MW ± 25,
  cLogP ± 1.5, rotB ± 2, HBD ± 1, HBA ± 2; ≥ 1 aromatic N; Morgan r=2 Tanimoto < 0.35 vs
  every external active; ≤ 45 heavy), fetched fresh from ChEMBL under seed 82 and disjoint
  by ChEMBL id from **every** set the project has used (all actives, both decoy sets, the
  verified inactives, both potency sets — 4,903 ids excluded), so no decoy can secretly be a
  known active or a reused decoy. Achieved match vs the actives: MW 429.8 / 443.1, cLogP
  4.0 / 3.8, heavy 30.4 / 31.0, rotB 5.9 / 6.7.
  **sha256 `32dc7d7a8a8fa411cdec36c880f914108a9878f3f630f0472291a6e4451d64f2`**

385 molecules total. Any active with no rankable pose is ranked **LAST**, never dropped —
the same anti-inflation rule as every prior gate.

## Protocol (fixed; identical to run 2)

- Receptor: rigid 5TZ1 chain A + heme, `work/receptor_A.pdb`
  (sha256 `a1a8410868d1b9ee2835117cf6327ba4721ca3830a4d5689b7e9ea0ee1eee22f`); docking into
  `work/receptor.pdbqt` prepared from it by `scripts/prep_receptor.py`.
- Box + search: from the hash-frozen `protocol.yaml`
  (sha256 `6eca16b1be2e5e33ae42d4b1b7a63a0e9b891a5dd6cbfb19642dd3e9a149aba4`) — box 26 Å @
  70.61/66.28/4.18, exhaustiveness 32, num_modes 20, seed 42. No per-molecule tuning.
- Ligand prep: `obabel … -p 7.4`, exactly as `scripts/screen.sh`.

## Primary criterion, endpoints and pass bar (fixed)

**Primary (gates): mode C** — rank by the minimum azole-nitrogen-to-heme-iron distance over
a molecule's poses (lower = better), the run-2 PRIMARY criterion, reproduced unchanged.
**Secondary (report only): mode A** (best Vina affinity) — the same worse-than-random
baseline contrast run 2 drew.

**Endpoints: AUC and EF@5%** (issue #82). EF@1% is deliberately **excluded** as an endpoint —
[RESULTS_reliability.md](RESULTS_reliability.md) established the top-1% figure is not robust
under seed sampling — though it is printed for continuity. Reported alongside for every mode:
BEDROC(α=20), a one-sided permutation test (20,000 shuffles, seed 42), and a percentile
bootstrap 95% CI over the actives (10,000 resamples, seed 42), from the same
`openafr.scoring` functions every prior gate uses.

**Pass bar: mode-C AUC ≥ 0.70** — the frozen protocol's `gate.min_auc`. (With n = 50 the
EF@1% ≥ 5.0× half of the run-2 bar is not an endpoint here, per the reliability finding.)

## Pre-committed interpretation (three outcomes, fixed before seeing numbers)

Read on mode-C **AUC** at the ≥ 0.70 bar, with the bootstrap CI qualifying it:

- **PASS — AUC ≥ 0.70 and bootstrap 95% CI lower bound ≥ 0.70.** The criterion generalizes
  to data it has never touched, at n = 50 (≈ 7× the run-2 holdout), with the whole interval
  above the bar. Overfitting to the run-2 set is excluded; the headline AUC ~0.79 result
  stands on a genuinely independent, well-powered holdout. This is the strongest available
  outcome.

- **SUPPORTED-BUT-WIDE — AUC ≥ 0.70 but CI lower bound < 0.70.** The point estimate passes
  and the permutation test rules out chance, but the interval still admits a true AUC below
  the bar. Generalization is supported, not proven at the bar — reported exactly as such, the
  same honesty run 2 applied to its n=7 CI.

- **FAIL — AUC < 0.70.** The run-2 ranking result does **not** transfer to a never-scored
  set; overfitting to the development sets cannot be excluded. Reported as a genuine negative.
  Per pre-registration the criterion is **not** re-tuned to recover a pass, and no alternative
  ranking is substituted and reported as a win.

## Honest scope (fixed in advance)

1. **Same curated population.** The actives are a never-*scored* subset of the same
   ChEMBL *C. albicans* azole population as `verified_actives`; this is a clean
   never-touched-by-any-gate holdout, **not** a different chemical world or a later time
   slice. A true temporal/prospective holdout would need a newer ChEMBL release (recorded as
   the stronger follow-up).
2. **Decoys are presumed inactive, not measured** — any real binder among the 335 depresses
   the measured AUC; any decoy-selection bias inflates it. Same limitation as run 2.
3. **One rigid *C. albicans* receptor**; rigid docking cannot model induced fit; domain
   ≤ 45 heavy atoms. Same blind spots as run 2.
4. **This validates a ranking criterion, not a drug.**

## Reproduce

    python scripts/make_external_holdout.py actives          # offline
    python scripts/make_external_holdout.py decoys           # fetches ChEMBL
    RECEPTOR=work/receptor.pdbqt scripts/screen.sh work/external_ligands work/screen_external
    python scripts/validate_gate_external.py work/screen_external work/receptor_A.pdb
