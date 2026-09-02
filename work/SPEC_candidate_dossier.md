# SPEC — Candidate Dossier generator

Status: draft spec (no code yet)
Pipeline stage: `report/` (runs after rank, alongside `report_candidates.py`)
Origin: the one item worth pulling forward from the 2026-09-02 "OpenAFR → drug program"
essay — turn a ranked library into a small set of machine-readable, per-molecule
**candidate dossiers**, each carrying an experimental-priority verdict and an explicit,
sourced reason for selection.

---

## 1. Why this exists (and what it is NOT)

Today a candidate's evidence is scattered across five stages that each emit their own
file, keyed by molecule name, that a human joins by hand:

| Evidence | Produced by | Output today |
|---|---|---|
| Geometry rank (min N–Fe, pose convergence, docking score) | `rank_candidates.py` / `analyze_poses.analyze` | TSV |
| Prior measured activity (ChEMBL / PubChem / BindingDB) | `check_precedent.py` → `precedent.render` | TSV |
| Physchem / drug-likeness / structural alerts | `admet_filter.py` → `admet.profile` | TSV/JSON |
| Human-CYP off-target liability | `cyp_offtarget_screen.py` → `cyp_offtarget.profile` | TSV/JSON |
| Synthesizability / availability | `synth.profile` | (in shortlist tools) |

`report_candidates.py` already fuses **geometry + precedent** into human-readable
Markdown. This spec adds the missing piece: **one machine-readable JSON dossier per
molecule** that fuses *all five*, assigns an A/B/C experimental priority from documented
rules, and states the selection hypothesis in one sentence.

**Non-goals (hard boundaries):**
- Computes **nothing new** about geometry. Reuses `analyze_poses.analyze` output; never
  re-docks, never re-runs the gate. Same honesty banner as `report_candidates.py`: a
  dossier is a *hypothesis for a wet lab, never a hit*, and means nothing unless
  `validate_gate2.py` has PASSED under the frozen protocol.
- Not a ranker. Ordering stays the validated criterion (min N–Fe, unrankable last).
  Priority A/B/C is an **orthogonal annotation**, not a re-rank.
- Not a filter/gate. It never drops a molecule; a bad profile lowers priority, it does
  not delete the row (same rule as the rest of the pipeline).
- No new third-party dependency. RDKit + stdlib only, matching the repo.

---

## 2. Module boundary

Two new files, mirroring the repo's split (pure logic in `openafr/`, CLI wiring in
`scripts/`), plus tests.

- `openafr/dossier.py` — **pure, network-free, deterministic.** Fuses already-computed
  per-molecule dicts into one dossier dict; assigns priority; writes the rationale. No
  I/O, no RDKit calls of its own beyond what the profilers already return. This is the
  unit-tested core.
- `scripts/build_dossier.py` — CLI. Reads the ranked screen output + the `.smi` library,
  calls `admet.profile` / `cyp_offtarget.profile` / `synth.profile` locally (all local,
  deterministic), optionally calls the network precedent path, and hands the assembled
  per-molecule dicts to `dossier.build`. Emits JSON.
- `tests/test_dossier.py` — priority rules + schema + provenance stamping.

`dossier.build` takes **already-computed** inputs so it stays pure and testable; the CLI
owns all the calling of profilers and the network.

---

## 3. The dossier schema (one JSON object per molecule)

```jsonc
{
  "name": "OAFR-000123",              // stable candidate id (from .smi name column)
  "smiles": "…",                      // canonical (RDKit) SMILES
  "inchikey": "…",                    // precedent.inchikey(smiles) — the join key

  "geometry": {                       // from analyze_poses.analyze (reused, not recomputed)
    "rank": 3,                         // position in the validated order (1 = best)
    "best_fe": 2.41,                   // closest N–Fe approach (Å); null = unrankable
    "n_passing": 7, "n_poses": 9,      // pose convergence
    "best_score": -8.1                 // docking score, carried DEMOTED (ranked worse than random here)
  },

  "precedent": {                      // from precedent.render (null if --no-precedent)
    "verdict": "no-precedent",         // tested-active | tested-inactive | tested-ambiguous | no-precedent | …
    "phenotypic": null,                // secondary whole-cell fungal verdict or null
    "n_off_target": 4,                 // novelty signal
    "nearest_tested": {"tanimoto": 0.72, "name": "…", "verdict": "tested-inactive"}  // or null
  },

  "developability": {                 // from admet.profile
    "mw": 341.4, "clogp": 3.1, "tpsa": 62.0,
    "ro5_violations": 0, "veber_ok": true,
    "esol_logs": -4.2,
    "pains": [], "structural_alerts": [["BRENK", "…"]], "n_alerts": 1,
    "herg_risk_factors": false
  },

  "selectivity": {                    // from cyp_offtarget.profile
    "azole_warhead": "triazole",       // the antifungal pharmacophore (or null)
    "flags": "type-II(triazole):HIGHx5",  // cyp_offtarget.flag_summary
    "n_high": 5, "n_moderate": 0,
    "flagged_isoforms": ["CYP3A4","CYP2C9","CYP2D6","CYP2C19","CYP1A2"]
  },

  "synthesis": {                      // from synth.profile
    "sa_score": 3.1, "sa_band": "easy",
    "availability": "catalogued-drug", // catalogued-drug | pubchem-record | no-pubchem-record | not-checked
    "pubchem_cid": 60795
  },

  "priority": "B",                    // A | B | C  — see §4
  "priority_reasons": [               // the individual rule hits, worst-first, human-readable
    "known azole warhead → pan-CYP HIGH liability (selectivity risk)",
    "no measured antifungal precedent (novel, but absence ≠ untested)"
  ],
  "selection_hypothesis": "Ranks #3 by N–Fe geometry with convergent poses, no prior fungal-CYP51 measurement, and clean drug-likeness; the triazole warhead is the selectivity liability to watch.",

  "provenance": {                     // the skeptical-reader block, stamped by the CLI
    "protocol_sha256": "…",            // protocol.yaml hash (protocol module, as report_candidates does)
    "protocol_unmodified": true,
    "receptor": "work/receptor_A.pdb",
    "openafr_version": "…",            // git describe / package version
    "generated_utc": "2026-09-02T…Z",
    "precedent_source": "both"         // or "none"
  }
}
```

Output is a **JSON array** written to `-o` (default stdout), ordered by geometry rank.
Optionally `--split DIR` writes one `<name>.json` per molecule for downstream provenance
(the `experiments/candidate_registry/` idea, without committing to a LIMS yet).

---

## 4. Priority rules (A / B / C)

Priority answers one question: *is this worth spending real bench money on, and why?*
It is a small, **documented, deterministic** decision table over fields already in the
dossier — never a black box, never a re-rank. Start conservative; every rule is a one-line
`priority_reasons` entry so the label is auditable.

**Automatic C (deprioritize — do not spend money re-testing / obvious liability):**
- `precedent.verdict == "tested-inactive"` — the file-drawer catch; re-testing a known
  negative wastes the experiment.
- `precedent.nearest_tested.verdict == "tested-inactive"` at Tanimoto ≥ 0.85 — a soft
  negative (near-identical compound already failed).
- `geometry.best_fe is null` — unrankable (no usable pose); no geometric hypothesis to test.
- `developability.pains` non-empty — assay-interference risk makes the readout unreliable.

**Automatic A (test first) requires ALL of:**
- `geometry.rank` in the top decile of the ranked set **and** `n_passing ≥ 2` (convergent).
- `precedent.verdict` in {`no-precedent`, `tested-ambiguous`} (not already a known active —
  a known active is real but not a discovery; route those to a separate "reference" bucket).
- `developability.ro5_violations ≤ 1` and `n_alerts == 0` (clean, buyable-drug-like).
- `synthesis.availability` in {`catalogued-drug`, `pubchem-record`} (physically obtainable).

**Everything else → B** (worth testing, with a caveat named in `priority_reasons`).

Two flags that annotate but do **not** by themselves change tier (they inform, they
don't gate — consistent with `cyp_offtarget_screen.py`'s stance):
- azole warhead / high pan-CYP liability → always logged as a `priority_reasons` line so
  the selectivity risk travels with the candidate.
- `precedent.verdict == "tested-active"` → tagged `reference` in reasons (a positive
  control for the panel, not a novel hit).

All thresholds live as named constants at the top of `dossier.py` (like `admet.RO5_*`,
`synth.SA_EASY`), documented inline, so the table is one place and easy to tune once real
experimental results come back and tell us which rules actually predicted hits.

---

## 5. CLI

```
python scripts/build_dossier.py SCREEN_DIR LIBRARY.smi [RECEPTOR_PDB] \
    [-o dossiers.json] [--split DIR] [-n TOP] \
    [--precedent both|chembl|pubchem|all|none] [--near TESTED.smi]

  SCREEN_DIR      directory of docked *.pdbqt poses (e.g. work/screen)  — geometry, reused
  LIBRARY.smi     the SAME library, `<SMILES>\t<name>` per line (admet_filter.read_smi format)
  RECEPTOR_PDB    receptor with the heme iron (default work/receptor_A.pdb)
  -o / --out      JSON array out (default stdout)
  --split DIR     also write one <name>.json per molecule
  -n / --top      cap dossiers to the top-N ranked (default: all)
  --precedent     network precedent source; 'none' skips the network entirely (default: both)
  --near          nearest-tested-neighbour annotation (soft negatives), passthrough to precedent
```

- Geometry comes from `collect()` (lift/reuse the one in `report_candidates.py` rather than
  re-implement — factor it out if cleaner).
- `admet` / `cyp_offtarget` / `synth` profiles are all **local + deterministic**; always run.
- Precedent is the only network stage; `--precedent none` produces a fully offline,
  reproducible dossier (precedent block = null, priority computed without it).
- Fails per-molecule, never per-run: an unparseable SMILES or a precedent lookup error is
  recorded on that molecule (`"error": "…"`, priority `C`) and the run continues — same
  contract as `check_precedent.py`.

---

## 6. Tests (`tests/test_dossier.py`)

Pure-function core makes this cheap. Cover:
- **Schema**: `build` on a hand-made set of profile dicts yields every documented key;
  JSON round-trips.
- **Priority table**: one case per rule — tested-inactive → C; PAINS → C; null best_fe →
  C; the full all-of-A conjunction → A; a single failed A condition → B; near-neighbour
  soft-negative → C.
- **Reasons travel**: azole warhead always appears in `priority_reasons` regardless of
  tier; tested-active is tagged `reference`, not A.
- **Ordering**: output array stays in geometry-rank order independent of priority.
- **Offline path**: `--precedent none` (precedent=null) still assigns a valid priority.
- **Determinism**: same inputs → byte-identical JSON (sorted keys, fixed float rounding).
- **Per-molecule failure isolation**: one bad SMILES doesn't sink the batch.

---

## 7. What this deliberately leaves for later

Named so scope stays honest (these are the essay's later phases, correctly sequenced
*after* a first experimental result exists — see memory `outreach-hold-build-first`):
- No `experiments/` provenance graph, no openBIS/LIMS, no assay definitions.
- No IP/patent screening (needs a data source we don't have wired).
- No learning loop / model retraining — there is nothing to learn from until a partner
  returns measurements.
- No commercial-vendor availability beyond the existing PubChem-record proxy in `synth`.

The dossier is the artifact those later stages would attach to; building it now costs
little because it is assembly over modules that already exist, and it makes the
"10–20 defensible candidates, each with a reason" shortlist a generated output instead of
a hand-built table.
```
