# Stage 0 result — the resistance-active pool with the domain ceiling lifted

**Run 2026-08-23, local, rdkit-only (no docking).** Reproduce with:

    conda run -n openafr python scripts/domain_lifted_pool_audit.py

This is the cheap, potentially-fatal gate on the flexible/ensemble-receptor work
(`docs/designs/flexible-receptor-domain-ceiling.md`, Stage 0). It asks one question, before
any receptor is built: **if a flexible receptor lifted the ≤45-heavy-atom domain cap, would
the resistance-active validation pool grow large enough to grade a resistance-breaker gate?**
The design's working bar is n≥4 (below that, a differential-activity permutation test has no
resolution — same reasoning that made n=7 "too thin" in T1).

## Method

For each candidate resistance-retaining CYP51 inhibitor
(`data/ligands/resistance_active_candidates.tsv`, SMILES copied verbatim from PubChem), the
audit computes — reusing the project's existing definitions, not new ones:

- **heavy-atom count** (the applicability-domain axis, `scripts/filter_library.py`);
- **max Morgan r=2 / 2048-bit Tanimoto** vs the 15-compound active panel (8 training + 7
  holdout) — the novelty rule from `PROVENANCE.md` constraint #2 (must be < 0.70);
- **C. auris evidence** — a primary-literature PMID for measured activity against
  azole-resistant *C. auris*, or `none`.

and classifies each into the pool under two regimes:

- **RIGID** = evidence AND novel AND ≤45 heavy atoms (today's tool);
- **DOMAIN-LIFTED** = evidence AND novel, no size cap (what a flexible receptor would buy).

**Built-in cross-check:** the audit re-derives the exact novelty values previously recorded by
hand in `PROVENANCE.md` — VT-1598 0.312, oteseconazole 0.333, quilseconazole 0.358,
ravuconazole 0.879. Matching those independently confirms the fingerprint method here is the
same one the prior curation used.

## Result

| compound | heavy | max Tc | nearest | C. auris evid. | novel <0.70 | rigid pool | lifted pool |
|---|---|---|---|---|---|---|---|
| VT-1598 | 43 | 0.312 | fluconazole | PMID:30530603 | yes | **IN** | **IN** |
| opelconazole (PC945) | 50 | 0.553 | posaconazole | PMID:31325309 | yes | — (over-domain) | **IN** |
| oteseconazole (VT-1161) | 37 | 0.333 | fluconazole | none | yes | — (no evidence) | — (no evidence) |
| quilseconazole (VT-1129) | 36 | 0.358 | fluconazole | none | yes | — (no evidence) | — (no evidence) |
| ravuconazole | 31 | 0.879 | isavuconazole | PMID:33609301 | **no** | — (not novel) | — (not novel) |

- **RIGID pool: n=1** — {VT-1598} (confirms the T1 finding).
- **DOMAIN-LIFTED pool: n=2** — {VT-1598, opelconazole}.
- **Unlocked by lifting the domain cap: n=1** — {opelconazole}.

## Verdict

**Lifting the domain ceiling does NOT reach a testable pool (n=2 < 4).** It buys exactly one
compound — opelconazole, which is excluded today *only* for being 50 heavy atoms. Everything
else stays excluded for a reason a bigger pocket cannot fix:

- oteseconazole, quilseconazole — no *C. auris*-specific MIC in the literature (their data is
  *C. albicans* / *C. glabrata* / *Cryptococcus*). Confirmed again 2026-08-23.
- ravuconazole — not novel (Tanimoto 0.879; it is essentially isavuconazole, differing only in
  the difluorophenyl substitution pattern — the PubChem SMILES makes this explicit).

**So the binding constraint is the compound universe, not (only) the rigid receptor.** The set
of CYP51 inhibitors that are simultaneously (a) verified structure + PubChem CID, (b) measured
active against azole-resistant *C. auris*, and (c) novel vs the azole panel is n=2 whether or
not we can dock big ligands. A flexible receptor is therefore **not justified on the
resistance-breaker-validation rationale alone** — it would double the pool from 1 to 2 and stop
there, still short of a gradeable differential-activity gate.

This is exactly the outcome Stage 0 was placed first to catch: it closes the flexible-receptor
line of work *before* any receptor code or cloud run was spent on it.

## Scope / honesty caveats

- **Not a proof of absence.** This is a focused literature pass (2026-08-23), not an exhaustive
  one. Additional verified *C. auris*-active CYP51 inhibitors may exist; if any are found, add a
  row to `resistance_active_candidates.tsv` and re-run — the verdict flips automatically at n≥4.
  But the *registered + C.auris-verified + novel* intersection is demonstrably thin: the
  standard reviews of anti-*C. auris* CYP51 agents surface the same short list, and most recent
  medicinal-chemistry hits are unregistered research compounds ("compound 7", "V23") that fail
  the verified-SMILES discipline (constraint #3).
- **Novelty and domain thresholds are the pre-registered ones** (0.70, 45 heavy), not tuned
  here. Changing either is a conscious, recorded protocol change.
- **This says nothing about the decoy ceiling** (the second, independent ceiling in the design
  §"two ceilings"). Even if the pool were large enough, a flexible protocol would still have to
  prove it separates actives from mechanism-matched decoys. Stage 0 fails before that question
  is reached.

## What this means for the next move

The design's Stage 0 was written with an explicit stop condition, and it fired. The honest
options now (for the user to choose — no receptor build is warranted by this result):

1. **Shelve the flexible-receptor build.** It cannot deliver a resistance-breaker gate while the
   verified compound pool is n=2. Record the negative; don't spend the build.
2. **Change what the pool is used for.** n=2 cannot support a permutation test, but VT-1598 +
   opelconazole *could* serve as a tiny confirmatory probe ("does the criterion rank the two
   known resistance-actives above the defeated panel in the Y132F pocket?") — a weak,
   anecdotal check, explicitly not a validation. Low yield; the T1 decision already
   de-prioritized the n=1 version of this.
3. **Reconsider the target/strategy** given that the resistance-breaker framing is
   compound-starved on CYP51 — e.g. the parked "method-as-a-public-benchmark" route
   (`refocus-resistant-cyp51-triage.md`, lateral option), which does not need a large
   in-house resistance-active set.

## Action taken (2026-08-23)

Following option 3's direction-agnostic first step, added **`scripts/applicability_triage.py`**
(+ `tests/test_applicability_triage.py`, 7 known-answer tests): the outward-facing front door
that answers, for any candidate SMILES and without a docking run, *does the validated method
even apply to this molecule, and is it novel?* It reuses the pre-registered envelope rules
(≤45 heavy atoms + ≥1 aromatic N, imported from `filter_library.py`) and the PROVENANCE novelty
rule (Morgan r=2, Tc < 0.70 vs the 15 validated actives), returning one honest per-molecule
verdict: `in-envelope-novel` / `near-known` / `out-of-domain` / `out-of-scope` / `unparseable`.
This is the atom of the method-as-benchmark route and is useful under every future (it is what
lets a collaborator try their own compound), while committing to no receptor build. It correctly
reports opelconazole (50 heavy) out-of-domain, consistent with this audit. Full suite green
(278 tests, coverage 97.23%).

**Sources (C. auris activity evidence):** VT-1598 — Wiederhold et al. 2018,
[PMID 30530603](https://pubmed.ncbi.nlm.nih.gov/30530603/). Opelconazole/PC945 — Rivero-Menendez
et al. 2019, [PMID 31325309](https://pubmed.ncbi.nlm.nih.gov/31325309/). Ravuconazole (C. auris
activity; excluded on novelty) — [PMID 33609301](https://pubmed.ncbi.nlm.nih.gov/33609301/).
