# Temporal holdout — pre-registered publication-date test (issue #82)

Date: 2026-08-31
Pre-registration: `work/PREREGISTRATION_temporal.md`
(sha256 `33ff9558...`, verified unmodified by the grading script at run time)
Protocol: `protocol.yaml` (sha256 `6eca16b1...`) — box 26 Å @ 70.61/66.28/4.18,
exhaustiveness 32, num_modes 20, seed 42, rigid 5TZ1 chain A + heme. No per-molecule tuning.
Test set: **50 newest-literature *C. albicans* azole actives (earliest antifungal document
year ≥ 2024) + 262 property-matched decoys** (312 total). The criterion was defined and
validated on the classic approved azoles (fluconazole … efinaconazole/luliconazole, approved
into the 2010s); this holds out azoles that entered the antifungal literature in 2024–2025 and
that no gate has ever docked. 13 molecules produced no dockable pose (**0 actives**, 13 decoys)
and are ranked **LAST**, never dropped.

## RESULT: SUPPORTED-BUT-WIDE — the criterion generalizes across time; point estimate clears the bar

Pre-registered bar: **mode-C AUC ≥ 0.70**, qualified by the bootstrap 95% CI per the
three-outcome interpretation fixed in advance (reused verbatim from #108).

| mode | AUC | EF@1% | EF@5% | BEDROC | role |
|---|---|---|---|---|---|
| **C — iron-approach distance** | **0.754** | 4.68x | **3.90x** | 0.553 | PRIMARY — point estimate passes |
| A — Vina score alone | 0.715 | 0.00x | 1.17x | 0.217 | secondary — worse on both endpoints |

The geometric criterion, **frozen** and applied with no re-tuning, ranks 50 azoles from the
2024–2025 literature — chemistry that did not exist when the criterion was built on old approved
drugs — above 262 property-matched decoys at **AUC 0.754**, above the 0.70 bar, with **10 of the
top 15** molecules being actives. The criterion generalizes forward in time, not only to a
random never-scored subset (#108, AUC 0.715): a criterion defined on decades-old approved azoles
still separates the newest research azoles from decoys. The point estimate is in fact slightly
higher than #108's and the run-2 headline is intact.

## Statistical support (mode C)

    permutation test (20,000 label shuffles): p < 0.0001   (0 of 20,000 shuffles matched)
    bootstrap 95% CI over actives:            AUC 0.682 - 0.823

Chance is decisively excluded. Per the pre-registered three-outcome rule this is
**SUPPORTED-BUT-WIDE**, not PASS: the point estimate (0.754) clears the bar and the permutation
test rules out chance, but the bootstrap CI lower bound (0.682) dips just below 0.70 at n=50, so
the interval still admits a true AUC marginally under the bar. This is the same honest read #108
applied to its n=50 CI (0.630–0.795), reported as such rather than rounded up to a pass. The CI
is narrower here (0.682–0.823 vs 0.630–0.795) because the actives separate slightly more cleanly.

## Geometry beats the scoring function — replicated on later chemistry

On the **identical poses**, mode C (0.754) beats raw Vina affinity (0.715) on AUC and beats it
decisively on the enrichment endpoints that matter for triage (EF@5% 3.90× vs 1.17×; 10 vs 3
actives in the top 15; BEDROC 0.553 vs 0.217). The "geometry beats the scoring function"
contrast run 2 drew (0.794 vs 0.471) survives on a temporally disjoint set — the margin on AUC
alone is narrower here, but the early-enrichment gap the tool is actually used for is large.

## Honest scope (from the pre-registration)

1. **Temporal within one release, not future-prospective.** This is a publication-year split of
   ChEMBL_37 (2026-05-01), the release the corpus was pulled from and — confirmed via the ChEMBL
   `/status` API on 2026-08-31 — the latest release. A set drawn from depositions that *post-date*
   the release is calendar-gated for everyone until ChEMBL_38 and cannot be built by any means
   today (the data does not exist yet; a compute host does not change that). Re-running
   `scripts/make_temporal_holdout.py years` and re-freezing on ChEMBL_38 lifts this to the fully
   future-prospective variant under the identical protocol.
2. **Decoys are presumed inactive, not measured** — any real binder among the 262 depresses the
   measured AUC; any decoy-selection bias inflates it. Same limitation as run 2 / #108.
3. **One rigid *C. albicans* receptor; domain ≤ 45 heavy atoms.** Same blind spots as run 2.
4. **This validates a ranking criterion, not a drug.**

## Reproduce

    python scripts/make_temporal_holdout.py years            # (re)pin the year map (network)
    python scripts/make_temporal_holdout.py actives          # offline, deterministic
    python scripts/make_temporal_holdout.py decoys           # fetches ChEMBL
    RECEPTOR=work/receptor.pdbqt scripts/screen.sh work/temporal_ligands work/screen_temporal
    python scripts/validate_gate_temporal.py work/screen_temporal work/receptor_A.pdb
