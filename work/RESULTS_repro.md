# Reproducibility / drift harness (#89)

Date: 2026-08-28.

The project's whole claim to credibility is *provable discipline*: a hash-frozen
protocol, hash-frozen pre-registrations, and a confidence stack whose every published
number derives from committed inputs. The candidate-confidence batch (#65–#88) stacked
eleven per-candidate axes into one scorecard — and with every axis added, the surface
that could silently drift grew: a pass bar quietly tuned in `protocol.yaml`, a
pre-registration reworded after it was frozen, a per-axis JSON regenerated under changed
logic without rebuilding the card. Any one of those would break a claim with no visible
signal.

[`scripts/repro.py`](../scripts/repro.py) re-checks that entire frozen surface in **one
command** and **fails loudly (exit 1) on any drift**, so CI can gate every PR on it
([`.github/workflows/repro.yml`](../.github/workflows/repro.yml)).

## What it checks

| # | check | what it catches |
|---|---|---|
| 1 | **Frozen protocol integrity** | `protocol.yaml` (docking box + search effort + the gate's pass bar — the numbers easiest to tune into a pass) matches its frozen SHA-256. |
| 2 | **Frozen pre-registration integrity** | Every `work/PREREG_*.sha256` verifies: the pre-registration markdown **and** the frozen input datasets each one pins (held-out SMILES, decoy sets). A drift here means the rules or the data changed after freezing. |
| 3 | **Docking-free derivation reproduces** | The consolidated scorecard (`work/candidates/shortlist_scorecard.json`, the capstone the whole stack feeds) re-derives **byte-for-byte** from its committed per-axis inputs. Catches a per-axis JSON that drifted, or merge logic that changed without rebuilding the card. |

Current state (2026-08-28): **all three green — 18 pre-registrations / 39 frozen files
unmodified, protocol frozen, scorecard reproduces byte-for-byte.**

```
==============================================================================
OPENAFR REPRODUCIBILITY / DRIFT HARNESS (#89)
==============================================================================
[PASS] 1. Frozen protocol integrity
  protocol.yaml unmodified (sha256 6eca16b1be2e5e33...)
[PASS] 2. Frozen pre-registration integrity
  18 pre-registrations / 39 frozen files all unmodified
[PASS] 3. Docking-free derivation reproduces
  scorecard re-derives byte-for-byte from committed per-axis inputs
------------------------------------------------------------------------------
PASS: all 3 checks green — no silent drift in the frozen surface.
```

## What it deliberately does NOT do

**Re-docking is out of scope, by design.** The ~43 MB of Vina pose trees are gitignored
(`work/screen*/`, `work/candidates/dock_*/`, …) and regenerate from `scripts/` — each
`RESULTS_*.md` carries its own "Reproduce" block. Re-running them needs the docking
toolkit and hours of compute, which is not a per-PR CI job. So the **pose-dependent**
per-axis JSONs (`shortlist_confidence.json`, `shortlist_resistance.json`, …) are treated
here as frozen committed inputs:

- check 3 proves the scorecard still derives from them **unchanged**, and
- check 2 proves the **rules that graded them** were not moved after the fact.

Together that guarantees no *silent* drift in anything that does not require docking. A
genuine re-dock — regenerating a pose tree and re-grading it through a
`scripts/validate_gate*.py` — remains a separate, documented, manual step. This harness
is honest about that boundary rather than implying it re-runs Vina.

The per-gate scripts already re-check their own single pre-registration + `protocol.yaml`
hash at run time (soft warnings). This harness is the **systematic superset**: every
frozen file at once, as a hard exit-1 gate, so nothing frozen can drift unnoticed between
runs of the individual gates.

## Reproduce

    conda run -n openafr python scripts/repro.py          # exits 1 on any drift
    conda run -n openafr python scripts/repro.py --quiet   # only the verdict + failures

Logic is unit-tested in [`tests/test_repro.py`](../tests/test_repro.py): each check passes
on the real tree and fails when its target is tampered (a drift-detector that cannot detect
drift is worse than none).
