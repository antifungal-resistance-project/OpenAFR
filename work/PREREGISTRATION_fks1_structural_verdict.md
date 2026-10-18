# Pre-registration — FKS1 structural contact label (moonshot #4, Rung 3)

> **STATUS: FROZEN 2026-10-18, before the cert grader was run and before any verdict
> code was written.** This document fixes the claim, the cert truth set, the pass bar,
> and the non-claims for `openafr.fks1_structural.structural_verdict`. It crosses the
> FKS1 track's declared "detection-only, no structural label" boundary by adding a
> spatial context label for drug-contact-zone VUS tokens, so under this project's
> proof-before-claim discipline it must declare its scope in advance. The repro gate
> (`scripts/repro.py`) re-checks this document's SHA-256 via
> `work/PREREG_fks1_structural_verdict.sha256` and fails on any drift.

## Why this exists — Rung 2 gives distance; Rung 3 gives a certifiable label

Rung 2 (`openafr/fks1_structural.py`) added the minimum heavy-atom distance from each
novel FKS1 hot-spot position to enfumafungin in the 9WZU cryo-EM structure. The proximity
tier (`drug_contact` ≤5 Å, `pocket_adjacent` 5–8 Å, etc.) is reported on the weather page
as a bare number — "3.7 Å from drug" — with no structural interpretation attached. The
`structural` field on FKS1 alerts stays null by design.

Rung 3 takes the next step: it attaches a certifiable **structural contact label** to
drug-contact-zone VUS variants — *"at drug-contact residue (9WZU cryo-EM, ±1–2 Å proxy)"*
— backed by a pre-registered cert that validates the contact map's structural relevance
anti-circularly against the known-resistance panel.

This is NOT a docking protocol and NOT a resistance prediction. The label claims spatial
co-location only. ACT-NOW remains unreachable for FKS1.

## The claim

`structural_verdict(token)` returns a **structural contact label dict** for FKS1 VUS
tokens in `drug_contact` or `pocket_adjacent` proximity tier, or `None` for all others.
It is a spatial context label, not a resistance prediction or a binding-affinity estimate.

Verdict shape (mirrors the ERG11 structural verdict keys for pipeline compatibility, but
with a weaker `evidence_class` that prevents it from ever driving ACT-NOW):

```python
{
  "evidence_class": "proximity_label",        # NOT "docking" or "measured"
  "proximity_tier": "drug_contact",           # or "pocket_adjacent"
  "min_dist_to_drug_A": 3.8,                  # float, from contact map
  "label": "at drug-contact residue",         # or "adjacent to drug-contact residues"
  "note": "<one non-claim sentence>",
}
```

Returns `None` for: core panel markers (already carry `RESISTANCE_MARKER_DETECTED`),
`pocket_proximate`, `distal`, `no_structure`, out-of-window positions, malformed tokens.

## The cert — anti-circular validation against the known-resistance panel

### Truth set

All **core-tier** panel POSITIONS (residue numbers, not specific mutations) from
`openafr.fks1_caller.FKS1_WINDOWS`, excluding the `LOW_PPV_PANEL` (D642, M690):

| Position | Window | Role | Source |
|---|---|---|---|
| 635 | HS1 | F635{C,Y} resistance marker | Perlin FKS1 HS1 literature |
| 639 | HS1 | S639{F,P,Y} resistance marker | Perlin FKS1 HS1 literature |
| 691 | HS3 | W691L resistance marker | AAC 2023 aac.00423-23 (CRISPR-confirmed) |
| 1354 | HS2 | R1354S resistance marker | Perlin FKS1 HS2 literature |

**Anti-circularity:** These four positions were identified as echinocandin-resistance
hotspots from clinical MIC studies on *C. auris* and other *Candida* spp. The 9WZU
contact map was derived entirely from cryo-EM drug coordinates in *S. cerevisiae* —
independent source, no MIC data used.

### Pre-committed predictions (before running the cert grader)

All four pre-committed values are derivable from the Rung 2 preregistration
(`work/PREREGISTRATION_fks1_structural_proximity.md`, frozen 2026-10-17), which pinned
the same contact map. Re-stating them here makes the Rung 3 pass bar concrete and
verifiable before the cert script runs.

| Position | Expected proximity tier | Expected dist (Å) |
|---|---|---|
| 635 | `drug_contact` | ~2.94 Å |
| 639 | `pocket_adjacent` | ~7.16 Å |
| 691 | `drug_contact` | ~3.66 Å |
| 1354 | `pocket_adjacent` | ~7.53 Å |

All four are in `drug_contact` or `pocket_adjacent` → predicted cert outcome: **100% (4/4)**.

### Pre-committed pass bar

| Outcome | Condition |
|---|---|
| **PASS** | ≥ 80% of core panel positions in `drug_contact` OR `pocket_adjacent` |
| **FAIL** | < 80% of core panel positions in `drug_contact` OR `pocket_adjacent` |
| **UNDERPOWERED** | Fewer than 3 core panel positions have a contact-map entry (absent/`no_structure`) |

The 80% bar is pre-committed above because: with 4 positions and a pre-registered
prediction of 100%, a single structural surprise (one core position turning out `distal`)
would still PASS; two would FAIL. This is the right gate — one anomaly is tolerable,
two indicate the alignment is wrong.

### Secondary observation (exploratory, not part of pass bar)

Low-PPV panel positions (D642, M690 — positions with uncertain clinical PPV, treated as
`UNCHARACTERIZED_VARIANT` abstentions in the verdict layer) are predicted to fall in
`pocket_proximate` or `distal` based on the same contact map. If confirmed, this would be
a structurally coherent result: the positions with strongest resistance association
(core panel) cluster at the drug-contact face, while positions with uncertain association
(low-PPV) sit in the outer binding cavity. This is descriptive only — the PPV
classification was made before this analysis and will not be changed based on it.

## The six things this module DOES NOT claim

1. **No resistance prediction.** The label states spatial co-location with the drug in the
   cryo-EM structure. It does not predict whether a substitution at this position will
   confer echinocandin resistance. A novel `drug_contact` VUS may or may not be resistant.
2. **No ΔΔG or binding-affinity estimate.** No FKS1 docking protocol exists; the cryo-EM
   proximity distance is not an energy score.
3. **No probability.** `evidence_class: "proximity_label"` is categorical; no confidence
   value is assigned. This is consistent with the blocked calibration track (#136).
4. **ACT-NOW remains unreachable for FKS1.** `_fks1_priority` is NOT extended to
   `act-now` in Rung 3. A drug-contact VUS alert may reach `watch` tier but never
   `act-now` — that tier requires a full fit-consequence verdict (`evidence_class:
   "measured"` or `"docking"`), which this preregistration explicitly does not claim.
5. **S. cerevisiae distances are a proxy.** 9WZU is a Sc structure; Ca distances carry
   ±1–2 Å uncertainty from species divergence (same caveat as Rung 2). Stated on every
   label.
6. **`pocket_adjacent` is a weaker claim than `drug_contact`.** 5–8 Å is adjacent to
   the drug-contact zone, not a direct drug contact. The label wording is calibrated to
   this difference.

## Verdict contract and alert-priority behavior

`structural_verdict` returns a dict with `evidence_class: "proximity_label"`. The ERG11
pipeline uses `evidence_class: "measured"` or `"docking"` to reach ACT-NOW. This new
class must be treated as insufficient to cross the ACT-NOW gate — `_fks1_priority` in
`openafr/alert.py` must guard on `evidence_class` when it eventually uses the structural
field, not just on the field's presence.

For non-backlog FKS1 signals with a drug-contact or pocket-adjacent structural verdict,
the alert tier stays `watch` (the maximum FKS1 tier). For backlog signals, tier stays
`context` regardless of structural verdict.

## Verification

No GCP required. Everything is deterministic:

1. `python scripts/certify_fks1_structural_verdict.py` → must print **VERDICT: PASS**
   and write `work/RESULTS_fks1_structural_verdict.md`.
2. `python -m pytest tests/test_fks1_structural.py` → new tests must pass:
   - `structural_verdict` returns non-None for a VUS at a pre-committed drug-contact
     position (e.g. a hypothetical `F635X` would be skipped as core; a novel VUS like
     `L638X` at a drug-contact position is the real test case).
   - `structural_verdict` returns None for distal/no_structure positions.
   - Core panel tokens return None (they already have a resistance verdict).
   - Returned dict carries required keys (`evidence_class`, `proximity_tier`,
     `min_dist_to_drug_A`, `label`, `note`).
3. `python scripts/render_weather.py --watching-only --out /tmp/rung3_test.html` →
   page renders without error.
4. `python scripts/repro.py` → drift gate passes including this prereg's SHA.

## Implementation (gated on PASS)

- `openafr/fks1_structural.py`: `structural_verdict(token) -> dict | None`
- `openafr/fks1_position.py`: add `structural_verdict` key to `position_context` output
- `openafr/alert.py`: set `base["structural"] = sv` for drug-contact FKS1 alerts
- `openafr/weather.py`: render structural label on VUS cards

## What stays for a future Rung 4

- **A docking protocol for FKS1**: no heme-iron analogue exists, so a new protocol must
  be designed, pre-registered, and its non-claims frozen. Would enable a `"docking"`
  evidence class and open the ACT-NOW tier for FKS1.
- **A novel-variant cert with MIC data (GCP)**: Rung 3 certifies spatial relevance via
  the known panel; Rung 4 would certify phenotype prediction via a held-out cohort of
  novel variants with measured susceptibility.
