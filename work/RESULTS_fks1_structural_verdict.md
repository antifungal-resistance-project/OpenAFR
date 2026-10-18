# Results — FKS1 structural contact label cert (moonshot #4, Rung 3)

> **Pre-registered bar:** ≥ 80% of core panel positions in `drug_contact` OR `pocket_adjacent`  
> **Pre-registration SHA-256:** `a6b26c5f0483aba6ce378e26bf8d7f54b68979953a3a45042119c41d5a3759e7`  
> **Contact map SHA-256:** `6a2033ac5bd28b821e5cc3a47b759932c5b728c7f85b69958a82f1f89902acea`

## Core panel position tiers

| Ca position | WT | Proximity tier | Distance | In contact zone | Matches prediction |
|---|---|---|---|---|---|
| 635 | F | `drug_contact` | 2.94 Å | YES | ✓ |
| 639 | S | `pocket_adjacent` | 7.16 Å | YES | ✓ |
| 691 | W | `drug_contact` | 3.66 Å | YES | ✓ |
| 1354 | R | `pocket_adjacent` | 7.53 Å | YES | ✓ |

## Secondary observation — low-PPV positions (exploratory, not in bar)

| Ca position | WT | Proximity tier | Distance |
|---|---|---|---|
| 642 | D | `pocket_proximate` | 11.34 Å |
| 690 | M | `pocket_proximate` | 9.46 Å |

## Summary

| Metric | Value |
|---|---|
| Core panel positions checked | 4 |
| Positions with contact-map entry | 4 |
| Positions in drug\_contact or pocket\_adjacent | 4 |
| Fraction | 100.0% (4/4) |
| Pre-registered bar | ≥ 80% |

## VERDICT: PASS

Core resistance panel positions cluster in the drug-contact / adjacent zone of the 9WZU cryo-EM structure. The contact map is structurally relevant to the known-resistance panel, anti-circularly validated.

**Consequence (per preregistration):** implement `structural_verdict()` and wire into the pipeline.
