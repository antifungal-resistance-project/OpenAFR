# Resistance Weather feed spec

Two machine-readable files are published next to the weather page, built from the same inputs
as the HTML (`openafr/weather_feed.py`), so they cannot disagree with the page. Research use
only; detection-only; not a clinical claim.

- `data.json` — full state (this document)
- `alerts.atom` — Atom 1.0, one entry per current alert

Contract is pinned by `tests/test_weather_feed.py`. Adding a key is backward compatible.
Renaming, removing, or changing the meaning of a key bumps `schema_version` (currently `1`).

## data.json

| key | type | meaning |
|---|---|---|
| `schema_version` | int | `1` |
| `generated` | string | UTC `YYYY-MM-DDTHH:MM:SSZ` build time |
| `notice` | string | research-use disclaimer |
| `gene` | string \| null | `FKS1` for the echinocandin arm |
| `page_url` | string \| null | canonical page URL |
| `as_of` | string \| null | as-of date of the data behind the page |
| `window_days` | int \| null | "recent" window used for emergence |
| `thresholds` | object | alert thresholds used |
| `coverage` | object | genotyping coverage (resolved isolates over total) |
| `alerts` | array | current alerts (below) |
| `release_trend` | array | one entry per committed re-call release (below) |
| `geography` | object | `{country: panel-carrier count}`, current release only, sorted by country |
| `vus_census` | array | per-variant census of non-panel hot-spot variants of uncertain significance (below) |
| `calls_provenance` | object | which committed call table(s) the page is built on |
| `certification` | object \| null | caller certification (benchmark, concordance) |

### alerts[]

`mutation`, `priority` (`watch` / `context` / `act now`), `headline`, `first_seen`, `regions[]`,
`n_carriers_total`, and `emergence` with: `first_appearance`, `rising`, `possible_backlog`,
`recent_count`, `recent_called`, `baseline_count`, `baseline_called`, `recent_prop`,
`baseline_prop`, `delta`. Optional: `position_context`, `position_note` (geometry-free
hot-spot context for novel variants; not a structural claim).

### release_trend[]

`release_tag`, `as_of`, `n_called`, `n_resistant`, `n_vus`, `prevalence`.

**How to read "resistant":** `n_resistant` counts isolates carrying a *known-resistance panel*
FKS1 substitution. `n_vus` counts resolved isolates with a hot-spot variant *not* on the panel
(uncertain significance); they are tracked but are **not** resistance. `prevalence =
n_resistant / n_called`, where `n_called` is resolved isolates. Do not sum `n_resistant` and
`n_vus` to get "resistant".

### vus_census[]

The per-variant breakdown behind `release_trend[].n_vus`. One entry per distinct non-panel
hot-spot FKS1 substitution in the latest committed call table, descending by `carriers`:
`token`, `carriers` (isolates carrying it), `position_context`, `position_note`. The last two
are the geometry-free hot-spot context (openafr/fks1_position): which hot-spot window, offset
from the anchor, and tier (`novel_substitution_at_panel_residue` / `in_hotspot_adjacent_to_panel`
/ `in_hotspot_background` / `known_low_ppv_marker`). These are **variants of uncertain
significance** — tracked, **not** counted as known resistance — and the position fields are
geography only, never a drug-fit or phenotype claim.

## alerts.atom

Title `C. auris Resistance Weather Report`; `<id>` `urn:openafr:weather:fks1`. Each entry has
title `PRIORITY: mutation`, a `<category term="priority">`, and a summary with recent counts and
regions. Entry `<id>` is stable per (mutation, priority): readers mark it read once, and a tier
change (e.g. WATCH to ACT NOW) appears as a new entry. No alerts yields a valid empty feed.
