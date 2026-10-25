# Resistance Weather feed spec

Two machine-readable files are published next to the weather page, built from the same inputs
as the HTML (`openafr/weather_feed.py`), so they cannot disagree with the page. Research use
only; detection-only; not a clinical claim.

- `data.json` — full state (this document)
- `alerts.atom` — Atom 1.0, one entry per current alert

Contract is pinned by `tests/test_weather_feed.py` and re-checked at publish time by
`scripts/check_published_page.py`. Adding a key is backward compatible. Renaming, removing, or
changing the meaning of a key bumps `schema_version` (currently `5`).

**Changelog.** `5` (2026-10): added `scope_limits` — the claims this service deliberately does
**not** make yet (chiefly a per-isolate calibrated probability), each with its reason and the
condition that would unlock it. `4` (2026-10): added `geography_prevalence` — panel-gated
prevalence (resistant / resolved, with a Wilson CI) per place plus the national headline rate,
so a consumer can see *where* resistance is proportionally elevated, not just where raw carriers
are. `3` (2026-10): added `panel_census` — the per-marker breakdown of the known-resistance
panel (which substitution each resistant carrier has) behind the single headline prevalence
number. `2` (2026-10): added `fill_age_days` + `fill_stale` so a consumer can read the overlaid
fill's real age — `generated` advances on every weekly publish even when the underlying calls
are frozen, so it is **not** a freshness signal on its own.

## data.json

| key | type | meaning |
|---|---|---|
| `schema_version` | int | `5` |
| `generated` | string | UTC `YYYY-MM-DDTHH:MM:SSZ` build time — advances on every weekly publish; **not** a data-freshness signal (see `fill_age_days`) |
| `notice` | string | research-use disclaimer |
| `gene` | string \| null | `FKS1` for the echinocandin arm |
| `page_url` | string \| null | canonical page URL |
| `as_of` | string \| null | as-of date of the data behind the page |
| `window_days` | int \| null | "recent" window used for emergence |
| `fill_age_days` | int \| null | age of the overlaid fill at build time (`generated` − fill `as_of`); `null` when unknown |
| `fill_stale` | bool | true once the fill is older than one `window_days` (the sampled window is behind the live feed; still older ⇒ aged fully out) |
| `thresholds` | object | alert thresholds used |
| `coverage` | object | genotyping coverage (resolved isolates over total) |
| `alerts` | array | current alerts (below) |
| `release_trend` | array | one entry per committed re-call release (below) |
| `forecast` | object \| null | prevalence trajectory + near-term forecast (below) |
| `geography` | object | `{country: panel-carrier count}`, current release only, sorted by country |
| `geography_history` | array | per-release `{country: panel-carrier count}` over the committed series (below) |
| `panel_census` | array | per-marker breakdown of the known-resistance panel (below) |
| `geography_prevalence` | object \| null | panel-gated prevalence per place + the national headline rate (below) |
| `vus_census` | array | per-variant census of non-panel hot-spot variants of uncertain significance (below) |
| `scope_limits` | array | claims the service deliberately does **not** make yet, with reasons (below) |
| `calls_provenance` | object | which committed call table(s) the page is built on |
| `certification` | object \| null | caller certification (benchmark, concordance) |

### alerts[]

`mutation`, `priority` (`watch` / `context` / `act now`), `headline`, `first_seen`, `regions[]`,
`n_carriers_total`, and `emergence` with: `first_appearance`, `rising`, `possible_backlog`,
`recent_count`, `recent_called`, `baseline_count`, `baseline_called`, `recent_prop`,
`baseline_prop`, `delta`. Optional: `position_context`, `position_note` (hot-spot context for
novel variants — where in the resistance geography this variant lands, plus structural proximity
when the contact map is available; not a resistance prediction, not a fit verdict).

**`position_context` sub-fields** (moonshot #4, Rungs 1–3; all fields derived, never invented):

| sub-field | type | meaning |
|---|---|---|
| `tier` | string | `novel_substitution_at_panel_residue` / `in_hotspot_adjacent_to_panel` / `in_hotspot_background` / `known_low_ppv_marker` |
| `window` | string | `HS1` / `HS2` / `HS3` |
| `anchor` | string | window anchor as `<wt><pos>` (e.g. `S639`) |
| `anchor_distance` | int | signed residue offset from anchor |
| `at_panel_residue` | bool | true if any panel letter falls at this residue |
| `panel_letters_at_position` | array | known panel mutant letters at this residue (empty if none) |
| `nearest_core_marker` | string \| null | position of the nearest core marker (e.g. `"639"`) |
| `distance_to_nearest` | int \| null | residue distance to nearest core marker |
| `proximity` | object \| absent | **(Rung 2)** structural proximity from 9WZU cryo-EM; present when the contact map covers this position. Keys: `proximity_tier` (`drug_contact` / `pocket_adjacent` / `pocket_proximate` / `distal`), `min_dist_to_drug_A` (float Å), `position`, `window`, `sc_position`, `sc_aa`, `reference_structure`, `note` |
| `structural_verdict` | object \| absent | **(Rung 3)** spatial co-location label; present only for `drug_contact` or `pocket_adjacent` positions. Keys: `evidence_class` (`"proximity_label"` — **NOT** `"docking"` or `"measured"`; must not be read as ACT-NOW), `proximity_tier`, `min_dist_to_drug_A`, `label`, `note` |

The `evidence_class: "proximity_label"` constraint is enforced: `structural_verdict` cannot trigger ACT-NOW (that tier requires `"docking"` or `"measured"`). This is a spatial co-location label, not a phenotype prediction.

### release_trend[]

`release_tag`, `as_of`, `n_called`, `n_resistant`, `n_vus`, `prevalence`.

**How to read "resistant":** `n_resistant` counts isolates carrying a *known-resistance panel*
FKS1 substitution. `n_vus` counts resolved isolates with a hot-spot variant *not* on the panel
(uncertain significance); they are tracked but are **not** resistance. `prevalence =
n_resistant / n_called`, where `n_called` is resolved isolates. Do not sum `n_resistant` and
`n_vus` to get "resistant".

### forecast

`series[]` — one entry per NCBI **deposit month**, ascending: `month` (`YYYY-MM`), `n_resolved`,
`n_resistant` (panel-gated), `prevalence`. `skill` (bool) — whether a pre-registered walk-forward
showed the trend model beats both naive baselines (`work/RESULTS_fks1_forecast.md`). `n_folds`
(int) — walk-forward folds evaluated. `model_mae`, `naive_last_mae`, `naive_pooled_mae` (float) —
held-out mean absolute error of the trend model and the two naive baselines the skill test
compares against. `projection` — the one-deposit-month-ahead `forecast_next` dict (`next_month`,
`point`, `lo`/`hi` 95% PI, `naive_last`, `naive_pooled`, `slope_logit`, `n_months`,
`has_interval`), or `null`.

**How to read it:** trust `projection.point`/`lo`/`hi` as a *forecast* **only** when `skill` is
`true`. When `skill` is `false` the series is still a real trajectory, but no projection is
claimed (prevalence is roughly stationary near `projection.naive_pooled`). The axis is **deposit
month, not epidemiological incidence** — deposit lags collection ~97 days and tracks sequencing
volume, not transmission. Detection-only, research use.

### geography_history[]

One entry per committed re-call release that carries geography, oldest release first:
`release_tag`, `as_of`, and `geography` (`{country: panel-carrier count}`, sorted by country,
same panel gate as the top-level `geography`). `geography` is **country-level** (as NCBI reported
it — no finer region), isolates with no recorded country are omitted, and a release whose call
table predates the geography column is absent entirely (not a zero). Use this to read the
geographic picture over time; the top-level `geography` is the current live pull only.

### panel_census[]

The per-marker breakdown behind `release_trend[].n_resistant` — the single headline prevalence
number split by which known-resistance substitution each carrier has. One entry per distinct
panel token in the latest committed call table, descending by `carriers`: `token`, `carriers`
(isolates carrying it), and `low_ppv` (true for the known-but-low-PPV positions D642Y / M690I,
which are counted in the headline but split R/S rather than tracking resistance). Per-token
counts sum to `n_resistant` except where one isolate carries more than one panel marker (then it
is counted under each). Detection only — the marker and its count, never a per-isolate phenotype.

### geography_prevalence

Panel-gated echinocandin **prevalence** by place for the current release (the proportional
counterpart to the raw-count `geography`), or `null` (azole arm, a snapshot render, a missing
table, or no placeable resolved isolates). Keys:

- `places[]` — one entry per place with at least `min_resolved` (5) resolved isolates, sorted by
  prevalence (then CI lower bound) descending: `place`, `n_resolved` (the denominator),
  `n_panel_positive` (the numerator), `event_frequency` (= positive / resolved, 0..1), `ci95`
  (Wilson 95% CI `[lo, hi]`).
- `pooled_small` — places below the `min_resolved` floor, pooled so the footnote math holds:
  place count + their summed `n_resolved` / `n_panel_positive`.
- `headline_prevalence` — the national panel-gated rate (`backtest.fks1_panel_prevalence`), the
  same single source of truth as the headline, so "above the national rate" never drifts.

The summed `n_panel_positive` across `places` + `pooled_small` equals the national count. Read a
place as *proportionally elevated* only when its `event_frequency` stands out, not its raw volume.
Detection only, research use.

### scope_limits[]

The claims this service deliberately does **not** make yet — surfaced so a consumer sees a
considered boundary, not a silent gap. Fixed product policy (`openafr.weather.SCOPE_LIMITS`),
shared with the HTML page so the two cannot disagree; the reasons are the honest ceilings fixed
in advance by this project's frozen pre-registrations / contract, not claims. One entry per
deferral, with keys:

| key | type | meaning |
|---|---|---|
| `claim` | string | the claim being withheld (e.g. "a per-isolate calibrated resistance probability") |
| `status` | string | `deferred` (earnable later) or `out of scope` (not this service's job) |
| `reason` | string | plain-language why it is withheld today |
| `unlock` | string \| null | the condition that would earn it, or `null` when out of scope |
| `source` | string | the frozen doc fixing the reason (e.g. `work/PREREGISTRATION_calibration.md`) |

The lead entry is the **calibrated-probability** deferral: the feed reports presence of a known
resistance marker (categorical), never a probability, because no public panel pairs genomes with
lab susceptibility at the size needed to *measure* that a stamped probability is reliable
(`work/PREREGISTRATION_calibration.md`, `docs/DIAGNOSTIC_VERDICT_CONTRACT.md`).

### vus_census[]

The per-variant breakdown behind `release_trend[].n_vus`. One entry per distinct non-panel
hot-spot FKS1 substitution in the latest committed call table, descending by `carriers`:
`token`, `carriers` (isolates carrying it), `position_context`, `position_note`. `position_context`
carries the same sub-fields documented under `alerts[]` above (Rungs 1–3 when available). These
are **variants of uncertain significance** — tracked, **not** counted as known resistance — and
the position fields are spatial context only, never a drug-fit or phenotype claim.

## alerts.atom

Title `C. auris Resistance Weather Report`; `<id>` `urn:openafr:weather:fks1`. Each entry has
title `PRIORITY: mutation`, a `<category term="priority">`, and a summary with recent counts and
regions. Entry `<id>` is stable per (mutation, priority): readers mark it read once, and a tier
change (e.g. WATCH to ACT NOW) appears as a new entry. No alerts yields a valid empty feed.
