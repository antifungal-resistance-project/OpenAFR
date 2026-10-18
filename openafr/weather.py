"""Resistance weather report -- an HTML view over the early-warning delivery output.

Pipeline stage: a SECOND renderer alongside delivery/alert, not a new pipeline. The
composed alert result (openafr.alert.compose_alerts) and the committed delivery state
(openafr.delivery) are the ONLY inputs; this module turns that same self-describing dict
into one public page a human can open. It re-derives nothing and pulls nothing -- if the
page and the alerts ever disagreed it would be a bug, so they read one source.

Why a page, next to the digest and the issue
---------------------------------------------
delivery.py already commits a Markdown digest and files an issue on new news, and stays
silent otherwise (the quiet no-op is the design). Neither is a thing you can send someone
a link to. This is: a "resistance weather report" at one URL that shows the current
picture on every scheduled run -- loud on a hit, honest and calm when quiet.

Three honest states (the page always shows exactly one)
-------------------------------------------------------
  * DAY-0 / WATCHING  -- no alerts over thresholds (today's real state: erg11 calls are
    pending the re-caller, so compose returns a `note` and zero alerts). The page says so
    plainly: "watching since X, nothing over threshold yet", never a false all-clear.
  * QUIET             -- same shape as day-0 but with a prior first-seen history; framed
    as "N days quiet, last new variant was ...".
  * HIT               -- one or more alerts. The page leads with the highest tier and
    renders one card per alert: the emergence numbers, geography, first-seen, the
    structural fit verdict, and (when a pre-rendered asset exists) the pocket image.

The pocket image is a repo asset, never rendered here
-----------------------------------------------------
PyMOL is out of the offline path (see structural.py); the wild-type-vs-mutant pocket PNG
is pre-rendered out-of-band and committed under pockets/<token>.png. The page only
REFERENCES it, with an onerror fallback so a missing asset degrades to text, never a
broken image. `pocket_asset` is the single place that names that path.
"""
import datetime
import html
import json

# Tier -> badge text. Keys match alert.PRIORITY exactly. The accent colours live once in
# CSS (`:root` tier tokens + `.tier-*` classes); a card names its tier by class, so the
# palette has a single source of truth and never gets duplicated into Python.
_TIER = {
    "act-now": "ACT-NOW",
    "watch": "WATCH",
    "context": "CONTEXT",
}
_TAXON = "Candida auris (Candidozyma auris)"


def pocket_asset(mutation):
    """The committed, pre-rendered pocket-image path for a mutation token.

    Single source of truth for where the out-of-band render commits its PNG and where the
    page looks for it. Filename-safe: only alphanumerics survive, so `F126L` -> pockets/
    F126L.png and a compound token like `TR34/L98H` -> pockets/TR34L98H.png.
    """
    safe = "".join(ch for ch in str(mutation) if ch.isalnum())
    return f"pockets/{safe}.png"


def _fmt_days_since(iso_date, now_date):
    """Whole days between an ISO date string and now; None when unparseable/empty."""
    if not iso_date:
        return None
    head = str(iso_date)[:10]
    try:
        then = datetime.date.fromisoformat(head)
    except ValueError:
        return None
    return (now_date - then).days


def _status_line(alerts, watching_since, now_date):
    """The one-line headline that names which of the three states the page is in."""
    if alerts:
        counts = {t: sum(1 for a in alerts if a["priority"] == t) for t in _TIER}
        parts = [f"{counts[t]} {_TIER[t].lower()}" for t in _TIER if counts[t]]
        return "RESISTANCE SEEN — " + ", ".join(parts)
    days = _fmt_days_since(watching_since, now_date)
    if days is None:
        return "WATCHING — no resistance variant over threshold yet"
    return f"QUIET — {days} day(s) with nothing new over threshold"


def _alert_card(a):
    """One alert -> an HTML card. Reads only fields the composer guarantees.

    Gene-agnostic: the ERG11 composer attaches a `structural` fit verdict and the page
    shows that row plus the pre-rendered pocket image; the FKS1 (echinocandin) composer
    is detection-only and attaches no `structural` block, so the card drops the fit row
    and the pocket image (there is no FKS1 pocket render) -- the detection-only caveat
    rides in on the composed `headline` instead.
    """
    tier_txt = _TIER.get(a["priority"], a["priority"].upper())
    tier_cls = f"tier-{a['priority']}"
    e = a["emergence"]
    s = a.get("structural")
    mut = html.escape(str(a["mutation"]))
    regions = ", ".join(a.get("regions") or []) or "undisclosed"
    rows = [
        ("Emergence",
         f"recent {e['recent_count']}/{e['recent_called']} "
         f"(&Delta; {e['delta']:+.1%})"),
        ("Backlog",
         "&#9888; possible archival" if e.get("possible_backlog") else "genuine spread"),
        ("Geography",
         f"{a.get('n_regions', 0)} region(s): {html.escape(regions)} "
         f"({a.get('n_carriers_total', 0)} carrier(s))"),
        ("First seen", html.escape(str(a.get("first_seen") or "unknown"))),
    ]
    if s:
        rows.append(
            ("Structural fit",
             f"[{html.escape(str(s['evidence_class']))}, conf={html.escape(str(s['confidence']))}] "
             f"{html.escape(str(s['fluconazole_fit_verdict']))}"))
    # FKS1 novel-variant position context (geometry-free; never a fit verdict). Present only
    # when the detection-only composer attached one for an uncharacterised hot-spot variant.
    # The composer pre-renders the sentence into `position_note`, so the page stays decoupled
    # from fks1_position (same discipline as it reading `headline`, not recomputing it).
    ctx = a.get("position_context")
    if ctx:
        rows.append(
            ("Position",
             f"[{html.escape(str(ctx['tier']))}] "
             f"{html.escape(str(a.get('position_note') or ''))}"))
    prox = ctx.get("proximity") if ctx else None
    if prox:
        tier_label = prox["proximity_tier"].replace("_", " ")
        dist = prox.get("min_dist_to_drug_A")
        dist_str = f"{dist:.1f}&thinsp;Å" if dist else "N/A"
        rows.append(
            ("Proximity",
             f"[{html.escape(tier_label)}] "
             f"{dist_str} from enfumafungin in {html.escape(prox.get('reference_structure',''))}"
             f" — {html.escape(prox.get('note','')[:120])}"))
    rows_html = "\n".join(
        f'      <div class="k">{k}</div><div class="v">{v}</div>' for k, v in rows)
    headline = html.escape(str(a.get("headline") or ""))
    img_html = ""
    if s:
        img = pocket_asset(a["mutation"])
        img_html = (
            f'\n      <img class="pocket" src="{img}" '
            f'alt="Mutant vs wild-type CYP51 pocket for {mut}"\n'
            f"           loading=\"lazy\" onerror=\"this.style.display='none'\">")
    return f"""    <article class="card {tier_cls}">
      <header><span class="badge">{tier_txt}</span>
        <code class="mut">{mut}</code></header>
      <p class="headline">{headline}</p>
      <div class="grid">
{rows_html}
      </div>{img_html}
    </article>"""


def _calm_copy(gene):
    """The two calm-state sentences, keyed to the gene the page is tracking.

    Default is the echinocandin (FKS1) surveillance product; ERG11/azole is the
    research arm, shown only when a result explicitly marks `gene == "ERG11"`.
    """
    if str(gene).upper() == "ERG11":
        return ("No azole-resistance variant is over threshold right now.",
                "When a resistance-conferring ERG11 substitution emerges, it appears "
                "here with its structural fit verdict.")
    return ("No echinocandin-resistance variant is over threshold right now.",
            "When an echinocandin-resistance FKS1 substitution emerges, it appears here "
            "with its emergence signal — detection only, no structural verdict.")


def _coverage_clause(coverage, provenance):
    """A freshness/coverage clause for the .meta line, or "" when there is nothing honest
    to say. For a public-health page, how much of the shown window is actually re-called —
    and how current those calls are — matters as much as the verdict itself.

    `coverage` -- the composer's coverage dict (recent_total / recent_called), or None.
    `provenance` -- {as_of, n_called} of the overlaid call table (the fill date), or None.
    Degrades honestly: zero re-called isolates in the window reads as "awaiting re-call",
    never as a confident current picture; missing data omits the clause rather than guessing.
    """
    if not coverage:
        return ""
    total = coverage.get("recent_total") or 0
    called = coverage.get("recent_called") or 0
    if total == 0:
        # No isolates in the window at all -- a genuinely quiet window, not a stale fill.
        # The calm "watching since" copy carries this; a coverage % would be meaningless.
        return ""
    through = (provenance or {}).get("as_of")
    through_html = (f" &middot; calls current through <code>{html.escape(str(through))}</code>"
                    if through else "")
    if called == 0:
        return f'<span class="stale">echinocandin calls: awaiting re-call{through_html}</span>'
    pct = round(100 * called / total)
    return (f"echinocandin calls: {called} of {total} window isolates re-called "
            f"({pct}%){through_html}")


def _is_awaiting(coverage):
    """True when a window exists but no isolate in it has been re-called yet -- a lapsed or
    not-yet-current fill, distinct from a genuinely never-started day-0 watch."""
    return bool(coverage) and (coverage.get("recent_total") or 0) > 0 \
        and (coverage.get("recent_called") or 0) == 0


def _pct(x):
    """A proportion (0..1) as a whole-percent string, or an em dash when absent."""
    return "—" if x is None else f"{round(100 * x)}%"


def _share_description(status, coverage, calls_provenance):
    """A one-line, plain-text summary for the <meta description>/OpenGraph card.

    Built from the same state the page shows -- the status headline plus, when a fill
    covers the window, how much of it is re-called -- so a shared link previews the real
    picture, not boilerplate. Never invents counts: falls back to the plain product line.
    """
    base = ("Echinocandin-resistance (FKS1) surveillance of Candida auris genomes from "
            "NCBI Pathogen Detection, genotyped by a certified caller and updated monthly.")
    total = (coverage or {}).get("recent_total") or 0
    called = (coverage or {}).get("recent_called") or 0
    if total and called:
        return (f"{status}. {called} of {total} recent isolates genotyped. {base}")
    return f"{status}. {base}"


def _cert_badge(certification):
    """A small 'certified caller' pill, or "" unless the concordance verdict is PASS.

    The product's core credibility claim -- the numbers come from a caller that passed a
    frozen concordance benchmark -- which was otherwise invisible on the page. Shown only
    on a genuine PASS; any other verdict (or no record) makes no claim at all.
    """
    if not certification or certification.get("verdict") != "PASS":
        return ""
    sens = _pct(certification.get("sensitivity"))
    spec = _pct(certification.get("specificity"))
    n = certification.get("n")
    return (f'  <p class="cert"><a href="#methods">&#10003; Calls from a certified FKS1 '
            f"caller &mdash; {sens} sensitivity / {spec} specificity on a "
            f"{html.escape(str(n))}-isolate benchmark</a></p>\n")


def _methods_section(certification, calls_provenance, window):
    """The collapsible 'how this works / methods' block, or "" when uncertified.

    Surfaces provenance the page otherwise dropped: the NCBI release the fill targeted, the
    re-call batch size, the window + cadence, and the full certification table scored from
    the committed concordance runlog. Gated on a certification record so it rides with the
    credibility claim; every measured number is passed in, never recomputed here.
    """
    if not certification:
        return ""
    prov = calls_provenance or {}
    release = html.escape(str(prov.get("release_tag") or "the latest NCBI release"))
    n_called = prov.get("n_called")
    fill_line = (f" The most recent batch re-called <code>{html.escape(str(n_called))}</code> "
                 f"isolates over release <code>{release}</code>."
                 if n_called else f" Calls are drawn from release <code>{release}</code>.")
    verdict = html.escape(str(certification.get("verdict")))
    bench = html.escape(str(certification.get("benchmark") or "a concordance panel"))
    n = html.escape(str(certification.get("n")))
    n_res = html.escape(str(certification.get("n_resolved")))
    sens, spec, ident = (_pct(certification.get("sensitivity")),
                         _pct(certification.get("specificity")),
                         _pct(certification.get("identity")))
    return f"""  <details class="about" id="methods">
    <summary>How this works &middot; methods &amp; provenance</summary>
    <p><b>Source.</b> {html.escape(_TAXON)} whole-genome deposits from NCBI Pathogen
      Detection, over a trailing <code>{html.escape(str(window))}</code>-day window.{fill_line}
      Public assemblies carry no echinocandin call (NCBI runs no AMR pipeline), so an
      out-of-band FKS1 re-call batch (cadence set by the ~97-day median deposit lag)
      produces the calls this page overlays.</p>
    <p><b>Certified caller.</b> The FKS1 re-caller is graded on a frozen external benchmark
      (<code>{bench}</code>, PMC12323592) on every change, through the same scorer this
      verdict uses:</p>
    <div class="grid cert-grid">
      <div class="k">Verdict</div><div class="v">{verdict} ({n_res}/{n} isolates resolved)</div>
      <div class="k">Sensitivity</div><div class="v">{sens}</div>
      <div class="k">Specificity</div><div class="v">{spec}</div>
      <div class="k">Exact-token identity</div><div class="v">{ident}</div>
    </div>
    <p class="note">Detection only &mdash; the echinocandin arm reports whether a known
      resistance variant is present, not a structural fit verdict. Research use only: a
      triage aid, not a clinical claim. Priority tiers and every number come straight from
      <code>openafr.alert</code>; this page renders that verdict, it never invents one.</p>
  </details>
"""


def _trend_section(release_trend):
    """The resistance-trend block: prevalence per committed re-call batch, with the direction
    of the latest move, or "" when there is no trajectory to show (fewer than two batches).

    `release_trend` -- the per-release series from render_weather (`as_of`, resolved
    `n_called`, `n_resistant` = KNOWN-resistance panel carriers, `n_vus` = non-panel
    variant-of-uncertain-significance carriers, `prevalence`), ascending by release.
    Reproducible from the committed call tables; passed in, never recomputed here. Degrades
    honestly: a single batch is a point, not a trend, so the section is omitted until a second
    batch gives it meaning.
    """
    if not release_trend or len(release_trend) < 2:
        return ""

    def _prev(x):
        return "&mdash;" if x is None else f"{100 * x:.1f}%"

    rows = "\n".join(
        f'      <div class="k"><code>{html.escape(str(e.get("as_of") or "—"))}</code></div>'
        f'<div class="v">{e.get("n_resistant", 0)} of {e.get("n_called", 0)} resolved '
        f'&middot; {_prev(e.get("prevalence"))}</div>'
        for e in release_trend)

    pts = [e["prevalence"] for e in release_trend if e.get("prevalence") is not None]
    delta_html = ""
    if len(pts) >= 2:
        dp = (pts[-1] - pts[-2]) * 100   # percentage points, latest vs previous batch
        if dp > 0.05:
            cls, line = "trend-up", f"&#9650; prevalence up {dp:.1f} pts since the previous batch."
        elif dp < -0.05:
            cls, line = "trend-down", (f"&#9660; prevalence down {abs(dp):.1f} pts since the "
                                       "previous batch.")
        else:
            cls, line = "trend-flat", "&#9644; prevalence held steady since the previous batch."
        delta_html = f'    <p class="trend-delta {cls}">{line}</p>\n'

    # Non-panel hot-spot variants (uncertain significance) in the latest batch -- surfaced so
    # the reader sees they are tracked but deliberately NOT counted as known resistance.
    n_vus = release_trend[-1].get("n_vus", 0)
    vus_html = (f" A further {n_vus} resolved isolate(s) carry a non-panel FKS1 variant of "
                "uncertain significance, tracked but not counted as known resistance."
                ) if n_vus else ""

    return f"""  <section class="trend">
    <h2>Resistance trend</h2>
    <div class="grid">
{rows}
    </div>
{delta_html}    <p class="note">Carriers of a <strong>known</strong> echinocandin-resistance FKS1
      substitution as a share of resolved isolates, per committed re-call batch &mdash;
      reproducible from the committed call tables.{vus_html}</p>
  </section>
"""


def _forecast_section(forecast):
    """The prevalence trajectory + near-term forecast block, or "" when there is no trajectory
    to show (fewer than two deposit-month bins).

    `forecast` -- the dict render_weather assembles from openafr.forecast: `series` (ascending
    per-deposit-month {month, n_resolved, n_resistant, prevalence}), `skill` (did the
    walk-forward backtest certify the model beats both naive baselines?), and `projection` (the
    one-step-ahead forecast_next dict, or None). Passed in, never recomputed here -- the same
    decoupling the trend/VUS sections use.

    Honest degradation is the whole point: the projected number renders ONLY when `skill` is
    True. With no demonstrated skill the section shows the monthly trajectory alone and says so
    plainly, naming the flat (pooled) prevalence the model could not beat -- never a forecast
    the backtest did not earn.
    """
    if not forecast:
        return ""
    series = forecast.get("series") or []
    if len(series) < 2:
        return ""

    def _prev(x):
        return "&mdash;" if x is None else f"{100 * x:.1f}%"

    rows = "\n".join(
        f'      <div class="k"><code>{html.escape(str(e.get("month") or "—"))}</code></div>'
        f'<div class="v">{e.get("n_resistant", 0)} of {e.get("n_resolved", 0)} resolved '
        f'&middot; {_prev(e.get("prevalence"))}</div>'
        for e in series)

    proj = forecast.get("projection") or {}
    if forecast.get("skill") and proj.get("point") is not None:
        last = series[-1].get("prevalence")
        point = proj["point"]
        if last is not None and point - last > 0.001:
            cls, arrow = "trend-up", "&#9650;"
        elif last is not None and last - point > 0.001:
            cls, arrow = "trend-down", "&#9660;"
        else:
            cls, arrow = "trend-flat", "&#9644;"
        pi = ""
        if proj.get("lo") is not None and proj.get("hi") is not None:
            pi = f" (95% PI {_prev(proj['lo'])}&ndash;{_prev(proj['hi'])})"
        fc_html = (f'    <p class="trend-delta {cls}">{arrow} Projected '
                   f'<code>{html.escape(str(proj.get("next_month") or "next month"))}</code>: '
                   f"{_prev(point)}{pi}.</p>\n")
        note = ("One-deposit-month-ahead projection, surfaced because a pre-registered "
                "walk-forward showed the trend model out-predicts naive persistence "
                "(work/RESULTS_fks1_forecast.md). Deposit-month axis, not epidemiological "
                "incidence.")
    else:
        pooled = proj.get("naive_pooled")
        pooled_str = _prev(pooled) if pooled is not None else "the series mean"
        fc_html = (f'    <p class="trend-delta trend-flat">&#9644; No forecast skill '
                   f"demonstrated over a flat baseline at this series depth; showing the "
                   f"trajectory only (prevalence is roughly stationary near "
                   f"{pooled_str}).</p>\n")
        note = ("A pre-registered walk-forward did not show the trend model beating a flat "
                "baseline on held-out months (work/RESULTS_fks1_forecast.md), so no projection "
                "is claimed. Deposit-month axis, not epidemiological incidence.")

    return f"""  <section class="trend forecast">
    <h2>Prevalence trajectory &amp; forecast</h2>
    <div class="grid">
{rows}
    </div>
{fc_html}    <p class="note">Panel-gated known-resistance carriers as a share of resolved
      isolates, per NCBI deposit month &mdash; reproducible from the committed call table.
      {html.escape(note)}</p>
  </section>
"""


def _geo_section(geo_breakdown):
    """Where the current release's resistant carriers are: a country-ranked list, or "" when
    there is nothing located. Labelled 'current release' because geography comes from the live
    NCBI pull, not the committed historical series (the call table does not yet carry country).

    `geo_breakdown` -- {country: carrier_count} from render_weather; passed in, not recomputed.
    """
    if not geo_breakdown:
        return ""
    items = sorted(geo_breakdown.items(), key=lambda kv: (-kv[1], kv[0]))
    rows = "\n".join(
        f'      <div class="k">{html.escape(str(country))}</div><div class="v">{n}</div>'
        for country, n in items)
    return f"""  <section class="geo">
    <h2>Where carriers are &middot; current release</h2>
    <div class="grid">
{rows}
    </div>
    <p class="note">Resistant FKS1 carriers by country in the current NCBI release's metadata
      (not the committed historical series). Isolates with no country are omitted.</p>
  </section>
"""


def _vus_section(vus_census):
    """The 'variants of uncertain significance' section: the non-panel hot-spot FKS1 variants
    the current table carries, each with how many isolates carry it and WHERE it sits relative
    to the known resistance geography -- or "" when there are none (honest degradation, same as
    trend/geo).

    `vus_census` -- the enriched per-variant list from render_weather: each entry carries
    `token`, `carriers`, `position_context` (geometry-free, openafr/fks1_position) and a
    pre-rendered `position_note` sentence. Passed in, never recomputed here. This turns the
    trend section's bare "a further N carry a variant of uncertain significance" count into the
    per-variant detail behind it -- strictly geography, never a fit or phenotype claim.
    """
    if not vus_census:
        return ""
    def _vus_row(v):
        token = html.escape(str(v["token"]))
        ctx = v.get("position_context") or {}
        tier_str = html.escape(str(ctx.get("tier", "")))
        note = html.escape(str(v.get("position_note") or ""))
        carriers = v.get("carriers", 0)
        prox = ctx.get("proximity")
        prox_html = ""
        if prox:
            tier_label = html.escape(prox["proximity_tier"].replace("_", " "))
            dist = prox.get("min_dist_to_drug_A")
            dist_str = f"{dist:.1f}&thinsp;Å" if dist else "N/A"
            prox_html = (
                f' <span class="prox-tier">[proximity: {tier_label}, '
                f'{dist_str} from drug]</span>'
            )
        return (
            f'      <div class="k"><code>{token}</code> '
            f'<span class="vus-tier">{tier_str}</span>{prox_html}</div>'
            f'<div class="v">{carriers} carrier(s) &mdash; {note}</div>'
        )

    rows = "\n".join(_vus_row(v) for v in vus_census)
    return f"""  <section class="vus">
    <h2>Variants of uncertain significance</h2>
    <div class="grid">
{rows}
    </div>
    <p class="note">Non-panel FKS1 hot-spot substitutions seen in the committed calls:
      <strong>tracked but NOT counted as known resistance</strong>. The position line states
      only where each residue sits relative to the known resistance geography (which hot-spot,
      how far from the anchor); the proximity line shows the minimum heavy-atom distance to
      enfumafungin (same FKS1 binding site as echinocandins) in the 9WZU cryo-EM structure.
      Neither is a drug-fit or phenotype claim. Detection only.</p>
  </section>
"""


def _structured_data(description, certification, provenance, page_url):
    """A schema.org JSON-LD <script> block describing the page + dataset, for search and
    link-preview engines. Built from the same inputs already in scope -- no I/O. Always
    emitted (the dataset description stands on its own); a PASS `certification` additionally
    records the caller's measured sensitivity/specificity as `variableMeasured`.

    The serialized JSON is escaped for safe embedding in HTML (`<`/`>`/`&`), so it can never
    break out of the <script> element even though every field here is page-controlled.
    """
    prov = provenance or {}
    dataset = {
        "@type": "Dataset",
        "name": "Candida auris FKS1 echinocandin-resistance surveillance",
        "description": description,
        "measurementTechnique":
            "Whole-genome FKS1 re-calling of NCBI Pathogen Detection deposits",
        "isBasedOn": "NCBI Pathogen Detection",
        "creator": {"@type": "Organization", "name": "OpenAFR"},
        "license": "Research use only",
    }
    if prov.get("release_tag"):
        dataset["version"] = str(prov["release_tag"])
    cert = certification or {}
    if cert.get("verdict") == "PASS":
        measured = [
            {"@type": "PropertyValue", "name": f"certified caller {label}",
             "value": round(100 * cert[key]), "unitText": "%"}
            for label, key in (("sensitivity", "sensitivity"),
                               ("specificity", "specificity"))
            if cert.get(key) is not None
        ]
        if measured:
            dataset["variableMeasured"] = measured
    doc = {
        "@context": "https://schema.org",
        "@type": "MedicalWebPage",
        "name": "C. auris Resistance Weather Report",
        "description": description,
        "about": {"@type": "MedicalCondition",
                  "name": "Echinocandin-resistant Candida auris infection"},
        "mainEntity": dataset,
    }
    if page_url:
        doc["url"] = str(page_url)
    safe = (json.dumps(doc, separators=(",", ":"))
            .replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026"))
    return f'<script type="application/ld+json">{safe}</script>\n'


def render_page(alert_result, *, now=None, watching_since=None, calls_provenance=None,
                certification=None, release_trend=None, geo_breakdown=None, vus_census=None,
                forecast=None, page_url=None, feeds=False):
    """Render the whole page as a self-contained HTML string (inline CSS, no framework).

    `alert_result` -- the dict from alert.compose_alerts (or compose_fks1_alerts).
    `now`          -- a datetime/ISO string for "last checked"; defaults to UTC now.
    `watching_since` -- ISO date the watch started / last-delivered date, for the quiet
                        framing. None renders the pure day-0 line.
    `calls_provenance` -- the overlaid FKS1 call table's INDEX row (as_of, n_called,
                        release_tag, window_days...), for the freshness clause and the
                        methods section; None omits them. Passed in (not read) to keep
                        this deterministic.
    `certification` -- the caller's concordance verdict dict (verdict, n, sensitivity,
                        specificity, identity, benchmark), scored upstream from the
                        committed runlog; a PASS adds the credibility badge + methods
                        section. None (or non-PASS) makes no certification claim.
    `release_trend` -- the per-release resistance series (as_of, n_called, n_resistant,
                        prevalence), ascending; two or more batches render the trend section.
                        None/one batch omits it.
    `geo_breakdown` -- {country: carrier_count} for the current release; renders the geography
                        section. None/empty omits it.
    `vus_census` -- the per-variant list of non-panel hot-spot variants (uncertain
                        significance), each with `token`, `carriers`, `position_context`,
                        `position_note`; renders the VUS section. None/empty omits it.
    `page_url` -- canonical URL for the shareable OpenGraph/Twitter card; None omits og:url.
    `feeds` -- True links the sibling machine-readable `data.json` / `alerts.atom` (written by
                        scripts/render_weather.py); off by default so a bare page has no dead links.
    Deterministic given its inputs -- no I/O, no network, no wall-clock unless `now` is None.
    """
    if now is None:
        now = datetime.datetime.now(datetime.timezone.utc)
    now_dt = (datetime.datetime.fromisoformat(str(now).replace("Z", "+00:00"))
              if isinstance(now, str) else now)
    now_date = now_dt.date()
    last_checked = now_dt.strftime("%Y-%m-%d %H:%M UTC")

    r = alert_result or {}
    alerts = r.get("alerts") or []
    status = _status_line(alerts, watching_since, now_date)
    # The page's headline accent is a CSS var reference (resolved from the :root tokens),
    # not a hex literal -- the top alert's tier colour on a hit, the default green otherwise.
    page_accent = (f"var(--tier-{alerts[0]['priority']})" if alerts else "var(--accent)")

    coverage = r.get("coverage")
    if alerts:
        body = "\n".join(_alert_card(a) for a in alerts)
    else:
        note = r.get("note")
        note_html = (f'<p class="note">{html.escape(str(note))}</p>' if note else "")
        calm_headline, calm_body = _calm_copy(r.get("gene"))
        if _is_awaiting(coverage):
            # A window exists but nothing in it is re-called yet: a lapsed / not-yet-current
            # fill, not a never-started watch. Frame it by the last batch, not day-0.
            through = (calls_provenance or {}).get("as_of")
            watch_line = ("Awaiting the next echinocandin re-call batch"
                          + (f"; the last batch covered through {html.escape(str(through))}."
                             if through else "."))
        else:
            since = html.escape(str(watching_since)) if watching_since else "the first run"
            watch_line = (f"Watching {html.escape(_TAXON)} genome deposits since {since}. "
                          f"{calm_body}")
        body = f"""    <article class="card calm">
      <p class="headline">{calm_headline}</p>
      <p>{watch_line}</p>
      {note_html}
    </article>"""

    as_of = html.escape(str(r.get("as_of", "—")))
    window = html.escape(str(r.get("window_days", "—")))
    coverage_clause = _coverage_clause(coverage, calls_provenance)
    # The freshness/coverage clause is the single most important fact on the meta line, so
    # it leads on its own readable line; the provenance detail (checked/as-of/window/source)
    # follows, muted. Same facts and the same honest degradation as before, better hierarchy.
    fresh_html = (f'  <p class="meta-fresh">{coverage_clause}</p>\n'
                  if coverage_clause else "")

    cert_badge = _cert_badge(certification)
    trend_html = _trend_section(release_trend)
    forecast_html = _forecast_section(forecast)
    geo_html = _geo_section(geo_breakdown)
    vus_html = _vus_section(vus_census)
    methods_html = _methods_section(certification, calls_provenance, window)
    share_desc_plain = _share_description(status, coverage, calls_provenance)
    share_desc = html.escape(share_desc_plain)
    jsonld_html = _structured_data(share_desc_plain, certification, calls_provenance, page_url)
    title = "C. auris Resistance Weather Report"
    feeds_html = ('    Machine-readable: <a href="data.json">data.json</a> &middot; '
                  '<a href="alerts.atom">alerts.atom</a> (Atom feed).\n' if feeds else "")
    feed_link_tag = ('\n<link rel="alternate" type="application/atom+xml" '
                     'title="Alerts" href="alerts.atom">' if feeds else "")
    og_url_tag = (f'\n<meta property="og:url" content="{html.escape(str(page_url))}">'
                  if page_url else "")
    canonical_tag = (f'\n<link rel="canonical" href="{html.escape(str(page_url))}">'
                     if page_url else "")
    favicon = ("data:image/svg+xml,<svg xmlns='http://www.w3.org/2000/svg' "
               "viewBox='0 0 32 32'><rect width='32' height='32' rx='7' fill='%230d1017'/>"
               "<circle cx='16' cy='16' r='5' fill='%234dd6a0'/><circle cx='16' cy='16' "
               "r='10' fill='none' stroke='%234dd6a0' stroke-opacity='.4' "
               "stroke-width='2'/></svg>")
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<link rel="icon" href="{favicon}">{canonical_tag}{feed_link_tag}
<meta name="description" content="{share_desc}">
<meta property="og:title" content="{title}">
<meta property="og:description" content="{share_desc}">
<meta property="og:type" content="website">{og_url_tag}
<meta name="twitter:card" content="summary">
<meta name="twitter:title" content="{title}">
<meta name="twitter:description" content="{share_desc}">
{jsonld_html}<style>
  :root {{
    color-scheme: dark;
    --bg:#0d1017; --surface:#151a23; --surface-accent:#132019; --text:#e6e9ef;
    --muted:#8a94a6; --muted-strong:#9aa4b2; --border:#232a36; --code:#c3cad6;
    --accent:#4dd6a0; --accent-border:#2a5a47; --link:#7aa2ff; --stale:#e0b050;
    --footer:#5c6675;
    --tier-act-now:#ff5c5c; --tier-watch:#ffcc4d; --tier-context:#8a94a6;
    --accent-page:{page_accent};
  }}
  body {{ margin:0; background:var(--bg); color:var(--text);
         font:16px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif; }}
  .wrap {{ max-width:760px; margin:0 auto; padding:2rem 1.25rem 4rem; }}
  h1 {{ font-size:1.15rem; letter-spacing:.02em; margin:0 0 .25rem;
        text-transform:uppercase; color:var(--muted-strong); font-weight:600; }}
  .status {{ font-size:1.6rem; font-weight:700; margin:.2rem 0 1rem;
             color:var(--accent-page); }}
  .meta {{ margin:0 0 2rem; }}
  .meta-fresh {{ color:var(--text); font-size:.9rem; margin:0 0 .35rem; }}
  .meta-detail {{ color:var(--muted); font-size:.82rem; margin:0; }}
  .meta code {{ color:var(--code); }}
  .meta .stale {{ color:var(--stale); }}
  .cert {{ margin:-.5rem 0 2rem; font-size:.85rem; }}
  .cert a {{ display:inline-block; color:var(--accent); text-decoration:none;
             border:1px solid var(--accent-border); background:var(--surface-accent);
             border-radius:999px; padding:.2rem .7rem; }}
  .cert a:hover {{ border-color:var(--accent); }}
  .about {{ margin-top:2.5rem; border-top:1px solid var(--border); padding-top:1rem;
            color:var(--muted); font-size:.85rem; }}
  .about summary {{ cursor:pointer; color:var(--muted-strong); font-weight:600; }}
  .about code {{ color:var(--code); }}
  .about .cert-grid {{ max-width:420px; margin:.6rem 0; }}
  h2 {{ font-size:.95rem; letter-spacing:.02em; text-transform:uppercase;
        color:var(--muted-strong); font-weight:600; margin:2.5rem 0 .6rem; }}
  .trend, .geo, .vus {{ margin:0 0 1rem; }}
  .vus-tier {{ color:var(--muted); font-size:.78rem; }}
  .trend-delta {{ font-size:.9rem; margin:.5rem 0 0; }}
  .trend-up {{ color:var(--tier-act-now); }}
  .trend-down {{ color:var(--accent); }}
  .trend-flat {{ color:var(--muted); }}
  .card {{ background:var(--surface); border:1px solid var(--border);
           border-left:4px solid var(--accent,var(--border));
           border-radius:10px; padding:1.1rem 1.25rem; margin:1rem 0; }}
  .tier-act-now {{ --accent:var(--tier-act-now); }}
  .tier-watch {{ --accent:var(--tier-watch); }}
  .tier-context {{ --accent:var(--tier-context); }}
  header {{ display:flex; align-items:center; gap:.6rem; margin-bottom:.5rem; }}
  .badge {{ font-size:.7rem; font-weight:700; letter-spacing:.05em;
            padding:.15rem .5rem; border-radius:999px; color:var(--bg);
            background:var(--accent); }}
  .mut {{ font-size:1.05rem; color:var(--text); }}
  .headline {{ margin:.3rem 0 .8rem; }}
  .grid {{ display:grid; grid-template-columns:max-content 1fr; gap:.3rem .9rem;
           font-size:.9rem; }}
  .grid .k {{ color:var(--muted); }}
  .pocket {{ display:block; width:100%; max-width:520px; margin:1rem 0 0;
             border-radius:8px; border:1px solid var(--border); }}
  .note {{ color:var(--muted); font-size:.85rem; font-style:italic; }}
  footer {{ margin-top:3rem; color:var(--footer); font-size:.78rem; }}
  a {{ color:var(--link); }}
</style>
</head>
<body>
<div class="wrap">
  <h1>C. auris Resistance Weather Report</h1>
  <p class="status">{html.escape(status)}</p>
  <div class="meta">
{fresh_html}    <p class="meta-detail">Last checked <code>{last_checked}</code> &middot;
       as of <code>{as_of}</code> &middot; window <code>{window}</code>d &middot;
       source: NCBI Pathogen Detection via OpenAFR early-warning.</p>
  </div>
{cert_badge}{body}
{trend_html}{forecast_html}{geo_html}{vus_html}{methods_html}  <footer>
    OpenAFR early-warning. Research use only — a triage aid, not a clinical claim.
    Echinocandin calls come from a monthly FKS1 re-call batch over the live NCBI feed
    (cadence set by the ~97-day median deposit lag); between batches the coverage line above
    states how current they are. Priority tiers and every number come straight from
    <code>openafr.alert</code>; this page renders that verdict, it never invents one.
{feeds_html}  </footer>
</div>
</body>
</html>
"""
