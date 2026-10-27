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


def _fmt_first_seen(fs):
    """Human-readable 'First seen' value from the {collected, visible} dict.

    collected keeps the raw NCBI granularity (bare year / year-month / ISO, on purpose);
    `visible` is when the isolate surfaced in NCBI, and the gap between the two is the
    backlog tell -- so we show both rather than collapse them. Falls back gracefully when
    only one is present, or 'unknown' for anything empty/unexpected.
    """
    if not isinstance(fs, dict):
        return str(fs) if fs else "unknown"
    collected = (fs.get("collected") or "").strip()
    visible = (fs.get("visible") or "").strip()
    if collected and visible:
        return f"{collected} (visible {visible})"
    if collected:
        return collected
    if visible:
        return f"visible {visible}"
    return "unknown"


def _status_line(alerts, watching_since, now_date):
    """The one-line headline that names which of the three states the page is in."""
    if alerts:
        counts = {t: sum(1 for a in alerts if a["priority"] == t) for t in _TIER}
        parts = [f"{counts[t]} {_TIER[t].lower()}" for t in _TIER if counts[t]]
        return "Resistance detected — " + ", ".join(parts)
    days = _fmt_days_since(watching_since, now_date)
    if days is None:
        return "Watching — no echinocandin-resistant variant has crossed the alert threshold yet"
    return f"Quiet — {_plural(days, 'day')} with nothing new above the alert threshold"


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
    # Only the ERG11 composer attaches a fluconazole fit verdict here. Guard on its shape so a
    # differently-shaped `structural` dict (e.g. a FKS1 contact label) can never crash the
    # render -- render_page runs outside render_weather's fail-open guard, so a KeyError here
    # would fail the deploy. The FKS1 Rung-3 label rides in via position_context instead.
    s = a.get("structural")
    if not (isinstance(s, dict) and "fluconazole_fit_verdict" in s):
        s = None
    mut = html.escape(str(a["mutation"]))
    regions = ", ".join(a.get("regions") or []) or "undisclosed"
    delta = e["delta"]
    delta_txt = f"{'up' if delta >= 0 else 'down'} {abs(delta):.1%}"
    rows = [
        ("Emergence",
         f"{e['recent_count']} of {_plural(e['recent_called'], 'recent isolate')} "
         f"({delta_txt})"),
        ("Backlog",
         "&#9888; possible archival" if e.get("possible_backlog") else "genuine spread"),
        ("Geography",
         f"{_plural(a.get('n_regions', 0), 'region')}: {html.escape(regions)} "
         f"({_plural(a.get('n_carriers_total', 0), 'carrier')})"),
        ("First seen", html.escape(_fmt_first_seen(a.get("first_seen")))),
    ]
    if s:
        rows.append(
            ("Structural fit",
             f"{html.escape(str(s['fluconazole_fit_verdict']))} "
             f"({html.escape(str(s['evidence_class']))}, "
             f"confidence {html.escape(str(s['confidence']))})"))
    # FKS1 novel-variant position context (geometry-free; never a fit verdict). Present only
    # when the detection-only composer attached one for an uncharacterised hot-spot variant.
    # The composer pre-renders the sentence into `position_note`, so the page stays decoupled
    # from fks1_position (same discipline as it reading `headline`, not recomputing it).
    ctx = a.get("position_context")
    if ctx:
        note = html.escape(str(a.get("position_note") or ""))
        tier = html.escape(str(ctx["tier"]))
        rows.append(
            ("Position", f"{note} ({tier})" if note else tier))
    prox = ctx.get("proximity") if ctx else None
    if prox:
        # The proximity note is already a complete, plain sentence (residue, hot-spot,
        # distance, structure); show it as-is rather than re-prefixing the distance.
        note = html.escape(str(prox.get("note") or ""))
        if note:
            note = note[0].upper() + note[1:]
        tier_label, _ = _prox_display(prox)
        rows.append(("Proximity", note or html.escape(tier_label)))
    sv = ctx.get("structural_verdict") if ctx else None
    if sv:
        label = html.escape(sv.get("label", ""))
        rows.append(
            ("Structural label",
             f"{label} &#x2014; spatial co-location only, not a resistance prediction"))
    rows_html = "\n".join(
        f"      <dt>{k}</dt><dd>{v}</dd>" for k, v in rows)
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
      <dl class="grid">
{rows_html}
      </dl>{img_html}
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
            "with its emergence signal — we flag the marker, not a structural verdict.")


def _fill_age_tier(provenance, window_days, now_date):
    """Classify the overlaid fill's age at render time -> (tier, age_days).

    The page's display window is ANCHORED to the committed fill (fill_status.display_window),
    so coverage reads a confident "100%" with a static "calls current through <fill date>"
    forever, no matter how old the fill is. That is honest about coverage but silent about age:
    a reader (or a consumer keying off the weekly-refreshed "Last checked") sees no signal that
    the underlying calls are months old. This supplies that signal.

    tier is 'fresh' while the fill is younger than one window, 'stale' once it is older than one
    window (its calls still show, but the sampled window is behind the live feed), and 'dark'
    once the fill's entire window has rolled past a trailing-today window -- i.e. age >= 2 *
    window_days, which is exactly fill_status._window_aged_out(as_of, window_days, today - window)
    (kept as local date math so this render layer stays I/O-free). Returns (None, None) when the
    fill date or window is unknown -- omit rather than guess, matching the clause's discipline.
    """
    as_of = (provenance or {}).get("as_of")
    age_days = _fmt_days_since(as_of, now_date) if as_of else None
    try:
        window = int(window_days)
    except (TypeError, ValueError):
        window = None
    if age_days is None or window is None or window <= 0:
        return None, None
    if age_days >= 2 * window:
        return "dark", age_days
    if age_days > window:
        return "stale", age_days
    return "fresh", age_days


def _age_clause(provenance, window_days, now_date):
    """The data-age half of the freshness line, or "" when age can't be stated honestly.

    Fresh fills say nothing (the "calls current through" date already carries it); a fill past
    one window gets a visible, coloured "fill is N days old" warning that strengthens once the
    window has fully aged out -- the page's own confidence degrading with its data."""
    tier, age_days = _fill_age_tier(provenance, window_days, now_date)
    if tier is None or tier == "fresh":
        return ""
    cls = "dark" if tier == "dark" else "stale"
    aged = (" &mdash; past its sampling window, refresh overdue" if tier == "dark"
            else " &mdash; older than the sampling window")
    return (f' &middot; <span class="{cls}">fill is {age_days} day'
            f'{"" if age_days == 1 else "s"} old{aged}</span>')


def _coverage_clause(coverage, provenance, now_date=None, window_days=None):
    """A freshness/coverage clause for the .meta line, or "" when there is nothing honest
    to say. For a public-health page, how much of the shown window is actually re-called —
    and how current those calls are — matters as much as the verdict itself.

    `coverage` -- the composer's coverage dict (recent_total / recent_called), or None.
    `provenance` -- {as_of, n_called} of the overlaid call table (the fill date), or None.
    `now_date`/`window_days` -- render date + displayed window, used only to append an honest
    age warning once the anchored fill is older than its window (see `_age_clause`); omitted
    when either is unknown.
    Degrades honestly: zero re-called isolates in the window reads as "awaiting re-call",
    never as a confident current picture; missing data omits the clause rather than guessing.
    """
    if not coverage:
        return ""
    total = coverage.get("recent_total") or 0
    called = coverage.get("recent_called") or 0
    if total == 0:
        # No isolates in the window -- the NCBI release predates the window (stale release).
        # Show fill provenance so the reader knows data exists, just not in this window.
        n_fill = (provenance or {}).get("n_called")
        through = (provenance or {}).get("as_of")
        if n_fill and through:
            return (f"echinocandin calls: {int(n_fill):,} isolates called &middot; "
                    f"release current through <code>{html.escape(str(through))}</code>")
        return ""
    through = (provenance or {}).get("as_of")
    through_html = (f" &middot; calls current through <code>{html.escape(str(through))}</code>"
                    if through else "")
    if called == 0:
        n_fill = (provenance or {}).get("n_called")
        fill_size = (f" &middot; last batch: {int(n_fill):,} isolates called"
                     if n_fill else "")
        return (f'<span class="stale">echinocandin calls: awaiting re-call'
                f'{fill_size}{through_html}</span>')
    pct = round(100 * called / total)
    age_html = _age_clause(provenance, window_days, now_date)
    return (f"echinocandin calls: {called} of {total} window isolates re-called "
            f"({pct}%){through_html}{age_html}")


def _is_awaiting(coverage):
    """True when a window exists but no isolate in it has been re-called yet -- a lapsed or
    not-yet-current fill, distinct from a genuinely never-started day-0 watch."""
    return bool(coverage) and (coverage.get("recent_total") or 0) > 0 \
        and (coverage.get("recent_called") or 0) == 0


def _prox_display(prox):
    """`(tier_label, dist_str)` for a structural-proximity dict: the tier with underscores
    spaced out, and the min heavy-atom distance as `"12.3&thinsp;Å"` (or `"N/A"` when absent).
    Shared by the alert card and the VUS row so the two renderers format proximity identically
    even though their surrounding markup differs."""
    tier_label = prox["proximity_tier"].replace("_", " ")
    dist = prox.get("min_dist_to_drug_A")
    dist_str = f"{dist:.1f}&thinsp;Å" if dist else "N/A"
    return tier_label, dist_str


def _pct(x, *, decimals=0, dash="—"):
    """A proportion (0..1) as a percent string, or `dash` when absent. `decimals=0` gives a
    whole percent (the certification figures); `decimals=1` the one-decimal prevalence used by
    the trend and forecast sections. The single proportion formatter for the page."""
    if x is None:
        return dash
    if decimals == 0:
        return f"{round(100 * x)}%"
    return f"{100 * x:.{decimals}f}%"


def _plural(n, singular, plural=None):
    """`"1 carrier"` / `"3 carriers"` -- a count with its noun, correctly pluralised, so the
    page never falls back to the machine-looking `"carrier(s)"`. `plural` defaults to
    `singular + "s"`."""
    word = singular if n == 1 else (plural or f"{singular}s")
    return f"{n} {word}"


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


# What this service deliberately does NOT claim yet, and why -- the single source of truth for
# both the page section (_scope_limits_section) and the data.json `scope_limits` key
# (weather_feed.build_data), so the two can never disagree. These are product-level policy, not
# data-derived: each `reason` is the honest ceiling named by a frozen pre-registration / contract
# (cited in `source`), NOT a new claim. Surfacing the deferrals turns a built-but-dark capability
# into an honesty signal rather than a silent gap.
SCOPE_LIMITS = [
    {
        "claim": "a per-isolate calibrated resistance probability",
        "status": "deferred",
        "reason": "This report tells you whether a known echinocandin-resistance FKS1 marker is "
                  "present -- a yes/no detection, not a probability. We do not attach a number "
                  "like \"80% likely resistant\", because saying so honestly means first "
                  "measuring that such a number holds up: across isolates we would stamp 80%, "
                  "about 80% must actually be resistant by a lab test. No public collection "
                  "pairs the genomes with lab susceptibility results at the size needed to "
                  "measure that, so the probability stays off the report until it can be earned.",
        "unlock": "a paired genotype + lab-susceptibility panel that clears the pre-registered "
                  "reliability bar (work/PREREGISTRATION_calibration.md, "
                  "docs/DIAGNOSTIC_VERDICT_CONTRACT.md)",
        "source": "work/PREREGISTRATION_calibration.md",
    },
    {
        "claim": "an azole (fluconazole) resistance read-out",
        "status": "out of scope",
        "reason": "This is the echinocandin (FKS1) service. We do not report fluconazole "
                  "resistance here: azole resistance is already near-ubiquitous in C. auris, so "
                  "a presence/absence signal would carry little early-warning value.",
        "unlock": None,
        "source": "TODOS.md (FKS1 surveillance redirection)",
    },
]


def _scope_limits_section(scope_limits=SCOPE_LIMITS):
    """The collapsible 'what this service does not claim yet' block, or "" when the list is empty.

    Renders the SCOPE_LIMITS policy (shared with weather_feed so page and feed agree). Static --
    it introduces no number or probability; it states what is deliberately withheld and why,
    which is the honesty signal the certified-caller / Wilson-CI / forecast-skill discipline
    already sets elsewhere on the page.
    """
    if not scope_limits:
        return ""
    items = []
    for lim in scope_limits:
        claim = html.escape(str(lim.get("claim") or ""))
        status = html.escape(str(lim.get("status") or "deferred"))
        reason = html.escape(str(lim.get("reason") or ""))
        unlock = lim.get("unlock")
        unlock_html = (f' <span class="note">Earned by: {html.escape(str(unlock))}.</span>'
                       if unlock else "")
        items.append(
            f"      <dt>{claim} &mdash; <em>{status}</em></dt>\n"
            f"      <dd>{reason}{unlock_html}</dd>")
    rows = "\n".join(items)
    return f"""  <details class="about" id="scope-limits">
    <summary>What this report does not claim yet &middot; and why</summary>
    <p>Credibility rests as much on what we withhold as on what we show. These are claims this
      service could make look plausible but has deliberately held back until each is earned; the
      reasons are fixed in advance in the pre-registrations and contract this project freezes,
      not decided after the fact.</p>
    <dl class="grid scope-grid">
{rows}
    </dl>
  </details>
"""


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
    <p><b>How to read this page.</b> It reports whether a <em>known</em>
      echinocandin-resistance FKS1 marker is present in public <em>C. auris</em> genomes.
      That is detection, not a lab susceptibility test and not a treatment recommendation
      &mdash; it is an early-warning triage aid for research use, a signal to look closer.
      Every number comes straight from the certified caller and is reproducible from the
      committed call tables (<code>openafr.alert</code> renders the verdict, it never
      invents one). Two terms recur: an isolate is <b>resolved</b> when the caller could
      read its FKS1 sequence cleanly (prevalence is always a share of resolved isolates,
      never of raw deposits), and a <b>95% CI</b> is the range the true rate most likely
      sits in given the sample size. The echinocandin arm makes no structural fit verdict
      because echinocandins do not coordinate a metal, so the CYP51/heme-iron geometry
      behind the azole-arm verdict does not transfer to FKS1.</p>
    <p><b>Source.</b> {html.escape(_TAXON)} whole-genome deposits from NCBI Pathogen
      Detection, over a trailing <code>{html.escape(str(window))}</code>-day window.{fill_line}
      Public assemblies carry no echinocandin call (NCBI runs no AMR pipeline), so an
      out-of-band FKS1 re-call batch (cadence set by the ~97-day median deposit lag)
      produces the calls this page overlays.</p>
    <p><b>Certified caller.</b> The FKS1 re-caller is graded on a frozen external benchmark
      (<code>{bench}</code>, PMC12323592) on every change, through the same scorer this
      verdict uses:</p>
    <dl class="grid cert-grid">
      <dt>Verdict</dt><dd>{verdict} ({n_res}/{n} isolates resolved)</dd>
      <dt>Sensitivity</dt><dd>{sens}</dd>
      <dt>Specificity</dt><dd>{spec}</dd>
      <dt>Exact-token identity</dt><dd>{ident}</dd>
    </dl>
    <p class="note">The echinocandin arm reports presence of a known resistance variant,
      not a structural fit verdict; the priority tiers come from the same
      <code>openafr.alert</code> verdict.</p>
  </details>
"""


def _prevalence_sparkline(points, *, projection=None, aria_label=""):
    """A tiny, deterministic, no-JS inline-SVG sparkline of a prevalence series.

    `points` -- oldest-first list of (label, value) where value is a proportion 0..1 or None
    (a gap: NCBI deposit-month with no resolved isolate). Returns "" unless at least two points
    carry a real value -- one point is not a trajectory, matching the sections' own guards.

    `projection` -- when the forecast backtest earned it: {point, lo, hi} proportions for the
    next step, drawn as a hollow marker past the series with a faint 95%-PI whisker. Omitted
    entirely otherwise, so the chart never shows a forecast the backtest did not demonstrate.

    All colour comes from CSS classes that read the page's :root tokens (no hex here), so the
    chart inherits light / dark / print exactly like the rest of the page. The y-axis scales to
    the data range (the signal lives in a ~0-7% band, not 0-100%) with honest min/max end-ticks.
    """
    vals = [v for _, v in points if v is not None]
    if len(vals) < 2:
        return ""

    # Reserve a slot to the right for the projected point so it reads as "next", not in-series.
    proj = projection if (projection and projection.get("point") is not None) else None
    slots = len(points) + (1 if proj else 0)
    if slots < 2:
        return ""

    # Scale to the data range, widened to include the projection point + its PI so nothing clips.
    lo_v = min(vals)
    hi_v = max(vals)
    if proj:
        for key in ("point", "lo", "hi"):
            pv = proj.get(key)
            if pv is not None:
                lo_v, hi_v = min(lo_v, pv), max(hi_v, pv)
    span = (hi_v - lo_v) or 1.0  # flat series -> avoid divide-by-zero; line sits mid-plot

    # viewBox geometry: a short plot band with room for the two end-ticks on the right.
    X0, X1, Y0, Y1 = 6.0, 286.0, 10.0, 50.0

    def _x(i):
        return X0 + (i * (X1 - X0) / (slots - 1))

    def _y(v):
        return round(Y1 - ((v - lo_v) / span) * (Y1 - Y0), 1)

    # One polyline per contiguous run of real values, so a gap breaks the line (never draws to 0).
    runs, cur = [], []
    for i, (_, v) in enumerate(points):
        if v is None:
            if len(cur) > 1:
                runs.append(cur)
            cur = []
        else:
            cur.append((round(_x(i), 1), _y(v)))
    if len(cur) > 1:
        runs.append(cur)
    polylines = "".join(
        f'<polyline class="spark-line" points="{" ".join(f"{x},{y}" for x, y in run)}"/>'
        for run in runs)
    dots = "".join(
        f'<circle class="spark-dot" cx="{round(_x(i), 1)}" cy="{_y(v)}" r="2"/>'
        for i, (_, v) in enumerate(points) if v is not None)

    proj_svg = ""
    if proj:
        px = round(_x(slots - 1), 1)
        py = _y(proj["point"])
        # Dashed connector from the last real point to the projection.
        last = next(((round(_x(i), 1), _y(v)) for i in range(len(points) - 1, -1, -1)
                     for v in (points[i][1],) if v is not None), None)
        connector = (f'<line class="spark-proj-link" x1="{last[0]}" y1="{last[1]}" '
                     f'x2="{px}" y2="{py}"/>' if last else "")
        whisker = ""
        if proj.get("lo") is not None and proj.get("hi") is not None:
            whisker = (f'<line class="spark-pi" x1="{px}" y1="{_y(proj["lo"])}" '
                       f'x2="{px}" y2="{_y(proj["hi"])}"/>')
        proj_svg = (f"{connector}{whisker}"
                    f'<circle class="spark-proj" cx="{px}" cy="{py}" r="2.6"/>')

    # Honest end-ticks: the max value at the top-right, the min at the bottom-right.
    ticks = (f'<text class="spark-tick" x="292" y="{Y0 + 3}">{hi_v * 100:.1f}%</text>'
             f'<text class="spark-tick" x="292" y="{Y1}">{lo_v * 100:.1f}%</text>')
    label = html.escape(str(aria_label))
    return (f'<svg class="spark" viewBox="0 0 320 60" preserveAspectRatio="xMidYMid meet" '
            f'role="img" aria-label="{label}"><title>{label}</title>'
            f'<line class="spark-base" x1="{X0}" y1="{Y1}" x2="{X1}" y2="{Y1}"/>'
            f"{polylines}{dots}{proj_svg}{ticks}</svg>")


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
        return _pct(x, decimals=1, dash="&mdash;")

    rows = "\n".join(
        f'      <dt><code>{html.escape(str(e.get("as_of") or "—"))}</code></dt>'
        f'<dd>{_prev(e.get("prevalence"))} &mdash; '
        f'{e.get("n_resistant", 0)} of {_plural(e.get("n_called", 0), "isolate")}</dd>'
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
    vus_html = (f" A further {_plural(n_vus, 'resolved isolate')} carry a non-panel FKS1 "
                "variant of uncertain significance &mdash; watched, not counted as resistant."
                ) if n_vus else ""

    spark = _prevalence_sparkline(
        [(e.get("as_of"), e.get("prevalence")) for e in release_trend],
        aria_label="Known-resistance prevalence across committed re-call batches")
    spark_html = f"    {spark}\n" if spark else ""

    return f"""  <section class="trend">
    <h2>Resistance trend</h2>
    <p class="gloss">&ldquo;Prevalence&rdquo; is the share of successfully called isolates
      carrying a known resistance marker &mdash; not a count of cases.</p>
{spark_html}    <dl class="grid">
{rows}
    </dl>
{delta_html}    <p class="note">Carriers of a <strong>known</strong> echinocandin-resistance FKS1
      marker as a share of resolved isolates, one point per re-call batch.{vus_html}</p>
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
        return _pct(x, decimals=1, dash="&mdash;")

    rows = "\n".join(
        f'      <dt><code>{html.escape(str(e.get("month") or "—"))}</code></dt>'
        f'<dd>{_prev(e.get("prevalence"))} &mdash; '
        f'{e.get("n_resistant", 0)} of {_plural(e.get("n_resolved", 0), "isolate")}</dd>'
        for e in series)

    proj = forecast.get("projection") or {}
    # The chart only carries the projection marker when the backtest earned it -- same gate as
    # the text line below, so the sparkline never shows a forecast the data did not demonstrate.
    spark_proj = None
    if forecast.get("skill") and proj.get("point") is not None:
        spark_proj = {"point": proj.get("point"), "lo": proj.get("lo"), "hi": proj.get("hi")}
    spark = _prevalence_sparkline(
        [(e.get("month"), e.get("prevalence")) for e in series],
        projection=spark_proj,
        aria_label="Panel-gated resistance prevalence by NCBI deposit month")
    spark_html = f"    {spark}\n" if spark else ""
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
        note = ("The projection (one month ahead) is shown only because a pre-registered "
                "walk-forward test showed the trend model out-predicts naive persistence.")
    else:
        pooled = proj.get("naive_pooled")
        pooled_str = _prev(pooled) if pooled is not None else "the series mean"
        fc_html = (f'    <p class="trend-delta trend-flat">&#9644; No forecast skill '
                   f"demonstrated over a flat baseline at this series depth; showing the "
                   f"trajectory only (prevalence is roughly stationary near "
                   f"{pooled_str}).</p>\n")
        note = ("A pre-registered walk-forward test did not show the trend model beating a "
                "flat baseline on held-out months, so no projection is claimed.")

    return f"""  <section class="trend forecast">
    <h2>How prevalence is trending</h2>
{spark_html}    <dl class="grid">
{rows}
    </dl>
{fc_html}    <p class="note">Known-resistance carriers as a share of resolved isolates, by the
      month each genome was deposited to NCBI. {html.escape(note)}</p>
  </section>
"""


def _geo_section(geo_breakdown):
    """Where the current release's resistant carriers are: a location-ranked list, or "" when
    there is nothing located. Labelled 'current release' to distinguish it from the per-release
    historical series below; the breakdown itself is read reproducibly from the committed call
    table's geography (falling back to the live pull only when a table predates the backfill).

    `geo_breakdown` -- {location: carrier_count} from render_weather, keyed by the finest
    geography NCBI gave us (region where available, else country); passed in, not recomputed.
    """
    if not geo_breakdown:
        return ""
    items = sorted(geo_breakdown.items(), key=lambda kv: (-kv[1], kv[0]))
    rows = "\n".join(
        f"      <dt>{html.escape(str(location))}</dt><dd>{n}</dd>"
        for location, n in items)
    return f"""  <section class="geo">
    <h2>Where carriers are &middot; latest batch</h2>
    <dl class="grid">
{rows}
    </dl>
    <p class="note">Resistant FKS1 carriers by location in the latest batch, at the finest
      geography reported (region where available, else country). Isolates with no recorded
      location are omitted.</p>
  </section>
"""


def _geo_history_section(geo_history):
    """Where resistant carriers have been across the committed call series: one location-ranked
    block per release, or "" when there is no historical geography to show.

    Distinct from `_geo_section`, which shows only the current live NCBI release. This reads the
    committed call tables' `country` column (via render_weather's `_geo_history`), so it is
    reproducible from committed data -- a release whose table predates that column is simply
    absent, so the section only shows releases it can honestly place.

    `geo_history` -- [{release_tag, as_of, geography: {country: carrier_count}}], oldest first;
    passed in, not recomputed. None/empty renders nothing (honest degradation, same as the other
    sections).
    """
    if not geo_history:
        return ""
    blocks = []
    for entry in geo_history:
        geo = entry.get("geography") or {}
        if not geo:
            continue
        items = sorted(geo.items(), key=lambda kv: (-kv[1], kv[0]))
        rows = "\n".join(
            f"      <dt>{html.escape(str(loc))}</dt><dd>{n}</dd>"
            for loc, n in items)
        tag = html.escape(str(entry.get("release_tag") or ""))
        as_of = html.escape(str(entry.get("as_of") or ""))
        label = f"{tag} &middot; as of {as_of}" if as_of else tag
        blocks.append(f"""    <div class="geo-release">
      <p class="geo-release-label">{label}</p>
      <dl class="grid">
{rows}
      </dl>
    </div>""")
    if not blocks:
        return ""
    return f"""  <section class="geo geo-history">
    <h2>Where carriers are &middot; over time</h2>
{chr(10).join(blocks)}
    <p class="note">Known-resistance carriers by country for each re-call batch, so the
      geographic picture can be read over time, not just for the latest batch above
      (country-level, as NCBI reported it; isolates with no recorded country are omitted).</p>
  </section>
"""


def _geo_prevalence_section(geo_prevalence):
    """The 'where resistance is proportionally elevated' section: panel-gated echinocandin
    PREVALENCE (resistant / resolved, with a 95% Wilson CI) per place -- the proportional
    counterpart to `_geo_section`, which shows only raw carrier counts. Raw counts conflate
    sequencing volume with resistance; this section answers the actual surveillance question
    (is resistance proportionally high anywhere?), flagging a place as 'elevated' only when its
    CI lower bound clears the national headline rate. "" when there is nothing to show (honest
    degradation, same as every other section).

    `geo_prevalence` -- the dict from render_weather's `_geo_prevalence`: `places` (each with
    `place`, `n_resolved`, `n_panel_positive`, `event_frequency`, `ci95`, and `markers` -- the
    per-place marker tally, each `{token, carriers, low_ppv}`, descending by carriers), already
    ranked by prevalence; `pooled_small` (count + summed resolved/positive of places held below
    the minimum resolved-isolate threshold, named in the footnote so nothing is hidden); and
    `headline_prevalence` (the national rate from the same resolved rows, the single source of
    truth the 'elevated' flag compares against). Passed in, never recomputed.
    """
    if not geo_prevalence or not geo_prevalence.get("places"):
        return ""
    headline = geo_prevalence.get("headline_prevalence")

    def _markers_clause(markers):
        # Name which marker(s) drive the place: the top one by carriers, plus a second only
        # when it is within half the top's count (a genuinely mixed place). A token-carrier
        # tally, never a fraction of the place's positives, so it can't read as >100%.
        if not markers:
            return ""
        top = markers[0]["carriers"]
        shown = [m for m in markers if m["carriers"] * 2 >= top][:2]
        parts = []
        for m in shown:
            tok = html.escape(str(m["token"]))
            low = ' <span class="low-ppv">weaker marker</span>' if m.get("low_ppv") else ""
            parts.append(f"<code>{tok}</code>&times;{m['carriers']}{low}")
        return f' &middot; driven by {", ".join(parts)}'

    def _row(p):
        place = html.escape(str(p["place"]))
        pct = f"{p['event_frequency'] * 100:.1f}%"
        ci = p.get("ci95")
        ci_part = (f" (95% CI {ci[0] * 100:.1f}&ndash;{ci[1] * 100:.1f}%)"
                   if ci else "")
        n = f"{p['n_panel_positive']} of {_plural(p['n_resolved'], 'isolate')}"
        # 'Elevated' is the strong claim -- the CI lower bound must clear the national rate, so
        # a place is only flagged when its resistance is statistically above the headline, not
        # merely numerically higher on a small sample.
        elevated = (ci and headline is not None and ci[0] > headline)
        tag = ' <span class="elevated">above national</span>' if elevated else ""
        markers = _markers_clause(p.get("markers") or [])
        return (f"      <dt>{place}{tag}</dt>"
                f"<dd>{pct}{ci_part} &mdash; {n}{markers}</dd>")

    rows = "\n".join(_row(p) for p in geo_prevalence["places"])
    has_low_ppv = any(m.get("low_ppv") for p in geo_prevalence["places"]
                      for m in (p.get("markers") or []))
    pooled = geo_prevalence.get("pooled_small") or {}
    n_small = pooled.get("n_places", 0)
    pooled_note = (
        f" {_plural(n_small, 'smaller place')} with too few resolved isolates to rate "
        f"{'is' if n_small == 1 else 'are'} pooled out "
        f"({_plural(pooled.get('n_panel_positive', 0), 'carrier')} over "
        f"{_plural(pooled.get('n_resolved', 0), 'resolved isolate')}) rather than shown as "
        f"unstable rates."
        if n_small else "")
    headline_note = (
        f" The national rate is {headline * 100:.1f}%; a place is flagged "
        f"<strong>above national</strong> only when its 95% CI lower bound clears it."
        if headline is not None else "")
    markers_note = (
        " &ldquo;Driven by&rdquo; names which known-resistance marker(s) the place's carriers "
        "hold &mdash; distinct dominant markers in two places point to independent emergence, "
        "not one spreading lineage."
        + (" A marker tagged a <strong>weaker marker</strong> (D642Y, M690I) still counts, but "
           "its carriers split resistant/susceptible rather than reliably tracking resistance."
           if has_low_ppv else ""))
    return f"""  <section class="geo-prevalence">
    <h2>Where resistance is highest (by share, not just counts)</h2>
    <p class="gloss">Echinocandin-resistance prevalence by place &mdash; resistant isolates as a
      share of those the caller could resolve, not raw carrier counts, so a heavily-sequenced
      region doesn't look hot from volume alone.</p>
    <dl class="grid">
{rows}
    </dl>
    <p class="note">Resistant isolates as a share of resolved isolates for each place, with a
      95% confidence interval.{headline_note}{markers_note}{pooled_note}</p>
  </section>
"""


def _panel_breakdown_section(panel_census):
    """The 'which known-resistance variant is circulating' section: the per-marker breakdown
    behind the single headline prevalence number -- each known-resistance panel substitution
    the current table carries, with how many isolates carry it -- or "" when there are none
    (honest degradation, same as trend/geo/VUS).

    `panel_census` -- the per-variant list from render_weather: each entry carries `token`,
    `carriers`, and `low_ppv` (True for the known-but-low-PPV positions D642Y/M690I, which are
    counted in the headline but split R/S rather than tracking resistance). Passed in, never
    recomputed here. Carrier counts sum to the headline `n_resistant`; this is detection detail,
    strictly the marker and its count, never a per-isolate phenotype claim.
    """
    if not panel_census:
        return ""
    def _panel_row(v):
        token = html.escape(str(v["token"]))
        carriers = v.get("carriers", 0)
        low_ppv = v.get("low_ppv")
        tag = (' <span class="low-ppv">weaker marker (low PPV)</span>' if low_ppv else "")
        return (
            f"      <dt><code>{token}</code>{tag}</dt>"
            f"<dd>{_plural(carriers, 'carrier')}</dd>"
        )

    rows = "\n".join(_panel_row(v) for v in panel_census)
    has_low_ppv = any(v.get("low_ppv") for v in panel_census)
    low_ppv_note = (
        " Markers tagged a <strong>weaker marker</strong> (D642Y, M690I) still count toward the "
        "headline, but their carriers split resistant/susceptible rather than reliably tracking "
        "resistance &mdash; read them with more caution than the core markers."
        if has_low_ppv else "")
    return f"""  <section class="panel-breakdown">
    <h2>Which resistance variant is circulating</h2>
    <p class="gloss">The known echinocandin-resistance FKS1 markers behind the headline
      prevalence &mdash; which specific substitution each resistant isolate carries, and how
      many carry it.</p>
    <dl class="grid">
{rows}
    </dl>
    <p class="note">How many isolates carry each known-resistance marker. An isolate carrying
      more than one marker is counted under each.{low_ppv_note}</p>
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
            tier_label, dist_str = _prox_display(prox)
            tier_label = html.escape(tier_label)
            prox_html = (
                f' <span class="prox-tier">{dist_str} from the drug ({tier_label})</span>'
            )
        sv = ctx.get("structural_verdict")
        sv_html = ""
        if sv:
            label = html.escape(sv.get("label", ""))
            sv_html = f' <span class="struct-label">structural: {label}</span>'
        return (
            f'      <dt><code>{token}</code> '
            f'<span class="vus-tier">{tier_str}</span>{prox_html}{sv_html}</dt>'
            f"<dd>{_plural(carriers, 'carrier')} &mdash; {note}</dd>"
        )

    rows = "\n".join(_vus_row(v) for v in vus_census)
    return f"""  <section class="vus">
    <h2>Variants we're watching but not yet counting (uncertain significance)</h2>
    <p class="gloss">Hot-spot DNA changes we have detected but that are <em>not yet</em> known
      to cause resistance &mdash; flagged and watched, deliberately not counted as resistant.</p>
    <dl class="vus-list">
{rows}
    </dl>
    <p class="note">Non-panel FKS1 hot-spot substitutions seen in the calls:
      <strong>watched, not counted as resistant</strong>. The position shows where the residue
      sits relative to the known-resistance hot spots; where available, the distance to the
      bound drug (enfumafungin &mdash; the same FKS1 site echinocandins act on, in the 9WZU
      structure) and a spatial co-location label are shown. These locate a variant; they are
      not a drug-fit or resistance prediction.</p>
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
                certification=None, release_trend=None, geo_breakdown=None, geo_history=None,
                geo_prevalence=None, panel_census=None, vus_census=None, forecast=None,
                page_url=None, feeds=False):
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
    `geo_history` -- [{release_tag, as_of, geography: {country: carrier_count}}], oldest first,
                        from the committed call series; renders the per-release historical
                        geography section. None/empty omits it.
    `geo_prevalence` -- {places, pooled_small, headline_prevalence} from the committed table:
                        panel-gated prevalence (resistant / resolved, with a Wilson CI) per
                        place; renders the "where resistance is proportionally elevated" section.
                        None/empty omits it.
    `panel_census` -- the per-variant breakdown of KNOWN-resistance panel markers, each with
                        `token`, `carriers`, `low_ppv`; renders the "which resistance variant
                        is circulating" section. Carrier counts sum to the headline
                        `n_resistant`. None/empty omits it.
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
        latest_fill = (release_trend or [{}])[-1] if release_trend else {}
        n_fill = latest_fill.get("n_called")
        n_res = latest_fill.get("n_resistant")
        prev = latest_fill.get("prevalence")
        if n_fill is not None and prev is not None:
            prev_pct = f"{prev * 100:.1f}%"
            n_res_part = (f" ({_plural(n_res, 'carrier')} of a known marker)"
                          if n_res is not None else "")
            fill_note = (
                f'<p class="meta-detail">Most recent batch: {int(n_fill):,} isolates called'
                f' &middot; {prev_pct} resistance prevalence{n_res_part}.</p>'
            )
        else:
            fill_note = ""
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
      {fill_note}
      {note_html}
    </article>"""

    as_of = html.escape(str(r.get("as_of", "—")))
    window = html.escape(str(r.get("window_days", "—")))
    coverage_clause = _coverage_clause(coverage, calls_provenance,
                                       now_date=now_date, window_days=r.get("window_days"))
    # The freshness/coverage clause is the single most important fact on the meta line, so
    # it leads on its own readable line; the provenance detail (checked/as-of/window/source)
    # follows, muted. Same facts and the same honest degradation as before, better hierarchy.
    fresh_html = (f'  <p class="meta-fresh">{coverage_clause}</p>\n'
                  if coverage_clause else "")

    cert_badge = _cert_badge(certification)
    trend_html = _trend_section(release_trend)
    forecast_html = _forecast_section(forecast)
    panel_html = _panel_breakdown_section(panel_census)
    geo_html = _geo_section(geo_breakdown)
    geo_prevalence_html = _geo_prevalence_section(geo_prevalence)
    geo_history_html = _geo_history_section(geo_history)
    vus_html = _vus_section(vus_census)
    scope_limits_html = _scope_limits_section()
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
    color-scheme: light dark;
    --bg:#0d1017; --surface:#151a23; --surface-accent:#132019; --text:#e6e9ef;
    --muted:#8a94a6; --muted-strong:#9aa4b2; --border:#232a36; --code:#c3cad6;
    --accent:#4dd6a0; --accent-border:#2a5a47; --link:#7aa2ff; --stale:#e0b050;
    --footer:#5c6675;
    --tier-act-now:#ff5c5c; --tier-watch:#ffcc4d; --tier-context:#8a94a6;
    --accent-page:{page_accent};
  }}
  /* Light palette for bright clinical screens / daytime reading. Every rule reads the
     tokens above, so re-defining them here is the whole theme -- no other rule changes.
     Values clear WCAG AA contrast on --text/--muted over --bg and keep the tier hues
     (act-now red / watch amber / context grey) legible on a light surface. */
  @media (prefers-color-scheme: light) {{
    :root {{
      --bg:#f7f8fa; --surface:#ffffff; --surface-accent:#eef7f1; --text:#1a1f28;
      --muted:#5a6475; --muted-strong:#3c4452; --border:#d9dee6; --code:#2d3748;
      --accent:#0f8a5f; --accent-border:#9fd6bd; --link:#2a5bd7; --stale:#9a6a00;
      --footer:#8a94a6;
      --tier-act-now:#d12b2b; --tier-watch:#9a6a00; --tier-context:#5a6475;
    }}
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
  .meta .dark {{ color:var(--tier-act-now); font-weight:600; }}
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
  /* Scope-limits read as stacked heading + paragraph, not a two-column key:value grid. */
  .about .scope-grid {{ display:block; max-width:640px; margin:.6rem 0; }}
  .about .scope-grid dt {{ color:var(--muted-strong); font-weight:600; margin:.9rem 0 .2rem; }}
  .about .scope-grid dt em {{ color:var(--stale); font-style:normal; font-weight:600; }}
  .about .scope-grid dd {{ margin:0 0 .3rem; }}
  .about .scope-grid dd .note {{ display:block; margin-top:.15rem; }}
  h2 {{ font-size:.95rem; letter-spacing:.02em; text-transform:uppercase;
        color:var(--muted-strong); font-weight:600; margin:2.5rem 0 .6rem; }}
  .trend, .geo, .vus {{ margin:0 0 1rem; }}
  .geo-release {{ margin:.3rem 0 1rem; }}
  .geo-release-label {{ color:var(--muted-strong); font-size:.82rem; font-weight:600;
                        margin:.6rem 0 .3rem; }}
  .vus-tier {{ color:var(--muted); font-size:.78rem; }}
  .struct-label {{ color:var(--accent); font-size:.78rem; }}
  .low-ppv {{ color:var(--muted); font-size:.72rem; font-weight:600; text-transform:uppercase;
              letter-spacing:.04em; border:1px solid var(--muted); border-radius:3px;
              padding:0 .3em; margin-left:.35em; vertical-align:middle; }}
  .elevated {{ color:var(--tier-watch); font-size:.72rem; font-weight:600;
               text-transform:uppercase; letter-spacing:.04em;
               border:1px solid var(--tier-watch); border-radius:3px;
               padding:0 .3em; margin-left:.35em; vertical-align:middle; }}
  .trend-delta {{ font-size:.9rem; margin:.5rem 0 0; }}
  .trend-up {{ color:var(--tier-act-now); }}
  .trend-down {{ color:var(--accent); }}
  .trend-flat {{ color:var(--muted); }}
  /* Inline prevalence sparkline. All colour via tokens so it tracks light/dark/print with the
     rest of the page; no JS, no external asset. */
  .spark {{ display:block; width:100%; max-width:420px; height:auto; margin:.2rem 0 .7rem;
            overflow:visible; }}
  .spark-base {{ stroke:var(--border); stroke-width:1; }}
  .spark-line {{ fill:none; stroke:var(--accent); stroke-width:1.6;
                 stroke-linejoin:round; stroke-linecap:round; }}
  .spark-dot {{ fill:var(--accent); }}
  .spark-proj-link {{ stroke:var(--muted); stroke-width:1.4; stroke-dasharray:3 2; }}
  .spark-pi {{ stroke:var(--muted); stroke-width:1.4; opacity:.5; }}
  .spark-proj {{ fill:var(--bg); stroke:var(--muted-strong); stroke-width:1.6; }}
  .spark-tick {{ fill:var(--muted); font-size:7px; }}
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
           font-size:.9rem; margin:0; }}
  .grid dt {{ color:var(--muted); }}
  .grid dd {{ margin:0; }}
  /* VUS rows carry a long token+tags header and a full-sentence note, so they stack
     (header line, note beneath) instead of the two-column label/value grid, which would
     squeeze the note into a cramped far-right column. */
  .vus-list {{ font-size:.9rem; margin:0; }}
  .vus-list dt {{ margin-top:.7rem; }}
  .vus-list dt:first-of-type {{ margin-top:0; }}
  .vus-list dd {{ margin:.1rem 0 0; }}
  .pocket {{ display:block; width:100%; max-width:520px; margin:1rem 0 0;
             border-radius:8px; border:1px solid var(--border); }}
  .note {{ color:var(--muted); font-size:.85rem; font-style:italic; }}
  footer {{ margin-top:3rem; color:var(--footer); font-size:.78rem; }}
  a {{ color:var(--link); }}
  .lede {{ color:var(--muted-strong); font-size:.95rem; margin:0 0 1.25rem;
           max-width:60ch; }}
  .gloss {{ color:var(--muted); font-size:.82rem; font-style:italic;
            margin:.15rem 0 0; }}
  /* Print: ink-light, readable on paper. Force the light tokens regardless of the
     reader's OS theme, drop surface fills, and open the collapsed methods details so a
     printed surveillance report is self-contained. */
  @media print {{
    :root {{
      --bg:#ffffff; --surface:#ffffff; --surface-accent:#ffffff; --text:#000000;
      --muted:#333333; --muted-strong:#000000; --border:#999999; --code:#000000;
      --accent:#006a44; --accent-border:#999999; --link:#000000; --stale:#664400;
      --footer:#333333;
      --tier-act-now:#a00000; --tier-watch:#664400; --tier-context:#333333;
    }}
    body {{ background:#fff; color:#000; }}
    .card {{ background:none; }}
    .cert a {{ background:none; }}
    .about {{ display:block; }}
    .about summary {{ display:none; }}
    a {{ color:#000; text-decoration:underline; }}
  }}
</style>
</head>
<body>
<main class="wrap">
  <h1>C. auris Resistance Weather Report</h1>
  <p class="lede">A free, weekly-updated public read on emerging echinocandin
    (first-line antifungal) resistance in <em>Candida auris</em>, a WHO critical-priority
    hospital pathogen. It is built for genomic-surveillance and public-health teams as an
    early-warning triage aid &mdash; a signal to look closer, not a clinical test result or a
    treatment recommendation.</p>
  <p class="status">{html.escape(status)}</p>
  <div class="meta">
{fresh_html}    <p class="meta-detail">Last checked <code>{last_checked}</code> &middot;
       as of <code>{as_of}</code> &middot; window <code>{window}</code>d &middot;
       source: NCBI Pathogen Detection via OpenAFR early-warning.</p>
  </div>
{cert_badge}{body}
{trend_html}{forecast_html}{panel_html}{geo_html}{geo_prevalence_html}{geo_history_html}{vus_html}{scope_limits_html}{methods_html}  <footer>
    OpenAFR early-warning. Research use only — a triage aid, not a clinical claim.
    Echinocandin calls come from a monthly FKS1 re-call batch over the live NCBI feed
    (cadence set by the ~97-day median deposit lag); between batches the coverage line above
    states how current they are. See the methods section for how the numbers are produced.
{feeds_html}  </footer>
</main>
</body>
</html>
"""
