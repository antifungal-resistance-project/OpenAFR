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

# Tier -> (badge text, accent colour). Keys match alert.PRIORITY exactly.
_TIER = {
    "act-now": ("ACT-NOW", "#ff5c5c"),
    "watch": ("WATCH", "#ffcc4d"),
    "context": ("CONTEXT", "#8a94a6"),
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
        parts = [f"{counts[t]} {_TIER[t][0].lower()}" for t in _TIER if counts[t]]
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
    tier_txt, accent = _TIER.get(a["priority"], (a["priority"].upper(), "#8a94a6"))
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
    return f"""    <article class="card" style="--accent:{accent}">
      <header><span class="badge" style="background:{accent}">{tier_txt}</span>
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


def render_page(alert_result, *, now=None, watching_since=None, calls_provenance=None,
                certification=None, page_url=None):
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
    `page_url` -- canonical URL for the shareable OpenGraph/Twitter card; None omits og:url.
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
    accent = _TIER[alerts[0]["priority"]][1] if alerts else "#4dd6a0"

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
    coverage_html = f"\n     {coverage_clause} &middot;" if coverage_clause else ""

    cert_badge = _cert_badge(certification)
    methods_html = _methods_section(certification, calls_provenance, window)
    share_desc = html.escape(_share_description(status, coverage, calls_provenance))
    title = "C. auris Resistance Weather Report"
    og_url_tag = (f'\n<meta property="og:url" content="{html.escape(str(page_url))}">'
                  if page_url else "")
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<meta name="description" content="{share_desc}">
<meta property="og:title" content="{title}">
<meta property="og:description" content="{share_desc}">
<meta property="og:type" content="website">{og_url_tag}
<meta name="twitter:card" content="summary">
<meta name="twitter:title" content="{title}">
<meta name="twitter:description" content="{share_desc}">
<style>
  :root {{ color-scheme: dark; }}
  body {{ margin:0; background:#0d1017; color:#e6e9ef;
         font:16px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif; }}
  .wrap {{ max-width:760px; margin:0 auto; padding:2rem 1.25rem 4rem; }}
  h1 {{ font-size:1.15rem; letter-spacing:.02em; margin:0 0 .25rem;
        text-transform:uppercase; color:#9aa4b2; font-weight:600; }}
  .status {{ font-size:1.6rem; font-weight:700; margin:.2rem 0 1rem;
             color:{accent}; }}
  .meta {{ color:#8a94a6; font-size:.85rem; margin-bottom:2rem; }}
  .meta code {{ color:#c3cad6; }}
  .meta .stale {{ color:#e0b050; }}
  .cert {{ margin:-.5rem 0 2rem; font-size:.85rem; }}
  .cert a {{ display:inline-block; color:#4dd6a0; text-decoration:none;
             border:1px solid #2a5a47; background:#132019; border-radius:999px;
             padding:.2rem .7rem; }}
  .cert a:hover {{ border-color:#4dd6a0; }}
  .about {{ margin-top:2.5rem; border-top:1px solid #232a36; padding-top:1rem;
            color:#8a94a6; font-size:.85rem; }}
  .about summary {{ cursor:pointer; color:#9aa4b2; font-weight:600; }}
  .about code {{ color:#c3cad6; }}
  .about .cert-grid {{ max-width:420px; margin:.6rem 0; }}
  .card {{ background:#151a23; border:1px solid #232a36; border-left:4px solid var(--accent,#232a36);
           border-radius:10px; padding:1.1rem 1.25rem; margin:1rem 0; }}
  .card.calm {{ --accent:#4dd6a0; }}
  header {{ display:flex; align-items:center; gap:.6rem; margin-bottom:.5rem; }}
  .badge {{ font-size:.7rem; font-weight:700; letter-spacing:.05em;
            padding:.15rem .5rem; border-radius:999px; color:#0d1017; }}
  .mut {{ font-size:1.05rem; color:#e6e9ef; }}
  .headline {{ margin:.3rem 0 .8rem; }}
  .grid {{ display:grid; grid-template-columns:max-content 1fr; gap:.3rem .9rem;
           font-size:.9rem; }}
  .grid .k {{ color:#8a94a6; }}
  .pocket {{ display:block; width:100%; max-width:520px; margin:1rem 0 0;
             border-radius:8px; border:1px solid #232a36; }}
  .note {{ color:#8a94a6; font-size:.85rem; font-style:italic; }}
  footer {{ margin-top:3rem; color:#5c6675; font-size:.78rem; }}
  a {{ color:#7aa2ff; }}
</style>
</head>
<body>
<div class="wrap">
  <h1>C. auris Resistance Weather Report</h1>
  <p class="status">{html.escape(status)}</p>
  <p class="meta">Last checked <code>{last_checked}</code> &middot;
     as of <code>{as_of}</code> &middot; window <code>{window}</code>d &middot;{coverage_html}
     source: NCBI Pathogen Detection via OpenAFR early-warning.</p>
{cert_badge}{body}
{methods_html}  <footer>
    OpenAFR early-warning. Research use only — a triage aid, not a clinical claim.
    Echinocandin calls come from a monthly FKS1 re-call batch over the live NCBI feed
    (cadence set by the ~97-day median deposit lag); between batches the coverage line above
    states how current they are. Priority tiers and every number come straight from
    <code>openafr.alert</code>; this page renders that verdict, it never invents one.
  </footer>
</div>
</body>
</html>
"""
