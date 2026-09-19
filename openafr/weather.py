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
    """One alert -> an HTML card. Reads only fields alert.compose_alerts guarantees."""
    tier_txt, accent = _TIER.get(a["priority"], (a["priority"].upper(), "#8a94a6"))
    e = a["emergence"]
    s = a["structural"]
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
        ("Structural fit",
         f"[{html.escape(str(s['evidence_class']))}, conf={html.escape(str(s['confidence']))}] "
         f"{html.escape(str(s['fluconazole_fit_verdict']))}"),
    ]
    rows_html = "\n".join(
        f'      <div class="k">{k}</div><div class="v">{v}</div>' for k, v in rows)
    img = pocket_asset(a["mutation"])
    headline = html.escape(str(a.get("headline") or ""))
    return f"""    <article class="card" style="--accent:{accent}">
      <header><span class="badge" style="background:{accent}">{tier_txt}</span>
        <code class="mut">{mut}</code></header>
      <p class="headline">{headline}</p>
      <div class="grid">
{rows_html}
      </div>
      <img class="pocket" src="{img}" alt="Mutant vs wild-type CYP51 pocket for {mut}"
           loading="lazy" onerror="this.style.display='none'">
    </article>"""


def render_page(alert_result, *, now=None, watching_since=None):
    """Render the whole page as a self-contained HTML string (inline CSS, no framework).

    `alert_result` -- the dict from alert.compose_alerts (or compose_fks1_alerts).
    `now`          -- a datetime/ISO string for "last checked"; defaults to UTC now.
    `watching_since` -- ISO date the watch started / last-delivered date, for the quiet
                        framing. None renders the pure day-0 line.
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

    if alerts:
        body = "\n".join(_alert_card(a) for a in alerts)
    else:
        note = r.get("note")
        since = html.escape(str(watching_since)) if watching_since else "the first run"
        note_html = (f'<p class="note">{html.escape(str(note))}</p>' if note else "")
        body = f"""    <article class="card calm">
      <p class="headline">No azole-resistance variant is over threshold right now.</p>
      <p>Watching {html.escape(_TAXON)} genome deposits since {since}. When a
         resistance-conferring ERG11 substitution emerges, it appears here with its
         structural fit verdict.</p>
      {note_html}
    </article>"""

    as_of = html.escape(str(r.get("as_of", "—")))
    window = html.escape(str(r.get("window_days", "—")))
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>C. auris Resistance Weather Report</title>
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
     as of <code>{as_of}</code> &middot; window <code>{window}</code>d &middot;
     source: NCBI Pathogen Detection via OpenAFR early-warning.</p>
{body}
  <footer>
    OpenAFR early-warning. Research use only — a triage aid, not a clinical claim.
    Priority tiers and every number come straight from
    <code>openafr.alert</code>; this page renders that verdict, it never invents one.
  </footer>
</div>
</body>
</html>
"""
