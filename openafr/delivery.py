"""Delivery + no-op decision for the C. auris early-warning track (issue #26).

Pipeline stage: the delivery layer -- makes the composed alert (#25) find a human,
on a cadence, minimally.

Gene-agnostic. The decision logic (new-vs-known, escalation, quiet no-op, committed
state) is identical for both tracks; only two rendered labels and which composer's
renderer to reuse depend on the gene, and those are read off the composed result's own
`gene` marker (ERG11/azole with a structural so-what, or FKS1/echinocandin which is
detection-only). Run the two tracks with separate state files so their memories do not
collide (scripts/deliver_alert.py --gene).

Why this exists
---------------
Everything upstream (#21-#25) turns an NCBI pull into a ranked, human-readable
alert. But an alert nobody sees is not an early warning. This module answers the
one question a scheduled runner needs: *given what we have already told people, is
there anything new enough this run to be worth interrupting a human for?* -- and if
so, renders the two artefacts the schedule delivers: a digest committed to the repo
and the body of an auto-filed GitHub issue.

The quiet no-op is the whole point
----------------------------------
A watch that cries every week is a watch nobody reads. The overwhelming majority of
runs will surface the same endemic mutations at the same tier they were last seen --
and those runs must produce *nothing*: no commit, no issue, no noise. We deliver only
on genuine news, defined narrowly and honestly:

  * a mutation we have never delivered an alert for before, or
  * a mutation whose so-what tier has *escalated* above the most urgent tier we have
    ever delivered for it (watch -> act-now, context -> watch/act-now).

A de-escalation (act-now -> watch) is deliberately NOT news: we already raised the
louder alarm; walking it back quietly is correct, not a fresh interruption. This
"most-urgent-ever" memory lives in a small committed state file so the decision is
reproducible and git-diffable, and so a fresh checkout picks up exactly where the
last delivery left off.

Why state is a committed file, not a database
---------------------------------------------
The schedule (a GitHub Action) is stateless between runs; the repo is the only
memory it carries forward. STATE.json is tiny (one line per mutation ever alerted),
so committing it costs nothing and its git history *is* the delivery log: every
escalation shows up as a one-line diff with the run that made it.

Stdlib only (json + datetime), matching the rest of the track.
"""
import datetime as _dt
import json

from openafr import alert

# Tier ordering, reused from the composer so "more urgent" means one thing in the
# whole track. Lower rank = more urgent (act-now=0, watch=1, context=2).
RANK = alert.PRIORITY


def _rank(priority):
    return RANK[priority]


def load_state(text):
    """Parse a STATE.json body into the delivery memory, tolerant of a fresh start.

    Returns {"delivered": {mutation: most_urgent_priority_ever}, "updated": iso}.
    An empty/None body (never delivered before) yields an empty, valid state rather
    than an error, so the first scheduled run needs no seed file.
    """
    if not text or not text.strip():
        return {"delivered": {}, "updated": ""}
    data = json.loads(text)
    delivered = data.get("delivered") or {}
    # Drop any tier token we no longer understand rather than trust it.
    delivered = {m: p for m, p in delivered.items() if p in RANK}
    return {"delivered": delivered, "updated": data.get("updated", "")}


def dumps_state(state):
    """Serialise state to a byte-stable JSON string (sorted keys, trailing newline)."""
    return json.dumps(state, indent=2, sort_keys=True) + "\n"


def new_alerts(alert_result, state):
    """The alerts in this run that are genuinely new news relative to `state`.

    New = mutation never delivered, OR its tier is strictly more urgent than the most
    urgent tier ever delivered for it. Returns the subset of `alert_result["alerts"]`
    (in the composer's existing rank order) that clears that bar. Empty list => no-op.
    """
    delivered = state.get("delivered", {})
    fresh = []
    for a in alert_result.get("alerts", []):
        mut, prio = a["mutation"], a["priority"]
        seen = delivered.get(mut)
        if seen is None or _rank(prio) < _rank(seen):
            fresh.append(a)
    return fresh


def is_noop(alert_result, state):
    """True when this run has no new news and should deliver absolutely nothing."""
    return not new_alerts(alert_result, state)


def next_state(alert_result, state, as_of=None):
    """State after a delivering run: remember the most urgent tier ever per mutation.

    Only current alerts update memory; mutations no longer present keep their recorded
    tier (history is never forgotten, so a mutation that vanishes and returns at the
    same tier stays quiet). `as_of` stamps the update for the git log; defaults to the
    alert's own as-of date so the state and the digest agree.
    """
    delivered = dict(state.get("delivered", {}))
    for a in alert_result.get("alerts", []):
        mut, prio = a["mutation"], a["priority"]
        seen = delivered.get(mut)
        if seen is None or _rank(prio) < _rank(seen):
            delivered[mut] = prio
    updated = as_of or alert_result.get("as_of") or _dt.date.today().isoformat()
    return {"delivered": delivered, "updated": updated}


# ------------------------------- rendering ------------------------------------
#
# Delivery is gene-agnostic except for two labels and which composer's renderer to
# reuse. The composed `result` already carries its own gene marker (compose_fks1_alerts
# sets gene="FKS1"; the ERG11 compose_alerts sets none), so we route off that rather
# than threading a parameter -- the digest can never disagree with the alert it wraps.
_GENE_META = {
    "ERG11": {"digest_title": "C. auris azole early-warning digest",
              "issue_gene": "ERG11", "render": alert.render_markdown},
    "FKS1": {"digest_title": "C. auris echinocandin (FKS1) early-warning digest",
             "issue_gene": "FKS1", "render": alert.render_fks1_markdown},
}


def _gene_meta(alert_result):
    """The per-gene labels + renderer for a composed result (default ERG11)."""
    return _GENE_META.get(alert_result.get("gene", "ERG11"), _GENE_META["ERG11"])


def _whats_new_block(fresh, state):
    """One-line-per-mutation summary of what made this run deliver.

    Distinguishes a brand-new mutation from an escalation of a known one, naming the
    prior tier for the escalation so the reader sees the movement, not just the level.
    """
    delivered = state.get("delivered", {})
    lines = []
    for a in fresh:
        mut, prio = a["mutation"], a["priority"]
        prior = delivered.get(mut)
        badge = alert._BADGE[prio]
        if prior is None:
            lines.append(f"- **`{mut}`** — new alert ({badge})")
        else:
            lines.append(f"- **`{mut}`** — escalated {alert._BADGE[prior]} → {badge}")
    return lines


def render_digest(alert_result, fresh, state):
    """The committed digest: a 'what's new this run' header over the full brief.

    The header names exactly what tripped delivery (so the git diff is self-explaining);
    the body is the composer's full Markdown so the committed file always shows the
    complete current picture, not just the delta.
    """
    meta = _gene_meta(alert_result)
    header = [f"# {meta['digest_title']}",
              "",
              f"_Delivered {alert_result.get('as_of', '')}_ — "
              f"{len(fresh)} new/escalated alert(s) this run.",
              "",
              "## What's new this run",
              ""]
    header += _whats_new_block(fresh, state)
    header += ["", "---", ""]
    return "\n".join(header) + "\n" + meta["render"](alert_result)


def issue_title(alert_result, fresh):
    """A scannable GitHub-issue title: count, tiers, date. No jargon."""
    n = len(fresh)
    tiers = sorted({a["priority"] for a in fresh}, key=_rank)
    badges = "/".join(alert._BADGE[t] for t in tiers)
    plural = "alert" if n == 1 else "alerts"
    gene = _gene_meta(alert_result)["issue_gene"]
    return (f"[early-warning] {n} new C. auris {gene} {plural} ({badges}) "
            f"— {alert_result.get('as_of', '')}")


def issue_body(alert_result, fresh, state):
    """The GitHub-issue body: the same digest a human would read in the repo."""
    return render_digest(alert_result, fresh, state)
