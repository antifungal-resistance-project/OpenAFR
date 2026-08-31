"""Alert composition for the C. auris azole early-warning track (issue #25).

Pipeline stage: the delivery tip -- the one screen a non-bioinformatician acts on.

Why this exists
---------------
The three stages before this each answer one question and stop:

  * emergence (#22) -- is this ERG11 substitution NEW or RISING, and is the spike
    real or an archival backlog artifact?
  * mapping (#23)   -- WHERE does the residue sit in the modeled azole pocket, and
    can we build the mutant?
  * structural (#24)-- so what: does fluconazole still FIT the pocket, measured or
    estimated, and with how much confidence?

None of them is an *alert*. An alert is the assembled, human-readable statement a
clinician or public-health analyst can act on without opening a single PDB file:

  "F126L is newly appearing on 4 isolates across 2 countries, first collected 2024;
   it lines the drug-contact surface and predicts weaker azole fit (estimate)."

This module joins the emergence signal to the structural verdict, enriches it with
the epidemiological facts the detector dropped (which regions, first-seen date, how
many carriers in total), assigns a so-what priority, and renders the result two ways:
Markdown for a human and JSON for a machine.

The join, honestly
------------------
The emergence signal carries only `carrier_keys` (the recent-window carriers). To
name regions and a first-seen date we re-join those keys against the full record set,
and we compute first-seen over *all* carriers of the mutation (baseline + recent +
undated), not just the recent window -- so "first seen 2019" reflects the earliest
sample we actually hold, not the window this run happened to scan.

Two dates, kept separate (the #22 discipline carried forward)
-------------------------------------------------------------
NCBI visibility (`target_creation_date`) is not sample age (`collection_date`). The
alert reports first-*collected* (biological first appearance) and first-*visible*
(when NCBI surfaced it) as distinct fields, because a years-old sample deposited last
month is exactly the backlog artifact the emergence guard flags -- and the reader
must be able to see that gap, not have it averaged away.

The so-what priority
--------------------
Three tiers, derived from the two upstream verdicts, never invented:

  * ACT-NOW  -- a genuine (non-backlog) signal whose structural verdict predicts a
    WEAKER azole fit, measured or estimated. The dangerous combination: real spread
    AND a drug-fit consequence we can defend.
  * WATCH    -- a genuine signal whose fit consequence is minimal, uncertain, or a
    no-confident-call. Real spread, but no defensible fit story yet.
  * CONTEXT  -- a possible-backlog signal. Reported in full (we never hide it), but
    de-prioritised because the emergence guard says the spike is likely archival.

Priority is a *triage aid*, not a clinical claim; every alert still carries the raw
emergence numbers and the structural caveats so the reader can overrule it.

Two genes, one loop
-------------------
`compose_alerts` is the ERG11/azole path described above -- emergence joined to the #24
structural fit verdict. `compose_fks1_alerts` is its FKS1/echinocandin (v2) sibling and
is DETECTION-ONLY: it runs the identical detector over the `fks1_call` column but claims
NO structural verdict, because echinocandins coordinate no metal, so the CYP51/heme-iron
geometry the ERG11 so-what rests on does not transfer to a glucan synthase (see
`openafr/fks1_caller.py`). It therefore never reaches ACT-NOW (that tier requires a
weaker-fit verdict) -- only WATCH/CONTEXT -- and states the detection-only caveat on
every alert. The two share the gene-agnostic half (`_base_alert_fields`, `_rank`).

Stdlib only. Consumes openafr.emergence signals + openafr.structural verdicts; runs
fully offline against the checked-in receptor/ligand.
"""
import json

from openafr import emergence, mapping, structural

# So-what priority tiers, strongest first. `rank` orders the rendered alert list.
PRIORITY = {"act-now": 0, "watch": 1, "context": 2}


def _min_date(dates):
    """Earliest non-empty date string, or '' if none.

    ISO dates and bare years both sort correctly lexically for a minimum
    ('2019' < '2024-03-01', '2023' < '2023-05'), so no parsing is needed. We keep
    the raw string -- a bare year stays a bare year -- rather than inventing a day.
    """
    usable = [d for d in dates if (d or "").strip()]
    return min(usable) if usable else ""


def _first_seen(carriers):
    """First-collected and first-visible dates over a mutation's carriers.

    Returns {"collected", "visible"} of earliest dates. Kept separate on purpose:
    the gap between them is the backlog tell the emergence guard reasons about, and
    the reader must be able to see it rather than have it collapsed into one date.
    """
    return {
        "collected": _min_date(c.get("collection_date", "") for c in carriers),
        "visible": _min_date(c.get("target_creation_date", "") for c in carriers),
    }


def _regions(carriers):
    """Sorted distinct countries among a set of carrier records (blanks dropped)."""
    return sorted({c["country"] for c in carriers if (c.get("country") or "").strip()})


def _priority(signal, verdict):
    """Map an (emergence signal, structural verdict) pair to a so-what tier.

    See the module docstring: backlog -> context; genuine + weaker-fit -> act-now;
    genuine but no defensible fit story -> watch. Derived, never invented.
    """
    if signal["possible_backlog"]:
        return "context"
    if verdict["direction"] == "weaker-fit" and verdict["confidence"] in ("measured", "estimate"):
        return "act-now"
    return "watch"


def _emergence_kind(signal):
    """'newly appearing' / 'rising' / 'newly appearing and rising' for prose."""
    if signal["first_appearance"] and signal["rising"]:
        return "newly appearing and rising"
    if signal["first_appearance"]:
        return "newly appearing"
    return "rising"


def _headline(signal, verdict, regions, first_seen):
    """One plain sentence: mutation, spread, geography, first-seen, fit so-what."""
    mut = signal["mutation"]
    n = signal["recent_count"]
    isolate = "isolate" if n == 1 else "isolates"
    geo = (f"across {len(regions)} region(s) ({', '.join(regions)})" if regions
           else "of undisclosed geography")
    seen = first_seen["collected"] or first_seen["visible"] or "an unknown date"
    fit = {
        "weaker-fit": "predicts WEAKER azole fit",
        "minimal-direct-effect": "predicts minimal direct fit change",
        "uncertain": "has no confident fit call",
    }[verdict["direction"]]
    tag = " (possible archival backlog)" if signal["possible_backlog"] else ""
    return (f"{mut} is {_emergence_kind(signal)} on {n} {isolate} {geo}, first seen "
            f"{seen}{tag}; structurally it {fit} ({verdict['confidence']}).")


def compose_alerts(records, as_of, window_days=180, min_count=3, min_delta=0.05,
                   backlog_frac=0.5, receptor=mapping.DEFAULT_RECEPTOR,
                   ligand=mapping.DEFAULT_LIGAND, cutoff=mapping.CONTACT_CUTOFF):
    """Run the whole detect -> map -> estimate -> compose loop into ranked alerts.

    Detects emergence over `records`, joins each flagged signal to its structural
    fit verdict, enriches it with regions / first-seen / total carriers from the
    full record set, assigns a so-what priority, and returns a self-describing dict:

      as_of, window_days, thresholds, coverage -- echoed from the detector
      note                                     -- set when nothing is callable yet
      alerts                                    -- ranked list, each carrying:
        mutation, canonical, priority, headline,
        emergence   -- first_appearance/rising/backlog + the raw counts & shares
        regions, n_regions, first_seen, n_carriers_total
        structural  -- the #24 fit verdict (evidence_class, confidence, direction,
                       fluconazole_fit_verdict, lost_contacts, empirical, caveats)
        carrier_keys

    The structure is intentionally flat and JSON-serialisable: render_markdown and
    render_json both consume exactly this dict.
    """
    result = emergence.detect_emergence(
        records, as_of, window_days=window_days, min_count=min_count,
        min_delta=min_delta, backlog_frac=backlog_frac)

    out = {
        "as_of": result["as_of"], "window_days": result["window_days"],
        "thresholds": result["thresholds"], "coverage": result["coverage"],
        "alerts": [],
    }
    if "note" in result:
        out["note"] = result["note"]
    if not result["signals"]:
        return out

    by_key = {r["isolate_key"]: r for r in records}
    all_carriers = emergence.mutation_carriers(records)  # over the full set, for first-seen
    verdicts = structural.estimate_signals(result["signals"], receptor, ligand, cutoff)

    alerts = []
    for signal, verdict in zip(result["signals"], verdicts):
        base = _base_alert_fields(signal, by_key, all_carriers)
        base["priority"] = _priority(signal, verdict)
        base["headline"] = _headline(signal, verdict, base["regions"], base["first_seen"])
        base["structural"] = {
            "evidence_class": verdict["evidence_class"],
            "confidence": verdict["confidence"],
            "direction": verdict["direction"],
            "fluconazole_fit_verdict": verdict["fluconazole_fit_verdict"],
            "lost_contacts": verdict["lost_contacts"],
            "empirical": verdict["empirical"],
            "caveats": verdict["caveats"],
        }
        alerts.append(base)

    out["alerts"] = _rank(alerts)
    return out


def _base_alert_fields(signal, by_key, all_carriers):
    """The gene-agnostic half of an alert: identity, the emergence block, and the
    epidemiology (regions / first-seen / total carriers) the detector dropped.

    Shared by the ERG11 (compose_alerts) and FKS1 (compose_fks1_alerts) composers;
    each then adds its own so-what layer -- a structural fit verdict for ERG11, none
    for the detection-only FKS1 track.
    """
    recent_carriers = [by_key[k] for k in signal["carrier_keys"] if k in by_key]
    regions = _regions(recent_carriers)
    mut_carriers = all_carriers.get(signal["mutation"], [])
    first_seen = _first_seen(mut_carriers)
    return {
        "mutation": signal["mutation"],
        "canonical": signal["canonical"],
        "emergence": {
            "first_appearance": signal["first_appearance"],
            "rising": signal["rising"],
            "possible_backlog": signal["possible_backlog"],
            "backlog_reason": signal["backlog_reason"],
            "recent_count": signal["recent_count"],
            "recent_called": signal["recent_called"],
            "baseline_count": signal["baseline_count"],
            "baseline_called": signal["baseline_called"],
            "recent_prop": signal["recent_prop"],
            "baseline_prop": signal["baseline_prop"],
            "delta": signal["delta"],
        },
        "regions": regions,
        "n_regions": len(regions),
        "first_seen": first_seen,
        "n_carriers_total": len(mut_carriers),
        "carrier_keys": signal["carrier_keys"],
    }


def _rank(alerts):
    """Rank by so-what tier, then loudest emergence within a tier (keeps the
    detector's genuine-before-backlog, biggest-share-first ordering)."""
    alerts.sort(key=lambda a: (PRIORITY[a["priority"]], -a["emergence"]["delta"],
                               -a["emergence"]["recent_count"], a["mutation"]))
    return alerts


# -------------------------- FKS1 (echinocandin) -------------------------------
#
# The FKS1/echinocandin early-warning track (v2) is DETECTION-ONLY by design:
# echinocandins do not coordinate a metal, so the CYP51/heme-iron geometry that
# powers the ERG11 structural so-what does not transfer to a glucan synthase (see
# openafr/fks1_caller.py, data/earlywarning/fks1_reference/PROVENANCE.md). So an FKS1
# alert reports the *emergence* signal honestly and stops there -- no structural fit
# verdict is invented. That also means no "act-now" tier: act-now is reserved for a
# genuine spike with a measured weaker-fit structural verdict, which FKS1 never has.

# Stated on every FKS1 alert so a reader never mistakes the absence of a fit call for
# an all-clear.
FKS1_DETECTION_ONLY_CAVEAT = (
    "detection only -- no structural verdict; echinocandins do not coordinate a "
    "metal, so the CYP51/heme-iron geometry does not transfer to FKS1.")


def _fks1_priority(signal):
    """Detection-only so-what tier: backlog -> context, otherwise watch. There is no
    act-now for FKS1 because that tier requires a structural weaker-fit verdict the
    detection-only track deliberately does not claim."""
    return "context" if signal["possible_backlog"] else "watch"


def _fks1_headline(signal, regions, first_seen):
    """One plain sentence for an FKS1 signal -- mutation, spread, geography, first-seen
    -- ending in the detection-only caveat instead of a fit so-what."""
    mut = signal["mutation"]
    n = signal["recent_count"]
    isolate = "isolate" if n == 1 else "isolates"
    geo = (f"across {len(regions)} region(s) ({', '.join(regions)})" if regions
           else "of undisclosed geography")
    seen = first_seen["collected"] or first_seen["visible"] or "an unknown date"
    tag = " (possible archival backlog)" if signal["possible_backlog"] else ""
    return (f"{mut} is {_emergence_kind(signal)} on {n} {isolate} {geo}, first seen "
            f"{seen}{tag}; {FKS1_DETECTION_ONLY_CAVEAT}")


def compose_fks1_alerts(records, as_of, window_days=180, min_count=3, min_delta=0.05,
                        backlog_frac=0.5):
    """The FKS1/echinocandin analog of compose_alerts -- DETECTION-ONLY.

    Runs the same emergence detector over the `fks1_call` column instead of
    `erg11_call`, enriches each flagged signal with regions / first-seen / carriers,
    and composes alerts WITHOUT a structural fit verdict (by design; see the section
    banner above). The result dict mirrors compose_alerts' shape -- as_of, window_days,
    thresholds, coverage, optional note, ranked alerts -- with two extra markers:

      gene            -- "FKS1"
      detection_only  -- True

    and each alert carries the standard identity/emergence/epidemiology fields but no
    `structural` block; render_fks1_markdown consumes exactly this dict.
    """
    result = emergence.detect_emergence(
        records, as_of, window_days=window_days, min_count=min_count,
        min_delta=min_delta, backlog_frac=backlog_frac, call_field="fks1_call")

    out = {
        "as_of": result["as_of"], "window_days": result["window_days"],
        "thresholds": result["thresholds"], "coverage": result["coverage"],
        "gene": "FKS1", "detection_only": True, "alerts": [],
    }
    if "note" in result:
        out["note"] = result["note"]
    if not result["signals"]:
        return out

    by_key = {r["isolate_key"]: r for r in records}
    all_carriers = emergence.mutation_carriers(records, call_field="fks1_call")

    alerts = []
    for signal in result["signals"]:
        base = _base_alert_fields(signal, by_key, all_carriers)
        base["priority"] = _fks1_priority(signal)
        base["headline"] = _fks1_headline(signal, base["regions"], base["first_seen"])
        base["detection_only"] = True
        alerts.append(base)

    out["alerts"] = _rank(alerts)
    return out


# ----------------------------- rendering --------------------------------------

_BADGE = {"act-now": "ACT-NOW", "watch": "WATCH", "context": "CONTEXT"}


def render_json(alert_result, indent=2):
    """The machine-readable form: the compose_alerts dict as JSON text."""
    return json.dumps(alert_result, indent=indent, sort_keys=True)


def _fmt_first_seen(fs):
    parts = []
    if fs["collected"]:
        parts.append(f"first collected {fs['collected']}")
    if fs["visible"]:
        parts.append(f"first visible in NCBI {fs['visible']}")
    return "; ".join(parts) if parts else "first-seen date unavailable"


def render_markdown(alert_result):
    """The human form: a 'so-what in one screen' Markdown brief.

    A short header (as-of, window, coverage), a one-line priority roll-call, then one
    compact block per alert -- headline, emergence numbers, geography/first-seen, the
    structural fit verdict, and its caveats. Signal, not a data dump.
    """
    r = alert_result
    c = r["coverage"]
    lines = ["# C. auris azole early-warning alert",
             "",
             f"- **As of** {r['as_of']}  (recent window {r['window_days']}d)",
             f"- **Coverage** baseline {c['baseline_called']}/{c['baseline_total']} "
             f"called, recent {c['recent_called']}/{c['recent_total']} called, "
             f"{c['undated']} undated",
             f"- **Thresholds** {r['thresholds']}"]

    if "note" in r:
        lines += ["", f"> **No calls to assess.** {r['note']}"]
    if not r["alerts"]:
        lines += ["", "_No emergence alerts over the configured thresholds._"]
        return "\n".join(lines) + "\n"

    counts = {t: sum(1 for a in r["alerts"] if a["priority"] == t) for t in PRIORITY}
    lines += ["",
              f"**{len(r['alerts'])} alert(s):** "
              f"{counts['act-now']} act-now, {counts['watch']} watch, "
              f"{counts['context']} context.",
              ""]

    for a in r["alerts"]:
        e = a["emergence"]
        s = a["structural"]
        lines.append(f"## `{a['mutation']}` — {_BADGE[a['priority']]}")
        lines.append("")
        lines.append(a["headline"])
        lines.append("")
        kind = _emergence_kind({"first_appearance": e["first_appearance"],
                                "rising": e["rising"]})
        lines.append(
            f"- **Emergence:** {kind} — baseline {e['baseline_count']}/"
            f"{e['baseline_called']} ({e['baseline_prop']:.1%}) → recent "
            f"{e['recent_count']}/{e['recent_called']} ({e['recent_prop']:.1%}), "
            f"Δ {e['delta']:+.1%}")
        if e["possible_backlog"]:
            lines.append(f"- **Backlog check:** ⚠ {e['backlog_reason']}")
        else:
            lines.append(f"- **Backlog check:** genuine — {e['backlog_reason']}")
        geo = ", ".join(a["regions"]) if a["regions"] else "undisclosed"
        lines.append(f"- **Geography:** {a['n_regions']} region(s): {geo} "
                     f"({a['n_carriers_total']} carrier(s) total in store)")
        lines.append(f"- **First seen:** {_fmt_first_seen(a['first_seen'])}")
        lines.append(f"- **Structural fit** [{s['evidence_class']}, "
                     f"conf={s['confidence']}]: {s['fluconazole_fit_verdict']}")
        if s["lost_contacts"]:
            lc = ", ".join(f"{x['atom']}@{x['dist_to_drug']:.2f}Å"
                           for x in s["lost_contacts"])
            lines.append(f"- **Lost drug contacts (geometry):** {lc}")
        for cav in s["caveats"]:
            lines.append(f"  - _caveat:_ {cav}")
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def render_fks1_markdown(alert_result):
    """The human form for a DETECTION-ONLY FKS1/echinocandin alert.

    Same 'so-what in one screen' shape as render_markdown -- header, priority
    roll-call, one compact block per alert -- but titled for echinocandins, with no
    structural-fit line (FKS1 claims none) and a standing banner that the absence of a
    fit verdict is by design, not an all-clear.
    """
    r = alert_result
    c = r["coverage"]
    lines = ["# C. auris echinocandin (FKS1) early-warning alert",
             "",
             f"> **Detection only.** {FKS1_DETECTION_ONLY_CAVEAT}",
             "",
             f"- **As of** {r['as_of']}  (recent window {r['window_days']}d)",
             f"- **Coverage** baseline {c['baseline_called']}/{c['baseline_total']} "
             f"called, recent {c['recent_called']}/{c['recent_total']} called, "
             f"{c['undated']} undated",
             f"- **Thresholds** {r['thresholds']}"]

    if "note" in r:
        lines += ["", f"> **No calls to assess.** {r['note']}"]
    if not r["alerts"]:
        lines += ["", "_No emergence alerts over the configured thresholds._"]
        return "\n".join(lines) + "\n"

    counts = {t: sum(1 for a in r["alerts"] if a["priority"] == t) for t in PRIORITY}
    lines += ["",
              f"**{len(r['alerts'])} alert(s):** "
              f"{counts['watch']} watch, {counts['context']} context.",
              ""]

    for a in r["alerts"]:
        e = a["emergence"]
        lines.append(f"## `{a['mutation']}` — {_BADGE[a['priority']]}")
        lines.append("")
        lines.append(a["headline"])
        lines.append("")
        kind = _emergence_kind({"first_appearance": e["first_appearance"],
                                "rising": e["rising"]})
        lines.append(
            f"- **Emergence:** {kind} — baseline {e['baseline_count']}/"
            f"{e['baseline_called']} ({e['baseline_prop']:.1%}) → recent "
            f"{e['recent_count']}/{e['recent_called']} ({e['recent_prop']:.1%}), "
            f"Δ {e['delta']:+.1%}")
        if e["possible_backlog"]:
            lines.append(f"- **Backlog check:** ⚠ {e['backlog_reason']}")
        else:
            lines.append(f"- **Backlog check:** genuine — {e['backlog_reason']}")
        geo = ", ".join(a["regions"]) if a["regions"] else "undisclosed"
        lines.append(f"- **Geography:** {a['n_regions']} region(s): {geo} "
                     f"({a['n_carriers_total']} carrier(s) total in store)")
        lines.append(f"- **First seen:** {_fmt_first_seen(a['first_seen'])}")
        lines.append(f"- **Structural fit:** none claimed — {FKS1_DETECTION_ONLY_CAVEAT}")
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"
