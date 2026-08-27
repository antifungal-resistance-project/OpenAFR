"""Consolidated per-candidate confidence scorecard (issue #88).

Rolls every confidence axis that has already been measured into ONE honest,
at-a-glance card per candidate — the artifact you actually reason over before
choosing what to send to a wet lab. It does NOT re-dock or re-measure anything:
it merges the frozen per-axis outputs and, crucially, marks the axes that have
*not* been checked for a candidate as `not-yet-checked` rather than passing them
over in silence. No composite score is invented — a single number would hide the
exact thing this card exists to show: where a molecule is strong vs unproven.

Inputs (all produced by earlier, separately-validated runs):
  work/repurposing/ranked.tsv                 single-dock geometry rank (#screen)
  work/candidates/shortlist_confidence.json   5-seed stability + human selectivity (#68,#73)
  work/candidates/shortlist_admet.tsv         ADMET / structural-alert flags (#75)
  work/candidates/shortlist_synth.json        synthesizability + availability (#86)
  work/candidates/domain_sensitivity.json     applicability-domain fragility (#87)
  work/candidates/top100_precedent.tsv        assay-precedent verdicts (#65)
  work/candidates/shortlist_resistance.json   C. auris Y132F retention (#69,#77)

Usage:  conda run -n openafr python work/candidates/build_scorecard.py
Writes: work/candidates/shortlist_scorecard.json  and prints a readable table.
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from scripts.scorecard import CONTROLS, FE_CUTOFF_A, NOT_YET_CHECKED, derive_concerns

C = os.path.join(ROOT, "work", "candidates")
FE_CUTOFF = FE_CUTOFF_A


def load_ranked():
    d = {}
    with open(os.path.join(ROOT, "work", "repurposing", "ranked.tsv")) as f:
        next(f)
        for line in f:
            rank, lig, nfe, vina = line.rstrip("\n").split("\t")
            d[lig] = {"rank": int(rank), "single_dock_nfe": float(nfe),
                      "vina_kcal_mol": float(vina)}
    return d


def load_confidence():
    data = json.load(open(os.path.join(C, "shortlist_confidence.json")))
    return {x["name"]: x for x in data["candidates"]}


def load_tsv(path, key="name"):
    d = {}
    with open(path) as f:
        hdr = next(f).rstrip("\n").split("\t")
        for line in f:
            row = dict(zip(hdr, line.rstrip("\n").split("\t")))
            d[row[key]] = row
    return d


def load_synth():
    return {x["name"]: x for x in json.load(open(os.path.join(C, "shortlist_synth.json")))}


def load_domain():
    data = json.load(open(os.path.join(C, "domain_sensitivity.json")))
    return {x["name"]: x for x in data["candidate_fragility"]}, data


def load_resistance():
    path = os.path.join(C, "shortlist_resistance.json")
    if not os.path.exists(path):
        return {}
    return {x["name"]: x for x in json.load(open(path))["candidates"]}


def build():
    ranked = load_ranked()
    conf = load_confidence()
    admet = load_tsv(os.path.join(C, "shortlist_admet.tsv"))
    synth = load_synth()
    domain, domain_meta = load_domain()
    prec = load_tsv(os.path.join(C, "top100_precedent.tsv"))
    resistance = load_resistance()

    # Candidate universe: everything the stability/selectivity run scored, in its order.
    order = [x["name"] for x in json.load(open(os.path.join(C, "shortlist_confidence.json")))["candidates"]]

    cards = []
    for name in order:
        cf = conf.get(name, {})
        am = admet.get(name, {})
        sy = synth.get(name, {})
        dm = domain.get(name, {})
        pr = prec.get(name, {})
        rs = resistance.get(name, {})
        is_ctrl = name in CONTROLS

        concerns = derive_concerns(cf, am, sy, dm, pr, is_ctrl, rs or None)
        fungal = cf.get("fungal_consensus")
        reaches = cf.get("reaches_fungal_iron")
        stable = cf.get("stable")
        margin = cf.get("selectivity_margin")
        admet_flags = am.get("flags", "")
        sa_band = sy.get("sa_band")
        dom_flags = dm.get("flags", [])
        enzyme = pr.get("verdict", "—")
        pheno = pr.get("phenotypic", "—")

        card = {
            "name": name,
            "is_control": is_ctrl,
            "geometry": {
                "status": "proven",
                "single_dock_rank": ranked.get(name, {}).get("rank"),
                "single_dock_nfe_A": ranked.get(name, {}).get("single_dock_nfe"),
                "consensus_nfe_A": fungal,
                "reaches_fungal_iron": reaches,
            },
            "pose_stability": {
                "status": "proven" if cf else "not-yet-checked",
                "seeds_coordinating": cf.get("fungal_seeds_coordinating"),
                "spread_A": cf.get("fungal_spread"),
                "stable": stable,
            },
            "human_selectivity": {
                "status": "proven" if margin is not None else "not-yet-checked",
                "human_consensus_nfe_A": cf.get("human_consensus"),
                "reaches_human_iron": cf.get("reaches_human_iron"),
                "selectivity_margin_A": margin,
            },
            "admet": {
                "status": "proven" if am else "not-yet-checked",
                "flags": admet_flags,
                "ro5_violations": am.get("ro5_violations"),
                "esol_logs": am.get("esol_logs"),
                "n_pains": am.get("n_pains"),
            },
            "synthesizability": {
                "status": "proven" if sy else "not-yet-checked",
                "sa_score": sy.get("sa_score"),
                "sa_band": sa_band,
                "availability": sy.get("availability"),
            },
            "applicability_domain": {
                "status": "proven" if dm else "not-yet-checked",
                "heavy_atoms": dm.get("heavy"),
                "aromatic_n": dm.get("aromatic_n"),
                "flags": dom_flags,
            },
            "precedent": {
                "status": "proven" if pr else "not-yet-checked",
                "enzyme_verdict": enzyme,
                "phenotypic": pheno,
                "chembl_id": pr.get("chembl_id"),
                "near_verdict": pr.get("near_verdict", "—"),
            },
            "resistance_retention": {
                "status": "proven" if rs else "not-yet-checked",
                "pocket": "C. auris Y132F",
                "verdict": rs.get("verdict"),
                "y132f_consensus_nfe_A": rs.get("y132f_consensus"),
                "seeds_coordinating": rs.get("y132f_seeds_coordinating"),
                "retention_shift_A": rs.get("retention_shift"),
            },
            "not_yet_checked": dict(NOT_YET_CHECKED),
            "concerns": concerns,
        }
        cards.append(card)

    return {
        "fe_cutoff_A": FE_CUTOFF,
        "n_candidates": len(cards),
        "n_controls": sum(c["is_control"] for c in cards),
        "domain_thresholds": domain_meta["frozen_thresholds"],
        "axes_measured": [
            "geometry (single-dock + 5-seed consensus)",
            "pose_stability (5-seed)",
            "human_CYP51_selectivity",
            "admet",
            "synthesizability",
            "applicability_domain",
            "assay_precedent",
            "resistance_retention (C. auris Y132F)",
        ],
        "axes_not_yet_checked": NOT_YET_CHECKED,
        "candidates": cards,
    }


def fmt(v, nd=3):
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.{nd}f}"
    return str(v)


def res_cell(c):
    """Compact C. auris Y132F retention verdict for the table."""
    return (c["resistance_retention"]["verdict"] or "—")[:8]


def print_table(sc):
    novel = [c for c in sc["candidates"] if not c["is_control"]]
    ctrl = [c for c in sc["candidates"] if c["is_control"]]
    cols = ("candidate", "N-Fe", "seeds", "sel-margin", "admet", "SA",
            "domain", "Y132F", "precedent", "#concern")
    w = "%-14s %-6s %-6s %-10s %-14s %-9s %-18s %-8s %-16s %s"
    print("\nNOVEL CANDIDATES (geometry-ranked)")
    print(w % cols)
    print("-" * 118)
    for c in novel:
        print(w % (
            c["name"][:14],
            fmt(c["geometry"]["consensus_nfe_A"]),
            c["pose_stability"]["seeds_coordinating"] or "—",
            fmt(c["human_selectivity"]["selectivity_margin_A"], 2),
            (c["admet"]["flags"] or "—")[:14],
            f"{c['synthesizability']['sa_band']}",
            (",".join(c["applicability_domain"]["flags"]) or "clean")[:18],
            res_cell(c),
            (c["precedent"]["phenotypic"] if c["precedent"]["phenotypic"] != "—"
             else c["precedent"]["enzyme_verdict"])[:16],
            len(c["concerns"]),
        ))
    print("\nCONTROLS (known azoles — validate the axes, not discoveries)")
    print(w % cols)
    print("-" * 118)
    for c in ctrl:
        print(w % (
            c["name"][:14],
            fmt(c["geometry"]["consensus_nfe_A"]),
            c["pose_stability"]["seeds_coordinating"] or "—",
            fmt(c["human_selectivity"]["selectivity_margin_A"], 2),
            (c["admet"]["flags"] or "—")[:14],
            f"{c['synthesizability']['sa_band']}",
            (",".join(c["applicability_domain"]["flags"]) or "clean")[:18],
            res_cell(c),
            (c["precedent"]["enzyme_verdict"])[:16],
            len(c["concerns"]),
        ))
    print(f"\nAxes NOT yet checked for any candidate: "
          f"{', '.join(sc['axes_not_yet_checked'].values())}")


if __name__ == "__main__":
    sc = build()
    out = os.path.join(C, "shortlist_scorecard.json")
    json.dump(sc, open(out, "w"), indent=2)
    print(f"wrote {out}  ({sc['n_candidates']} candidates, "
          f"{sc['n_controls']} controls)")
    print_table(sc)
