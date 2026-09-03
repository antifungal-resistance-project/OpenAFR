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
  work/candidates/shortlist_resistance_K143R.json  C. auris K143R retention, modeled rotamer (#77)
  work/candidates/shortlist_resistance_F126L.json  C. auris F126L retention, modeled rotamer (#77)
  work/candidates/pose_convergence.json       within-search pose convergence (#71)
  work/candidates/shortlist_stereo.json       stereochemistry / isomer audit (#85)
  work/candidates/shortlist_protomer.json     protonation/tautomer microstate audit (#84)
  work/candidates/shortlist_ensemble.json     crystallographic ensemble consistency (#70)
  work/candidates/shortlist_cyp.json          human-CYP off-target liability panel (#74)

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


def load_resistance_addatom():
    """{name: {mutant: verdict-row}} for the add-atom C. auris mutants K143R/F126L (#77).

    Each mutant is a separate frozen per-axis file; a missing file just leaves that pocket
    not-yet-checked, so the scorecard degrades honestly rather than erroring."""
    out = {}
    for mutant in ("K143R", "F126L"):
        path = os.path.join(C, f"shortlist_resistance_{mutant}.json")
        if not os.path.exists(path):
            continue
        for x in json.load(open(path))["candidates"]:
            out.setdefault(x["name"], {})[mutant] = x
    return out


def load_convergence():
    path = os.path.join(C, "pose_convergence.json")
    if not os.path.exists(path):
        return {}
    return {x["name"]: x for x in json.load(open(path))["candidates"]}


def load_stereo():
    path = os.path.join(C, "shortlist_stereo.json")
    if not os.path.exists(path):
        return {}
    return {x["name"]: x for x in json.load(open(path))["candidates"]}


def load_protomer():
    path = os.path.join(C, "shortlist_protomer.json")
    if not os.path.exists(path):
        return {}
    return {x["name"]: x for x in json.load(open(path))["candidates"]}


def load_ensemble():
    path = os.path.join(C, "shortlist_ensemble.json")
    if not os.path.exists(path):
        return {}
    return {x["name"]: x for x in json.load(open(path))["candidates"]}


def load_cyp():
    path = os.path.join(C, "shortlist_cyp.json")
    if not os.path.exists(path):
        return {}
    return {x["name"]: x for x in json.load(open(path))["candidates"]}


def load_coordinator():
    path = os.path.join(C, "coordinator_identity.json")
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
    resistance_addatom = load_resistance_addatom()
    convergence = load_convergence()
    stereo = load_stereo()
    protomer = load_protomer()
    ensemble = load_ensemble()
    cyp = load_cyp()
    coordinator = load_coordinator()

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
        rsa = resistance_addatom.get(name, {})
        pc = convergence.get(name, {})
        st = stereo.get(name, {})
        pt = protomer.get(name, {})
        en = ensemble.get(name, {})
        cy = cyp.get(name, {})
        co = coordinator.get(name, {})
        is_ctrl = name in CONTROLS

        concerns = derive_concerns(cf, am, sy, dm, pr, is_ctrl, rs or None,
                                   st or None, pt or None, en or None, cy or None,
                                   rsa or None, co or None)
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
            "resistance_retention_addatom": {
                "status": "proven" if rsa else "not-yet-checked",
                "note": "MODELED-ROTAMER pockets (side-chain repacker; minimum-strain "
                        "rotamer) — a weaker structural claim than the Y132F deletion",
                "mutants": {
                    m: {
                        "clade": {"K143R": "I", "F126L": "III"}[m],
                        "verdict": rsa[m].get("verdict"),
                        "consensus_nfe_A": rsa[m].get("mutant_consensus"),
                        "seeds_coordinating": rsa[m].get("mutant_seeds_coordinating"),
                        "retention_shift_A": rsa[m].get("retention_shift"),
                    }
                    for m in ("K143R", "F126L") if m in rsa
                },
            },
            "pose_convergence": {
                "status": "proven" if pc else "not-yet-checked",
                "verdict": pc.get("verdict"),
                "top_cluster_fraction": pc.get("top_cluster_fraction"),
                "n_clusters": pc.get("n_clusters"),
                "coordinating_pose_rank": pc.get("coordinating_pose_rank"),
                "coordinating_pose_in_top_cluster": pc.get("coordinating_pose_in_top_cluster"),
            },
            "stereochemistry": {
                "status": "proven" if st else "not-yet-checked",
                "verdict": st.get("verdict"),
                "stereo_explicit": st.get("stereo_explicit"),
                "n_isomers": st.get("n_isomers"),
                "n_unassigned_centers": st.get("n_unassigned"),
                "iso_spread_A": st.get("iso_spread"),
            },
            "protomer_tautomer": {
                "status": "proven" if pt else "not-yet-checked",
                "verdict": pt.get("verdict"),
                "protomer_explicit": pt.get("protomer_explicit"),
                "n_states": pt.get("n_states"),
                "n_protonation_states": pt.get("n_protonation_states"),
                "state_spread_A": pt.get("state_spread"),
            },
            "ensemble_consistency": {
                "status": "proven" if en else "not-yet-checked",
                "ensemble": "C. albicans CYP51 {5TZ1, 5FSA, 5V5Z}",
                "verdict": en.get("verdict"),
                "ensemble_min_nfe_A": en.get("ensemble_min"),
                "ensemble_mean_nfe_A": en.get("ensemble_mean"),
                "conformers_coordinating": en.get("n_conf_coord"),
            },
            "cyp_offtarget": {
                "status": "proven" if cy else "not-yet-checked",
                "panel": "human CYP3A4/2C9/2D6/2C19/1A2",
                "azole_warhead": cy.get("azole_warhead"),
                "type_ii_ligator": cy.get("type_ii_ligator"),
                "n_high": cy.get("n_high"),
                "n_moderate": cy.get("n_moderate"),
                "flags": cy.get("flags"),
                "liability": ({k: v["liability"] for k, v in cy.get("panel", {}).items()}
                              or None),
            },
            "coordinator_identity": {
                "status": "proven" if co else "not-yet-checked",
                "verdict": co.get("verdict"),
                "reported_kind": co.get("reported_kind"),
                "in_mechanism": co.get("reported_in_mechanism"),
                "aromatic_nfe_A": co.get("in_mechanism_distance"),
                "reported_nfe_A": co.get("reported_distance"),
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
            "resistance_retention_addatom (C. auris K143R/F126L, modeled-rotamer repack)",
            "pose_convergence (within-search RMSD clustering)",
            "stereochemistry (supplied-SMILES isomer audit + isomer dock)",
            "protomer_tautomer (pH-7.4 protonation/tautomer microstate dock)",
            "ensemble_consistency (crystallographic C. albicans ensemble dock)",
            "cyp_offtarget (human CYP3A4/2C9/2D6/2C19/1A2 type-II + pharmacophore liability)",
            "coordinator_identity (which nitrogen sets the rank: aromatic ring vs nitrile/amine/azide)",
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


def res_addatom_cell(c):
    """Compact K143R/F126L (modeled-rotamer) retention: one initial + verdict per mutant.

    e.g. 'K:RET F:LOST'. RET/SHF/LOST keep the column narrow; the full verdict + numbers live
    in resistance_retention_addatom. '—' until measured."""
    ra = c["resistance_retention_addatom"]
    if ra["status"] != "proven":
        return "—"
    short = {"RETAINED": "RET", "LOST": "LOST", "SHIFTED-away": "SHFa",
             "SHIFTED-closer": "SHFc", "UNRANKABLE": "?"}
    return " ".join("%s:%s" % (m[0], short.get(ra["mutants"][m]["verdict"],
                                              (ra["mutants"][m]["verdict"] or "—")[:4]))
                    for m in ("K143R", "F126L") if m in ra["mutants"])


def conv_cell(c):
    """Compact within-search pose-convergence cell: verdict + coordinating-pose rank.

    Every candidate is DISPERSED under the frozen big-box protocol, so the discriminating
    signal is the coordinating pose's affinity rank (rk1 = the reported geometry is also the
    best-energy pose; ✓ = it sits in the dominant cluster)."""
    pc = c["pose_convergence"]
    if pc["status"] != "proven":
        return "—"
    v = "conv" if pc["verdict"] == "CONVERGED" else "disp"
    rk = pc["coordinating_pose_rank"]
    tick = "✓" if pc["coordinating_pose_in_top_cluster"] else ""
    return f"{v}/rk{rk}{tick}" if rk is not None else v


def stereo_cell(c):
    """Compact stereochemistry cell: EXPLICIT (single defined isomer) or the isomer verdict."""
    st = c["stereochemistry"]
    if st["status"] != "proven":
        return "—"
    if st.get("stereo_explicit"):
        return "explicit"
    v = st.get("verdict")
    return f"{'robust' if v == 'ROBUST' else 'iso-dep'}/{st.get('n_isomers')}iso"


def protomer_cell(c):
    """Compact protonation/tautomer cell: EXPLICIT (single microstate) or the microstate verdict."""
    pt = c["protomer_tautomer"]
    if pt["status"] != "proven":
        return "—"
    if pt.get("protomer_explicit"):
        return "explicit"
    v = pt.get("verdict")
    tag = "robust" if v == "ROBUST" else "state-dep"
    return f"{tag}/{pt.get('n_states')}st"


def ens_cell(c):
    """Compact ensemble-consistency cell: verdict + conformers coordinating.

    ROBUST = coordinates the iron in all 3 crystallographic conformers; frag = only 1
    (single-snapshot artifact). ens-min in Angstrom is in the full JSON card."""
    en = c["ensemble_consistency"]
    if en["status"] != "proven":
        return "—"
    tag = {"ROBUST": "robust", "CONSISTENT": "consist", "CONFORMER-FRAGILE": "FRAGILE",
           "NON-COORDINATING": "NONE", "UNRANKABLE": "unrank"}.get(en["verdict"], en["verdict"])
    return f"{tag}/{en['conformers_coordinating']}"


def cyp_cell(c):
    """Compact human-CYP off-target cell: type-II HIGH pan-panel, or the moderate profile.

    'II:HIGHx5' == an azole warhead flags all five human CYPs HIGH (the same heme ligation
    the on-target criterion selected for). 'modx5' == a weaker azine/pharmacophore profile.
    The per-isoform liability + basis live in the full JSON card."""
    cy = c["cyp_offtarget"]
    if cy["status"] != "proven":
        return "—"
    if cy.get("n_high"):
        return "II:HIGHx%d" % cy["n_high"]
    if cy.get("n_moderate"):
        return "modx%d" % cy["n_moderate"]
    return "clean"


def coord_cell(c):
    """Compact coordinator-identity cell: which nitrogen sets the rank.

    'in' = aromatic ring N (validated mechanism); 'dual/<kind>' = non-aromatic N sets the rank
    but the aromatic N still reaches the band; 'OUT/<kind>' = non-aromatic N and no aromatic N
    in band (the rank is an artifact). Full distances live in the JSON card."""
    co = c["coordinator_identity"]
    if co["status"] != "proven":
        return "—"
    v = co.get("verdict")
    if v == "in-mechanism":
        return "in"
    kind = (co.get("reported_kind") or "?")[:6]
    return ("dual/%s" % kind) if v == "dual-mode" else ("OUT/%s" % kind
            if v == "out-of-mechanism" else (v or "—"))


def print_table(sc):
    novel = [c for c in sc["candidates"] if not c["is_control"]]
    ctrl = [c for c in sc["candidates"] if c["is_control"]]
    cols = ("candidate", "N-Fe", "seeds", "sel-margin", "admet", "SA",
            "domain", "Y132F", "K143R/F126L", "ensemble", "pose-conv", "stereo", "protomer",
            "hCYP", "coord", "precedent", "#concern")
    w = ("%-14s %-6s %-6s %-10s %-14s %-9s %-18s %-8s %-13s %-11s %-10s %-11s %-12s "
         "%-10s %-12s %-16s %s")
    print("\nNOVEL CANDIDATES (geometry-ranked)")
    print(w % cols)
    print("-" * 130)
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
            res_addatom_cell(c),
            ens_cell(c),
            conv_cell(c),
            stereo_cell(c),
            protomer_cell(c),
            cyp_cell(c),
            coord_cell(c),
            (c["precedent"]["phenotypic"] if c["precedent"]["phenotypic"] != "—"
             else c["precedent"]["enzyme_verdict"])[:16],
            len(c["concerns"]),
        ))
    print("\nCONTROLS (known azoles — validate the axes, not discoveries)")
    print(w % cols)
    print("-" * 130)
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
            res_addatom_cell(c),
            ens_cell(c),
            conv_cell(c),
            stereo_cell(c),
            protomer_cell(c),
            cyp_cell(c),
            coord_cell(c),
            (c["precedent"]["enzyme_verdict"])[:16],
            len(c["concerns"]),
        ))
    print(f"\nAxes NOT yet checked for any candidate: "
          f"{', '.join(sc['axes_not_yet_checked'].values())}")


if __name__ == "__main__":
    sc = build()
    # --out lets the reproducibility harness (scripts/repro.py) regenerate to a scratch
    # path and byte-compare against the committed card without clobbering it.
    out = (sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv
           else os.path.join(C, "shortlist_scorecard.json"))
    json.dump(sc, open(out, "w"), indent=2)
    print(f"wrote {out}  ({sc['n_candidates']} candidates, "
          f"{sc['n_controls']} controls)")
    print_table(sc)
