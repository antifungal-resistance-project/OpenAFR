"""Fuse a molecule's already-computed evidence into one machine-readable candidate dossier.

Pipeline stage: report/  (runs after rank, alongside report_candidates.py)

Why this exists
---------------
A candidate's evidence is scattered across five stages that each emit their own file keyed
by molecule name: geometry rank (rank_candidates.py / analyze_poses), assay precedent
(check_precedent.py -> precedent.lookup), drug-likeness (admet.profile), human-CYP
off-target liability (cyp_offtarget.profile), and synthesizability (synth.profile).
report_candidates.py already fuses geometry + precedent into human-readable Markdown; this
module adds the missing piece — ONE machine-readable dossier per molecule that fuses all
five, assigns an A/B/C experimental priority from a documented decision table, and states
the selection hypothesis in one sentence.

What this is NOT (hard boundaries, matching the rest of the pipeline):
- It computes NOTHING new. Every field is passed in already computed; this file only
  arranges them and applies the priority table. It never re-docks and never re-runs the
  gate. A dossier is a hypothesis for a wet lab, never a hit, and means nothing unless
  validate_gate2.py has PASSED under the frozen protocol.
- It is not a ranker. Ordering stays the validated criterion (min N-Fe, unrankable last).
  Priority A/B/C is an ORTHOGONAL annotation, never a re-rank.
- It never drops a molecule. A bad profile lowers priority; it does not delete the row.

Pure and network-free by construction: callers (scripts/build_dossier.py) own all profiler
calls and the one network stage (precedent), and hand the results here.
"""

import math

from openafr import precedent as _precedent
from openafr import cyp_offtarget as _cyp

# --- priority thresholds (one place, tune once real experimental results come back) -------
TOP_FRACTION_A = 0.10        # "top decile" for automatic-A geometry
CONVERGENCE_MIN_A = 2        # >= this many pose-passing contacts to count as convergent
RO5_MAX_A = 1                # <= this many Lipinski violations for automatic A
ALERTS_MAX_A = 0             # structural-alert count allowed for automatic A (0 = clean)
NEAR_TANIMOTO = _precedent.NEIGHBOUR_TANIMOTO   # 0.85 — soft-negative neighbour cutoff

# precedent verdicts that mean "already measured active against the enzyme" — real, but a
# known inhibitor, not a discovery. Routed to a 'reference' tag, never automatic A.
_KNOWN_ACTIVE = {"tested-active"}
# verdicts that leave a novel geometric hypothesis on the table (eligible for A)
_NOVEL_OK = {"no-precedent", "tested-ambiguous", "record-inconclusive", None}


def top_rank_cutoff(n_rankable):
    """Largest geometry rank still inside the top TOP_FRACTION_A decile (1-based, >= 1)."""
    if n_rankable <= 0:
        return 0
    return max(1, math.ceil(TOP_FRACTION_A * n_rankable))


def _availability_ok(synth):
    return synth is not None and synth.get("availability") in {"catalogued-drug", "pubchem-record"}


def assign_priority(dossier, top_cutoff):
    """Return (priority, reasons) for one assembled dossier dict.

    priority is 'A' | 'B' | 'C'; reasons is a worst-first list of the individual rule hits,
    each a short human-readable string so the label is fully auditable. See module docstring
    and SPEC_candidate_dossier.md section 4 for the table.
    """
    geom = dossier["geometry"]
    prec = dossier["precedent"]
    dev = dossier["developability"]
    sel = dossier["selectivity"]
    syn = dossier["synthesis"]
    coord = dossier.get("coordinator")
    reasons = []

    # informational flags — logged so they travel with the candidate, but do not gate tier
    if sel and sel.get("azole_warhead"):
        reasons.append("azole warhead (%s) -> pan-CYP liability %s (selectivity risk)"
                       % (sel["azole_warhead"], sel.get("flags", "?")))
    prec_verdict = prec.get("verdict") if prec else None
    if prec_verdict in _KNOWN_ACTIVE:
        reasons.append("reference: already measured active against the enzyme (positive control, not a novel hit)")
    # Coordinator identity (RESULTS_coordinator.md / RESULTS_aromatic_criterion.md): the rank is
    # only azole-like if the nitrogen reaching the iron is the validated aromatic ring N. Both a
    # dual-mode and an out-of-mechanism rank BLOCK automatic-A (see the A-conditions) and are
    # flagged here — but neither forces C. The pre-registered gate re-validation showed the active
    # luliconazole itself coordinates via its nitrile in-pose (aromatic N 5.3 A), so a non-aromatic
    # coordinator is NOT proof of an artifact: it is unvalidated, not disproven. Hard-demoting it to
    # C would risk rejecting a real active, so the coordinator axis caps priority at B, never C.
    cv = coord.get("verdict") if coord else None
    if cv == "dual-mode":
        ad = coord.get("aromatic_distance")
        reasons.append("dual-mode coordination: reported rank set by a %s N; the aromatic "
                       "ring N reaches %s A (azole-like mode available, but the rank is inflated)"
                       % (coord.get("reported_kind"),
                          "%.3f" % ad if ad is not None else "?"))
    elif cv == "out-of-mechanism":
        reasons.append("out-of-mechanism coordination: reported rank set by a %s N, not the "
                       "validated aromatic ring N, and no aromatic N reaches the actives' band "
                       "— the geometric evidence is unvalidated (blocks A, not disproven: the "
                       "active luliconazole also coordinates via a non-aromatic N)"
                       % coord.get("reported_kind"))

    # --- automatic C (deprioritise: don't spend money, or unreliable readout) --------------
    demote = False
    if prec_verdict == "tested-inactive":
        reasons.append("tested-inactive: enzyme measured inactive — don't re-run a lost experiment")
        demote = True
    near = dossier.get("_near")
    if near and near.get("verdict") == "tested-inactive" and near.get("tanimoto", 0) >= NEAR_TANIMOTO:
        reasons.append("soft negative: near-identical compound (Tanimoto %.2f) measured inactive"
                       % near["tanimoto"])
        demote = True
    if geom.get("best_fe") is None:
        reasons.append("unrankable: no usable pose — no geometric hypothesis to test")
        demote = True
    if dev and dev.get("pains"):
        reasons.append("PAINS alert: assay-interference risk makes the readout unreliable")
        demote = True
    if dossier.get("error"):
        reasons.append("profile error: %s" % dossier["error"])
        demote = True
    if demote:
        return "C", reasons

    # --- automatic A (test first) requires ALL of the following ----------------------------
    a_conditions = {
        "top-decile geometry rank": geom.get("rank") is not None and top_cutoff >= 1
            and geom["rank"] <= top_cutoff,
        "convergent poses": (geom.get("n_passing") or 0) >= CONVERGENCE_MIN_A,
        "not a known active": prec_verdict in _NOVEL_OK,
        "clean drug-likeness": dev is not None
            and dev.get("ro5_violations", 99) <= RO5_MAX_A
            and dev.get("n_alerts", 99) <= ALERTS_MAX_A,
        "physically obtainable": _availability_ok(syn),
        # the rank must reflect the validated azole mechanism (an aromatic ring N reaching
        # the iron), not a nitrile/amine/azide artifact. Unsupplied coordinator == pass, so
        # a dossier built without pose access is unchanged (RESULTS_coordinator.md).
        "in-mechanism coordination": coord is None or coord.get("verdict") == "in-mechanism",
    }
    if all(a_conditions.values()):
        reasons.insert(0, "priority A: " + ", ".join(sorted(a_conditions)))
        return "A", reasons

    # --- everything else is B, with the failed A-conditions named as the caveat ------------
    missing = [k for k, ok in a_conditions.items() if not ok]
    reasons.append("priority B: not automatic-A (%s)" % ", ".join(sorted(missing)))
    return "B", reasons


def selection_hypothesis(dossier):
    """One-sentence, honest 'reason for selection' assembled from the dossier fields."""
    geom = dossier["geometry"]
    prec = dossier["precedent"]
    sel = dossier["selectivity"]
    if geom.get("best_fe") is None:
        return ("No usable docked pose — carried for completeness (the pipeline ranks "
                "unrankable molecules last, never drops them), not a live hypothesis.")
    rank_bit = "ranks #%s by N-Fe geometry" % geom.get("rank", "?")
    conv_bit = ("convergent poses (%d/%d contact the iron)"
                % (geom.get("n_passing", 0), geom.get("n_poses", 0)))
    verdict = prec.get("verdict") if prec else None
    if verdict in _KNOWN_ACTIVE:
        prec_bit = "already a measured enzyme active (reference/positive control)"
    elif verdict == "no-precedent":
        prec_bit = "no prior fungal-CYP51 measurement (novel, but absence is not proof of untested)"
    elif verdict is None:
        prec_bit = "precedent not checked"
    else:
        prec_bit = "prior-assay verdict '%s'" % verdict
    tail = ""
    if sel and sel.get("azole_warhead"):
        tail = "; the %s warhead is the selectivity liability to watch" % sel["azole_warhead"]
    coord = dossier.get("coordinator")
    if coord and coord.get("verdict") == "out-of-mechanism":
        tail += ("; WARNING the rank is set by a %s nitrogen, not the validated aromatic "
                 "ring N — outside the method's mechanism" % coord.get("reported_kind"))
    elif coord and coord.get("verdict") == "dual-mode":
        ad = coord.get("aromatic_distance")
        tail += ("; the rank is inflated by a %s N (aromatic ring N at %s A)"
                 % (coord.get("reported_kind"), "%.3f" % ad if ad is not None else "?"))
    return "%s with %s, %s%s." % (rank_bit.capitalize(), conv_bit, prec_bit, tail)


def build_one(name, smiles, geometry, *, admet=None, cyp=None, synth=None,
              prec=None, near=None, coord=None, error=None, top_cutoff=0):
    """Assemble + prioritise ONE dossier from already-computed pieces.

    geometry: {rank, best_fe, n_passing, n_poses, best_score}. admet/cyp/synth: the
    respective profile() dicts, or None if that profiler failed/was skipped. prec:
    precedent.lookup() dict or None (network stage; None == not checked). near: nearest
    tested-neighbour dict {tanimoto, name, verdict} or None. error: per-molecule failure
    string (unparseable SMILES, etc.) — recorded, priority forced to C.
    """
    inchikey = prec.get("inchikey") if prec else _precedent.inchikey(smiles)

    selectivity = None
    if cyp is not None:
        selectivity = {
            "azole_warhead": cyp.get("azole_warhead"),
            "flags": cyp.get("flags") or _cyp.flag_summary(cyp),
            "n_high": cyp.get("n_high"), "n_moderate": cyp.get("n_moderate"),
            "flagged_isoforms": cyp.get("flagged_isoforms", []),
        }

    developability = None
    if admet is not None:
        developability = {
            "mw": admet.get("mw"), "clogp": admet.get("clogp"), "tpsa": admet.get("tpsa"),
            "ro5_violations": admet.get("ro5_violations"), "veber_ok": admet.get("veber_ok"),
            "esol_logs": admet.get("esol_logs"),
            "pains": admet.get("pains", []), "structural_alerts": admet.get("structural_alerts", []),
            "n_alerts": admet.get("n_alerts"), "herg_risk_factors": admet.get("herg_risk_factors"),
        }

    synthesis = None
    if synth is not None:
        synthesis = {
            "sa_score": synth.get("sa_score"), "sa_band": synth.get("sa_band"),
            "availability": synth.get("availability"), "pubchem_cid": synth.get("pubchem_cid"),
        }

    precedent_block = None
    if prec is not None and prec.get("verdict") not in (None, "not-in-library"):
        precedent_block = {
            "verdict": prec.get("verdict"),
            "phenotypic": prec.get("phenotypic"),
            "n_off_target": prec.get("n_off_target"),
            "nearest_tested": near,
        }

    d = {
        "name": name,
        "smiles": smiles,
        "inchikey": inchikey,
        "geometry": geometry,
        "precedent": precedent_block,
        "developability": developability,
        "selectivity": selectivity,
        "synthesis": synthesis,
        "coordinator": coord,   # {reported_kind, in_mechanism, aromatic_distance, verdict} or None
        "_near": near,          # transient input for the priority table; stripped below
    }
    if error:
        d["error"] = error

    priority, reasons = assign_priority(d, top_cutoff)
    d["priority"] = priority
    d["priority_reasons"] = reasons
    d["selection_hypothesis"] = selection_hypothesis(d)
    d.pop("_near", None)
    return d


def build_all(records, provenance):
    """Assemble a whole ranked batch into dossiers, preserving geometry-rank order.

    records: an iterable of dicts, each with the build_one keyword pieces plus 'name',
    'smiles', 'geometry' (already in rank order — rank field 1-based over rankable rows).
    provenance: the stamped-once block copied onto every dossier. Returns the list of
    dossiers in the same order they came in (i.e. geometry rank).
    """
    records = list(records)
    n_rankable = sum(1 for r in records if r["geometry"].get("best_fe") is not None)
    cutoff = top_rank_cutoff(n_rankable)
    out = []
    for r in records:
        d = build_one(
            r["name"], r["smiles"], r["geometry"],
            admet=r.get("admet"), cyp=r.get("cyp"), synth=r.get("synth"),
            prec=r.get("prec"), near=r.get("near"), coord=r.get("coord"),
            error=r.get("error"), top_cutoff=cutoff,
        )
        d["provenance"] = provenance
        out.append(d)
    return out
