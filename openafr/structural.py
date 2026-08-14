"""Structural so-what: the azole-fit consequence of a flagged mutation (issue #24).

Pipeline stage: the actionable tip of the C. auris azole early-warning loop.

Why this exists
---------------
The emergence detector (#22) says a mutation is *new* or *rising*; the mapper (#23)
says *where* it sits in the modeled azole pocket and whether we can build it. Neither
answers the question a clinician actually asks:

  If this mutation spreads, does fluconazole still fit the pocket -- and how much worse?

This module turns a mapping verdict into a fit-consequence verdict. It is deliberately
narrow about what it can honestly claim, because a wrong "the drug still works" is the
most dangerous output the whole track could emit.

The three evidence tiers (strongest to weakest)
-----------------------------------------------
1. EMPIRICAL-DOCK. For a mutation we have *actually run* through the #13 C. auris port
   pipeline (build the mutant, dock the held-out azoles, grade the geometric gate), we
   report the *measured, pre-registered* outcome -- not a re-run. Right now that is
   Y132F alone (EMPIRICAL_DOCKS). This is the only tier that carries a real number, and
   it still ships every caveat the #13 result declared.

2. GEOMETRIC-ESTIMATE. For a deletion-buildable mutation we have not yet docked, we give
   a *direction*, never a number: a deletion at a residue whose removed atoms contact
   the drug takes favorable contacts away -> predicted WEAKER fit; a deletion that only
   removes second-shell atoms -> minimal predicted *direct* effect. The magnitude waits
   on the dock, and we print the exact command to get it. The estimate is deterministic
   (the removed atoms are fixed, their distance to the co-crystal azole is measured).

3. NO-CALL. For a mutation the pinned env cannot even build -- an add-atom substitution
   needing a side-chain repacker (K143R, F126L) -- or a token that does not map to a
   residue at all (TR34), we refuse to call the fit. A pocket-contact unbuildable
   mutation is flagged as *plausible but unresolved*; a second-shell one as *low prior*.
   We never invent a ΔΔG we cannot compute.

Honesty bar (same as the main repo)
-----------------------------------
* Direction without a dock is labeled an estimate, never a measurement.
* "The drug still fits" is the claim we are most reluctant to make: even the empirical
  Y132F tier reports that broad iron-reach *survives* while the top-rank *advantage*
  collapses, and flags that #13 docked an azole *class*, not fluconazole itself.
* Every branch resolves to measured | estimate | none -- there is always an explicit
  confidence, and "no confident call" is a first-class, expected outcome.

Stdlib only. Consumes openafr.mapping verdicts + the checked-in work/receptor_A.pdb and
work/ref_VT1.pdb; runs fully offline.
"""
from openafr import mapping

# Atoms a deterministic deletion mutation removes from the wild-type side chain. MUST
# mirror scripts/mutate_receptor.py SUPPORTED and openafr.mapping.DELETION_MUTATIONS --
# test_structural asserts all three agree, so they cannot silently drift. The atom names
# (not just the residue pair) are what lets us measure exactly which pocket contacts a
# deletion takes away.
DELETION_ATOMS = {
    ("TYR", "PHE"): {"OH"},          # Phe = Tyr minus the para-hydroxyl
    ("VAL", "ALA"): {"CG1", "CG2"},  # Ala = Val minus the two gamma methyls
}

# Mutations that have actually been run through the #13 port pipeline, with the measured,
# pre-registered outcome. This is the honest "reuse #13": we record the graded result and
# its provenance rather than re-running a multi-hour, pre-registered docking exercise. Add
# a mutation here only when a real, pre-registered dock exists for it.
EMPIRICAL_DOCKS = {
    "Y132F": {
        "source": "work/RESULTS_auris_Y132F.md",
        "prereg_sha256": "689de2fd",       # PREREGISTRATION_auris_Y132F.md, verified at grading
        "protocol_sha256": "6eca16b1",     # frozen docking protocol, verified unmodified
        "mutant_sha256": "924cf7b8",       # work/receptor_auris_Y132F.pdb (rigid OH deletion)
        "gate": "FAIL",                    # bar: AUC >= 0.70 AND EF@1% >= 5.0x
        "auc": 0.805, "auc_wildtype": 0.794,
        "ef1_wildtype": 12.68, "ef1_mutant": 0.00, "ef5_mutant": 5.63,
        "ligand_set": "7 held-out azole actives (class) + property-matched decoys",
        # The texture the number alone would hide (verbatim intent of RESULTS_auris_Y132F.md):
        "direction": "weaker-fit",
        "finding": ("broad iron-reach separation survives (AUC 0.805 vs wild-type 0.794) -- "
                    "azoles as a class still coordinate the iron -- but the top-rank fit "
                    "advantage collapses: property-matched decoys tie and edge past the best "
                    "real drugs at the extreme tip, dropping EF@1% from 12.68x to 0.00x"),
        "caveats": [
            "docked a held-out azole *class*, NOT fluconazole specifically -- the "
            "fluconazole-level claim is an inference from the class result, not a measurement",
            "the top-rank collapse is a sub-0.1 A reshuffle among near-tied molecules, right "
            "at the tool's per-pose noise floor",
            "a rigid single-atom OH deletion under-models clinical Y132F (no conformational "
            "relaxation, no H-bond-network rearrangement, no co-occurring efflux/expression)",
        ],
    },
}


def lost_pocket_contacts(mv, residues, ligand_atoms, cutoff=mapping.CONTACT_CUTOFF):
    """Removed atoms of a deletion mutation that were within `cutoff` of the co-crystal drug.

    `mv` -- an openafr.mapping.map_mutation verdict. Returns a list of
    {"atom", "dist_to_drug"} for each deleted heavy atom whose closest approach to the
    reference ligand is <= cutoff, sorted nearest-first. Empty for a non-deletion mutation,
    an unmapped verdict, or a deletion whose removed atoms are all second-shell. This is the
    deterministic core of the geometric estimate: the specific favorable contacts a deletion
    takes away from the drug.
    """
    if not mv.get("mappable") or not mv.get("buildable"):
        return []
    removed = DELETION_ATOMS.get((mv["wt"], mv["mut"]), set())
    _, atoms = residues[mv["position"]]
    out = []
    for name, xyz in atoms:
        if name in removed:
            d = mapping._min_dist([(name, xyz)], ligand_atoms)
            if d is not None and d <= cutoff:
                out.append({"atom": name, "dist_to_drug": d})
    return sorted(out, key=lambda c: c["dist_to_drug"])


def estimate_fit_consequence(mv, residues, ligand_atoms, cutoff=mapping.CONTACT_CUTOFF):
    """Estimate the azole-fit consequence of one mapping verdict `mv`.

    Returns a self-describing dict carrying the mapping facts plus:
      evidence_class -- 'empirical-dock' | 'geometric-estimate' | 'no-call'
      confidence     -- 'measured' | 'estimate' | 'none'
      direction      -- 'weaker-fit' | 'minimal-direct-effect' | 'uncertain'
      fluconazole_fit_verdict -- one plain-English sentence a non-bioinformatician can act on
      empirical      -- the EMPIRICAL_DOCKS entry when one exists, else None
      lost_contacts  -- geometric lost-contact list (see lost_pocket_contacts)
      caveats        -- the honesty ledger for this call

    The branching mirrors the three evidence tiers in the module docstring; every path
    ends in an explicit confidence, and 'no-call' is an expected, first-class outcome.
    """
    out = {
        "token": mv.get("token"),
        "mappable": mv.get("mappable", False),
        "pocket_contact": mv.get("pocket_contact", False),
        "buildable": mv.get("buildable", False),
        "empirical": None,
        "lost_contacts": [],
        "caveats": [],
    }

    # Tier 3a: nothing to place in the pocket -> no structural fit call at all.
    if not mv.get("mappable"):
        reason = mv.get("reason", "token does not map to a pocket residue")
        return {**out, "evidence_class": "no-call", "confidence": "none",
                "direction": "uncertain",
                "fluconazole_fit_verdict":
                    f"no structural fit call -- {reason}",
                "caveats": ["not a mappable point substitution; there is no residue to "
                            "reshape the pocket, so azole fit cannot be reasoned about here"]}

    label = f"{mv['wt']}{mv['position']}{mv['mut']}"
    lost = lost_pocket_contacts(mv, residues, ligand_atoms, cutoff)
    out["lost_contacts"] = lost

    # Tier 1: a real, pre-registered dock exists -> report the measured outcome.
    emp = EMPIRICAL_DOCKS.get(out["token"])
    if emp is not None:
        out["empirical"] = emp
        out["caveats"] = list(emp["caveats"])
        verdict = (
            f"{label}: DOCKED (#13, gate {emp['gate']}). {emp['finding']}. "
            f"Read: fluconazole is not structurally *excluded* -- azoles still reach the "
            f"iron -- but the best-fitting azoles lose their top-rank edge. Measured, not "
            f"estimated; see {emp['source']}.")
        return {**out, "evidence_class": "empirical-dock", "confidence": "measured",
                "direction": emp["direction"], "fluconazole_fit_verdict": verdict}

    # Tier 2: deletion-buildable, no dock yet -> deterministic *direction*, never a number.
    if mv.get("buildable"):
        cmd = (f"python scripts/mutate_receptor.py {out['token']} "
               f"-o work/receptor_auris_{out['token']}.pdb  (then prep_receptor + screen.sh, "
               f"reusing the #13 pipeline) to measure the magnitude")
        if lost:
            near = lost[0]
            names = ", ".join(c["atom"] for c in lost)
            verdict = (
                f"{label}: ESTIMATE (buildable, not yet docked). Deletion removes "
                f"{len(lost)} pocket-wall atom(s) [{names}] that contact the co-crystal "
                f"azole (closest {near['dist_to_drug']:.2f} A) -> predicted WEAKER azole "
                f"fit; magnitude unmeasured. {cmd}.")
            return {**out, "evidence_class": "geometric-estimate", "confidence": "estimate",
                    "direction": "weaker-fit", "fluconazole_fit_verdict": verdict,
                    "caveats": ["direction only -- a rigid atom deletion predicts the sign of "
                                "the effect, not its size; dock to quantify",
                                "a lost static contact need not abolish binding; azoles are "
                                "multi-contact and can re-pose"]}
        verdict = (
            f"{label}: ESTIMATE (buildable, not yet docked). The deletion removes only "
            f"second-shell atoms (none within {cutoff:.1f} A of the drug) -> minimal "
            f"predicted DIRECT fit change; any real effect would be indirect (pocket "
            f"dynamics) and a static model cannot call it. {cmd} if you want the number.")
        return {**out, "evidence_class": "geometric-estimate", "confidence": "estimate",
                "direction": "minimal-direct-effect", "fluconazole_fit_verdict": verdict,
                "caveats": ["static, direct-contact reasoning only; does not see allosteric "
                            "or dynamic effects a second-shell change might have"]}

    # Tier 3b: mapped but the pinned env cannot build it (add-atom substitution).
    if mv.get("pocket_contact"):
        verdict = (
            f"{label}: NO CONFIDENT CALL (plausible). Sits directly on the drug-contact "
            f"surface ({mv['dist_to_ligand']:.2f} A), so an azole-fit effect is plausible -- "
            f"but this substitution ADDS side-chain atoms and the pinned env has no repacker, "
            f"so the mutant cannot be built or docked. Resolve with repacker tooling.")
        return {**out, "evidence_class": "no-call", "confidence": "none",
                "direction": "uncertain", "fluconazole_fit_verdict": verdict,
                "caveats": ["unbuildable in the pinned env (needs a side-chain repacker); "
                            "pocket-contact makes a fit effect plausible but unresolved -- "
                            "we do not guess an added side chain"]}
    verdict = (
        f"{label}: NO CONFIDENT CALL (low prior). Second-shell ({mv['dist_to_ligand']:.2f} A, "
        f"no drug contact) AND unbuildable (adds atoms, needs a repacker). Low prior for a "
        f"direct fit effect, and no way to dock it in the pinned env.")
    return {**out, "evidence_class": "no-call", "confidence": "none",
            "direction": "uncertain", "fluconazole_fit_verdict": verdict,
            "caveats": ["second-shell and unbuildable; both a low prior for a direct effect "
                        "and no dockable path -- reported as unresolved, not as 'no effect'"]}


def estimate_signals(signals, receptor=mapping.DEFAULT_RECEPTOR,
                     ligand=mapping.DEFAULT_LIGAND, cutoff=mapping.CONTACT_CUTOFF):
    """Map + estimate a list of emergence signals (or bare tokens) in one call.

    Each element is an emergence signal dict (uses its 'mutation' key) or a bare token
    string. Returns fit-consequence verdicts in input order, closing the whole
    detect -> map -> estimate loop. Loads the structure once.
    """
    from openafr.pdbqt import iron_position
    residues = mapping.read_residues(receptor)
    ligand_atoms = mapping.read_ligand_atoms(ligand)
    fe = iron_position(receptor)
    tokens = [s["mutation"] if isinstance(s, dict) else s for s in signals]
    verdicts = [mapping.map_mutation(t, residues, ligand_atoms, fe, cutoff) for t in tokens]
    return [estimate_fit_consequence(mv, residues, ligand_atoms, cutoff) for mv in verdicts]
