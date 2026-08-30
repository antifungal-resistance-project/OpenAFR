"""Pure decision logic for the consolidated per-candidate confidence scorecard (#88).

The scorecard rolls every *already-measured* confidence axis onto one row per
candidate and marks every *unmeasured* axis `not-yet-checked` — so a card can
never look complete by omission, and no axis is ever blended into a single score.
This module holds the two rules a test can pin:

  - `derive_concerns` maps one candidate's per-axis records to a flat list of
    human-readable concerns (one string per triggered axis);
  - `NOT_YET_CHECKED` is the fixed set of axes no candidate has been scored on yet
    (open issues), which the builder stamps onto every card verbatim.

The runnable merge that reads the frozen per-axis files lives in
`work/candidates/build_scorecard.py`; it imports from here.
"""

FE_CUTOFF_A = 3.0  # N-Fe within this = coordinates the heme iron (frozen protocol)

# Known azole controls — present to validate the axes, not discoveries to test.
CONTROLS = {
    "flutrimazole", "ravuconazole", "voriconazole", "letrozole",
    "fluconazole", "ketoconazole", "fosfluconazole", "luliconazole",
    "croconazole", "clotrimazole",
}

# Axes NOT yet measured per candidate (open issues) — named so a card can never
# look complete by omission.
NOT_YET_CHECKED = {
    # ensemble/flexible docking (#70) is now measured (ensemble_consistency axis,
    # RESULTS_candidate_ensemble.md): each candidate is scored across the exhaustive
    # crystallographic C. albicans ensemble {5TZ1, 5FSA, 5V5Z} and its coordination is
    # ROBUST / CONSISTENT / CONFORMER-FRAGILE / NON-COORDINATING. This closes the
    # "single-snapshot artifact" hole for the RIGID crystal ensemble, not the induced-fit /
    # MD-flexibility one (no MD stack in the pinned env — see the prereg limitations).
    # pose-convergence (#71) is now measured (pose_convergence axis,
    # RESULTS_pose_convergence.md): within-search RMSD clustering of the frozen seed-42 poses.
    # Y132F retention is now measured (resistance_retention axis, RESULTS_candidate_resistance.md);
    # the OTHER C. auris mutants still need a side-chain repacker not in the pinned env.
    "auris_other_mutants": "auris K143R/F126L mutant retention (#77; Y132F now measured)",
    # protonation/tautomer (#84) is now measured (protomer axis, RESULTS_candidate_protomer.md):
    # each candidate's pH-7.4 microstates (obabel-band protonation x RDKit tautomers) are docked;
    # 6/16 are single-microstate (EXPLICIT), the rest ROBUST/MICROSTATE-DEPENDENT. This closes the
    # "rank silently tied to one protomer/tautomer" hole, not the "obabel pKa is wrong" one.
    # stereochemistry (#85) is now measured (stereochemistry axis, RESULTS_candidate_stereo.md):
    # every novel candidate's supplied SMILES pins a single defined isomer (EXPLICIT); the two
    # ambiguous controls are docked per-isomer. This closes the "unspecified stereo" hole, not
    # the "supplied the wrong configuration" one (see the prereg limitations).
    # broad human-CYP off-target panel (#74) is now measured (cyp_offtarget axis,
    # RESULTS_cyp_offtarget.md): each candidate carries a per-isoform inhibition-liability
    # profile for CYP3A4/2C9/2D6/2C19/1A2 from local mechanistic (type-II heme ligation) +
    # substrate-pharmacophore reasoning. This closes the "silent on the human-CYP DDI liability
    # the azole warhead carries" hole, not the "measured IC50" one (it is a flag, not a docking
    # or QSAR result — see RESULTS_cyp_offtarget.md).
}


def derive_concerns(cf, am, sy, dm, pr, is_ctrl, res=None, st=None, pt=None, en=None,
                    cyp=None):
    """Given a candidate's per-axis records, list its open concerns.

    Each rule maps one measured axis to a human-readable concern string, so the
    scorecard never blends axes into a single number. A control is exempt only from
    the "known CYP51 inhibitor", "phenotypic active elsewhere", "loses coordination
    in the resistant pocket", "conformer-fragile" and "human-CYP type-II" demotions, which
    are the point of a control, not a strike against it (controls are the Y132F-defeated
    azoles; how they behave on the mutant / across conformers, and that they carry the
    azole warhead that inhibits human CYPs, is a robustness check, not a mark against a
    discovery — ketoconazole IS the textbook CYP3A4 inhibitor).

    Args mirror the per-axis file rows, keyed by candidate name:
      cf   shortlist_confidence.json    (seed stability + human selectivity)
      am   shortlist_admet.tsv          (ADMET flags)
      sy   shortlist_synth.json         (synthesizability band)
      dm   domain_sensitivity.json      (candidate_fragility entry)
      pr   top100_precedent.tsv         (assay-precedent verdicts)
      res  shortlist_resistance.json    (C. auris Y132F retention verdict; None = not-yet-checked)
      st   shortlist_stereo.json        (stereochemistry verdict; None = not-yet-checked)
      pt   shortlist_protomer.json      (protonation/tautomer verdict; None = not-yet-checked)
      en   shortlist_ensemble.json      (ensemble-consistency verdict; None = not-yet-checked)
      cyp  shortlist_cyp.json           (human-CYP off-target liability; None = not-yet-checked)
    """
    concerns = []
    if cf.get("reaches_fungal_iron") is False:
        concerns.append("does not reach fungal iron on consensus")
    if cf.get("stable") is False:
        concerns.append("seed-unstable (single-dock artifact)")
    margin = cf.get("selectivity_margin")
    if margin is not None and margin <= 0:
        concerns.append("no fungal-vs-human selectivity window")
    admet_flags = am.get("flags", "")
    if admet_flags and admet_flags != "clean":
        concerns.append("ADMET: " + admet_flags)
    if sy.get("sa_band") == "hard":
        concerns.append("hard to synthesize (SAscore>6)")
    for fl in dm.get("flags", []):
        concerns.append("domain: " + fl)
    if pr.get("phenotypic") == "tested-active" and not is_ctrl:
        concerns.append("phenotypic tested-active elsewhere (read assay)")
    if pr.get("verdict") == "tested-active":
        concerns.append("known CYP51 inhibitor (control)")
    # Resistance retention (#69,#77): a molecule whose coordination collapses in the
    # C. auris Y132F pocket is demoted — the geometry that shortlisted it is fragile in
    # the target the mission actually cares about. Control-exempt (see docstring).
    if res is not None and res.get("verdict") == "LOST" and not is_ctrl:
        concerns.append("loses coordination in C. auris Y132F pocket")
    # Stereochemistry (#85): a candidate whose rank flips with an UNSPECIFIED stereocenter is
    # demoted — the shortlist SMILES did not pin the isomer the geometry rests on. Control-exempt:
    # the ambiguous molecules are the racemic azole controls, docked per-isomer to exercise the
    # axis, not discoveries to strike. An EXPLICIT candidate never triggers this.
    if st is not None and st.get("verdict") == "ISOMER-DEPENDENT" and not is_ctrl:
        concerns.append("rank depends on unspecified stereochemistry")
    # Protonation/tautomer (#84): a candidate whose coordination collapses or shifts across the
    # protonation/tautomer microstates populated at pH 7.4 is demoted — the shortlist docked one
    # microstate and the rank does not survive the others. Control-exempt: the multi-state
    # molecules include the racemic/tautomeric azole controls, docked per-microstate to exercise
    # the axis, not discoveries to strike. A protomer-EXPLICIT candidate never triggers this.
    if pt is not None and pt.get("verdict") == "MICROSTATE-DEPENDENT" and not is_ctrl:
        concerns.append("rank depends on protonation/tautomer state at pH 7.4")
    # Ensemble consistency (#70): a candidate whose coordination appears in only one crystal
    # conformer (CONFORMER-FRAGILE) or in none (NON-COORDINATING) across the exhaustive
    # {5TZ1,5FSA,5V5Z} ensemble is demoted — the geometry that shortlisted it is a
    # single-snapshot artifact, not a robust binding hypothesis. Control-exempt (see docstring):
    # a control's conformer behavior validates the axis, it is not a strike. A ROBUST or
    # CONSISTENT candidate never triggers this.
    if en is not None and en.get("verdict") in ("CONFORMER-FRAGILE", "NON-COORDINATING") \
            and not is_ctrl:
        concerns.append("coordination is conformer-fragile (single-snapshot artifact)")
    # Human-CYP off-target (#74): a candidate carrying an azole warhead (imidazole/triazole/
    # tetrazole) is flagged HIGH type-II across the human CYP panel — it inhibits the human
    # drug-metabolizing CYPs by the SAME heme-ligation the geometry criterion selected for,
    # inheriting the azole class's CYP3A4 drug-interaction / hepatotoxicity liability. This is
    # the concern the whole panel exists to surface. Control-exempt (see docstring): the
    # controls ARE the azole CYP3A4 inhibitors — ketoconazole firing HIGH validates the axis,
    # it is not a strike. Weaker (MODERATE) azine/pharmacophore hits are reported on the card
    # but do not demote (nearly every heme-seeking scaffold carries one — too broad to strike).
    if cyp is not None and cyp.get("n_high", 0) > 0 and not is_ctrl:
        concerns.append("human-CYP type-II inhibition liability (azole warhead; CYP3A4 DDI risk)")
    return concerns
