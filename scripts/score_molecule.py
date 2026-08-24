"""Score your own molecule: SMILES -> geometry verdict, the whole validated pipeline in one command.

Pipeline stage: the outward-facing front door (end-to-end). This is what lets a
collaborator try *their* compound against the validated CYP51 method without knowing the
repo's internals. It chains the four stages that already exist as separate scripts —
applicability triage, ligand prep, docking, and the N-to-Fe geometry read — behind one
honest command, and never fakes the parts that need a cloud VM.

    SMILES
      |  (1) applicability triage        rdkit only, runs anywhere
      v      does the validated method even claim to cover this molecule, and is it novel?
    in-envelope? --no--> STOP with the reason (out-of-scope / out-of-domain / unparseable).
      |                  The method makes no claim here; docking it would imply one.
     yes
      |  (2) prep_ligands.py             SMILES -> deterministic 3D SDF
      |  (3) screen.sh                   vina under the hash-frozen protocol   [needs openafr env]
      |  (4) analyze_poses.analyze       closest N-Fe, iron-bound pose count
      v
    per-molecule geometry readout, contextualized against the validated actives' iron-bound
    band, carrying every standing caveat (see the header of the printed report).

Two honest gates, both load-bearing:
  * Only molecules the triage puts *inside* the validated envelope proceed to docking. A
    molecule the method has no validated basis for is reported and stopped, never docked.
  * Docking needs vina + obabel and a prepared receptor, provided by the pinned `openafr`
    conda env (these install on Apple Silicon and Linux alike — unlike the recaller's
    Linux/x86-only reads tools, docking has no cloud-VM requirement). If they are absent on
    this machine this prints the *exact* command to finish the run and stops cleanly — it
    never silently skips the geometry or invents a number.

What a geometry readout does and does NOT mean (carried into the printed report):
  * A high rank / tight N-Fe is a *hypothesis for a wet lab*, never a hit.
  * It is meaningful only if `validate_gate2.py` has PASSED under the frozen protocol.
  * Reaching the iron like the validated actives is necessary, NOT sufficient: property-
    matched decoys also reach it (work/RESULTS_axial.md), which is the decoy ceiling the
    preprint records. This tool triages; it does not predict activity.
  * The durable validated metric is AUC (~0.79); the top-1% EF was a single-run artifact.

Usage:
    conda run -n openafr python scripts/score_molecule.py --smiles "<SMILES>" [--name X]
    conda run -n openafr python scripts/score_molecule.py LIBRARY.smi   # SMILES<TAB>name
    # options: --receptor-pdbqt work/receptor.pdbqt  --receptor-pdb work/receptor_A.pdb
    #          --workdir work/score  --triage-only  --json
    # --json emits one machine-readable object on stdout (human text -> stderr).
"""
import argparse
import json
import os
import pathlib
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from scripts.analyze_poses import FE_CUTOFF, analyze
from scripts.applicability_triage import VERDICTS, classify, load_panel, read_library

# The validated actives' iron-bound band: the min/max closest N-Fe over every azole that
# produced an iron-bound pose (work/RESULTS_all_azoles.md: 2.47 A clotrimazole .. 2.88 A
# itraconazole). A candidate that lands in this band coordinates the heme iron the way the
# real drugs do. This is a reference for interpretation only — not a pass bar.
VALIDATED_FE_BAND = (2.47, 2.88)

# Triage verdicts that earn a docking run: the molecule is inside the validated envelope.
# Everything else (out-of-scope / out-of-domain / unparseable) is reported and stopped —
# the method has no validated basis for it, so it must not be docked.
DOCKABLE = {"in-envelope-novel", "near-known"}

REQUIRED_TOOLS = ("obabel", "vina")


def plan(rows, panel):
    """Triage every molecule and decide which proceed to docking. Pure; no compute spent.

    rows: [(name, smiles)]. Returns [(name, smiles, triage_dict, dock_bool)], preserving
    input order, so nothing is silently dropped — a molecule the method cannot cover still
    gets a row, it just carries dock=False.
    """
    out = []
    for name, smi in rows:
        t = classify(smi, panel)
        out.append((name, smi, t, t["verdict"] in DOCKABLE))
    return out


def classify_geometry(best_fe, n_passing):
    """Interpret one molecule's geometry against the validated actives' band.

    best_fe: closest nitrogen-to-iron approach over all poses (A), or None if no pose had a
    nitrogen. n_passing: how many poses reached coordination (N-Fe < FE_CUTOFF).
    Returns (label, coordinating_bool). Pure — the interpretation the report renders.
    """
    if best_fe is None:
        return ("no pose reached the pocket with a nitrogen — nothing to interpret", False)
    lo, hi = VALIDATED_FE_BAND
    if n_passing:
        if best_fe < lo:
            return (f"iron-bound, tighter than the validated actives' band "
                    f"({lo:.2f}-{hi:.2f} A)", True)
        if best_fe <= hi:
            return (f"iron-bound, within the validated actives' band "
                    f"({lo:.2f}-{hi:.2f} A) — coordinates like the real drugs", True)
        return (f"iron-bound but at the loose edge (beyond the {hi:.2f} A actives band, "
                f"within the {FE_CUTOFF:.1f} A cutoff)", True)
    return (f"NOT iron-bound: closest N is {best_fe:.2f} A, beyond the {FE_CUTOFF:.1f} A "
            f"coordination cutoff", False)


def render_molecule(name, triage, geom, pending):
    """One molecule's readout as text. triage: classify() dict. pending=True means the
    molecule is inside the envelope but was not docked (triage-only, or the VM step is still
    to run). geom: analyze() dict, or None when docking ran but produced no usable pose. Pure."""
    L = [f"## {name}"]
    heavy = "-" if triage["heavy"] is None else triage["heavy"]
    tc = "-" if triage["max_tc"] is None else f"{triage['max_tc']:.3f}"
    L.append(f"- triage   : {triage['verdict']} — {VERDICTS[triage['verdict']]}")
    L.append(f"             heavy atoms {heavy}, nearest validated active "
             f"{triage['nearest'] or '-'} (Tanimoto {tc})")
    if triage["verdict"] not in DOCKABLE:
        L.append("- geometry : not docked — the validated method makes no claim about this "
                 "molecule, so no pose was computed.")
        return "\n".join(L)
    if pending:
        L.append("- geometry : PENDING — inside the envelope; run the dock to get the "
                 "N-Fe read (scripts/score_molecule.py without --triage-only in the openafr "
                 "env, with vina + a prepared receptor).")
        return "\n".join(L)
    if geom is None:
        L.append("- geometry : docking produced NO usable pose — no ranking evidence for "
                 "this molecule.")
        return "\n".join(L)
    label, _ = classify_geometry(geom["best_fe"], geom["n_passing"])
    best_fe = "no nitrogen pose" if geom["best_fe"] is None else f"{geom['best_fe']:.3f} A"
    sc = "NA" if geom["best_score"] is None else f"{geom['best_score']:.1f} kcal/mol"
    L.append(f"- geometry : {label}")
    L.append(f"             closest N-Fe {best_fe} | iron-bound poses "
             f"{geom['n_passing']}/{geom['n_poses']} | Vina {sc} (secondary, does not decide)")
    return "\n".join(L)


def molecule_record(name, smiles, triage, geom, pending):
    """One molecule's readout as a JSON-serializable dict (the machine mirror of
    render_molecule). Pure. `geometry_status` is the single stable cause a caller keys on:
    not-docked (outside the envelope) / pending (in-envelope, dock still to run) / no-pose
    (docked, no usable pose) / scored. `geometry` is populated only when scored."""
    rec = {
        "name": name,
        "smiles": smiles,
        "triage": triage["verdict"],
        "heavy_atoms": triage["heavy"],
        "nearest_active": triage["nearest"],
        "max_tanimoto": triage["max_tc"],
        "dockable": triage["verdict"] in DOCKABLE,
        "geometry": None,
    }
    if triage["verdict"] not in DOCKABLE:
        rec["geometry_status"] = "not-docked"
    elif pending:
        rec["geometry_status"] = "pending"
    elif geom is None:
        rec["geometry_status"] = "no-pose"
    else:
        label, coordinating = classify_geometry(geom["best_fe"], geom["n_passing"])
        rec["geometry_status"] = "scored"
        rec["geometry"] = {
            "closest_n_fe": geom["best_fe"],
            "iron_bound_poses": geom["n_passing"],
            "n_poses": geom["n_poses"],
            "vina_score": geom["best_score"],
            "coordinates_iron": coordinating,
            "interpretation": label,
        }
    return rec


HEADER_CAVEATS = (
    "OpenAFR — score your own molecule\n"
    "=================================\n"
    "This ranks a molecule by ONE validated criterion: how closely a nitrogen reaches the\n"
    "heme iron (smaller = better). Read every line below as a hypothesis for a wet lab,\n"
    "never a hit. It is meaningful only if validate_gate2.py has PASSED under the frozen\n"
    "protocol. Reaching the iron like the validated actives is necessary but NOT sufficient:\n"
    "property-matched decoys reach it too (the decoy ceiling, work/RESULTS_axial.md). The\n"
    "durable validated metric is AUC (~0.79); the top-1% enrichment was a single-run artifact.\n"
)


def _dock_instructions(smi_path, workdir, receptor_pdbqt, receptor_pdb):
    """The exact commands to finish an in-envelope run on the VM, when tools are absent."""
    lig_dir = os.path.join(workdir, "ligands")
    screen_dir = os.path.join(workdir, "screen")
    return (
        "Docking needs vina + obabel + a prepared receptor from the pinned openafr env.\n"
        "They are not available on this machine, so the geometry read is PENDING. To finish:\n\n"
        f"  # (once) prepare the receptor if {receptor_pdbqt} is missing:\n"
        f"  conda run -n openafr python scripts/prep_receptor.py\n\n"
        f"  conda run -n openafr python scripts/prep_ligands.py {smi_path} {lig_dir}\n"
        f"  RECEPTOR={receptor_pdbqt} scripts/screen.sh {lig_dir} {screen_dir}\n"
        f"  conda run -n openafr python scripts/score_molecule.py {smi_path} \\\n"
        f"      --from-screen {screen_dir} --receptor-pdb {receptor_pdb}\n"
    )


def _run_docking(rows, workdir, receptor_pdbqt, receptor_pdb):
    """Prep -> dock the given [(name, smiles)] rows and return {name: analyze() dict}.

    The only part that shells out to external tools (prep_ligands, screen.sh). Kept thin so
    the decision + interpretation logic above stays pure and offline-testable.
    """
    os.makedirs(workdir, exist_ok=True)
    smi_path = os.path.join(workdir, "candidates.smi")
    with open(smi_path, "w") as fh:
        for name, smi in rows:
            fh.write(f"{smi}\t{name}\n")
    lig_dir = os.path.join(workdir, "ligands")
    screen_dir = os.path.join(workdir, "screen")
    subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "prep_ligands.py"),
                    smi_path, lig_dir], check=True)
    subprocess.run([os.path.join(ROOT, "scripts", "screen.sh"), lig_dir, screen_dir],
                   check=True, env={**os.environ, "RECEPTOR": receptor_pdbqt})
    return _collect_geometry(screen_dir, receptor_pdb)


def _collect_geometry(screen_dir, receptor_pdb):
    """{name: analyze() dict} for every docked pose file in screen_dir."""
    from openafr.pdbqt import iron_position
    fe = iron_position(receptor_pdb)
    geom = {}
    for f in sorted(os.listdir(screen_dir)):
        if f.endswith(".pdbqt"):
            geom[f[:-6]] = analyze(os.path.join(screen_dir, f), fe)
    return geom


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("library", nargs="?", help="SMILES<TAB>name library file")
    ap.add_argument("--smiles", help="a single SMILES to score")
    ap.add_argument("--name", default="query", help="name for --smiles (default: query)")
    ap.add_argument("--receptor-pdbqt", default="work/receptor.pdbqt",
                    help="receptor for docking (default work/receptor.pdbqt)")
    ap.add_argument("--receptor-pdb", default="work/receptor_A.pdb",
                    help="receptor with the heme iron for geometry (default work/receptor_A.pdb)")
    ap.add_argument("--workdir", default="work/score",
                    help="where prepped ligands + poses land (default work/score)")
    ap.add_argument("--triage-only", action="store_true",
                    help="stop after triage; never dock (runs anywhere)")
    ap.add_argument("--from-screen", metavar="SCREEN_DIR",
                    help="skip docking; read poses already docked into SCREEN_DIR")
    ap.add_argument("--json", action="store_true",
                    help="emit one JSON object on stdout (caveats + per-molecule records); "
                         "human/progress text goes to stderr")
    args = ap.parse_args()

    if bool(args.library) == bool(args.smiles):
        ap.error("give exactly one of: a LIBRARY file, or --smiles")

    # In --json mode, stdout must be pure JSON; route all human/progress text to stderr.
    info = (lambda *a, **k: print(*a, file=sys.stderr, **k)) if args.json else print

    rows = [(args.name, args.smiles)] if args.smiles else read_library(args.library)
    panel = load_panel()
    planned = plan(rows, panel)

    info(HEADER_CAVEATS)
    info(f"validated panel: {len(panel)} actives | domain <= 45 heavy | "
         f"actives iron-bound band {VALIDATED_FE_BAND[0]:.2f}-{VALIDATED_FE_BAND[1]:.2f} A\n")

    dock_rows = [(n, s) for n, s, _, d in planned if d]
    geom = {}
    docked = False   # did a docking step actually run (vs. pending on the VM)?

    if dock_rows and not args.triage_only:
        if args.from_screen:
            geom = _collect_geometry(args.from_screen, args.receptor_pdb)
            docked = True
        else:
            missing = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
            receptor_ready = pathlib.Path(args.receptor_pdbqt).is_file()
            if missing or not receptor_ready:
                smi_path = os.path.join(args.workdir, "candidates.smi")
                os.makedirs(args.workdir, exist_ok=True)
                with open(smi_path, "w") as fh:
                    for n, s in dock_rows:
                        fh.write(f"{s}\t{n}\n")
                why = (f"missing tools: {', '.join(missing)}" if missing
                       else f"receptor not prepared: {args.receptor_pdbqt}")
                info(f"[{len(dock_rows)} molecule(s) inside the envelope — geometry PENDING] "
                     f"{why}\n")
                info(_dock_instructions(smi_path, args.workdir,
                                        args.receptor_pdbqt, args.receptor_pdb))
                info()
            else:
                geom = _run_docking(dock_rows, args.workdir,
                                    args.receptor_pdbqt, args.receptor_pdb)
                docked = True

    records = []
    for name, smi, triage, dock in planned:
        # pending = in-envelope but no docking ran yet (triage-only, or VM step still to do).
        # When docking DID run, a name absent from geom means it produced no pose (None).
        pending = dock and not docked
        g = geom.get(name)
        if args.json:
            records.append(molecule_record(name, smi, triage, g, pending))
        else:
            print(render_molecule(name, triage, g, pending))
            print()

    if args.json:
        print(json.dumps({
            "tool": "openafr/score_molecule",
            "caveats": HEADER_CAVEATS.strip(),
            "validated_band_angstrom": list(VALIDATED_FE_BAND),
            "molecules": records,
        }, indent=2))


if __name__ == "__main__":
    main()
