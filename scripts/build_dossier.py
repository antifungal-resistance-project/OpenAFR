#!/usr/bin/env python
"""Turn a screened, ranked library into one machine-readable candidate dossier per molecule.

Pipeline stage: report/  (runs after screen.sh, alongside report_candidates.py)

Why this exists
---------------
report_candidates.py fuses geometry + precedent into human-readable Markdown. This CLI
produces the machine-readable sibling: for each molecule it assembles ALL of the evidence
the pipeline already computes — geometry rank (reused from report_candidates.collect),
drug-likeness (admet.profile), human-CYP off-target liability (cyp_offtarget.profile),
synthesizability (synth.profile), and — over the network — assay precedent
(precedent.lookup) — into one JSON dossier carrying an A/B/C experimental-priority verdict
and an explicit, sourced reason for selection. The fusion + priority logic lives in the
pure, unit-tested openafr.dossier; this file only does the I/O and the profiler calls.

It computes nothing about geometry and does NOT re-run the gate: a dossier is a hypothesis
for a wet lab, never a hit, and means nothing unless validate_gate2.py has PASSED under the
frozen protocol first. It never drops a molecule — a bad profile lowers priority, the row
stays. It fails per-molecule, never per-run: an unparseable SMILES or a precedent lookup
error is recorded on that molecule (priority C) and the batch continues.

The local profilers (admet/cyp/synth) always run and are deterministic. Precedent is the
only network stage; --precedent none produces a fully offline, reproducible dossier.

Usage:
  python scripts/build_dossier.py SCREEN_DIR LIBRARY.smi [RECEPTOR_PDB] \\
      [-o dossiers.json] [--split DIR] [-n TOP] \\
      [--precedent both|chembl|pubchem|none] [--near TESTED.smi]

  SCREEN_DIR    directory of docked *.pdbqt poses (e.g. work/screen)   — geometry, reused
  LIBRARY.smi   the SAME library, `<SMILES>\\t<name>` per line (admet_filter.read_smi format)
  RECEPTOR_PDB  receptor with the heme iron (default work/receptor_A.pdb)
  -o / --out    write the JSON array here (default: stdout)
  --split DIR   also write one <name>.json per molecule (candidate registry)
  -n / --top    cap dossiers to the top-N ranked molecules (default: all)
  --precedent   network precedent source; 'none' skips the network (default: both)
  --near        nearest previously-tested neighbour annotation (soft negatives)
"""
import argparse
import datetime
import hashlib
import json
import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import admet, cyp_offtarget, dossier, precedent, protocol, synth
from scripts.admet_filter import read_smi
from scripts.report_candidates import collect


def _ranked_geometry(screen_dir, receptor_pdb):
    """Reuse report_candidates.collect and attach a 1-based rank (rankable rows only)."""
    fe, rows = collect(screen_dir, receptor_pdb)
    rank = 0
    for r in rows:
        if r["best_fe"] is not None:
            rank += 1
            r["rank"] = rank
        else:
            r["rank"] = None
    return rows


def _profiles(smiles, name, catalogued):
    """Local, deterministic profilers. Returns (admet, cyp, synth, error)."""
    try:
        a = admet.profile(smiles)
        c = cyp_offtarget.profile(smiles)
        c["flags"] = cyp_offtarget.flag_summary(c)
        s = synth.profile(smiles, name=name, catalogued=catalogued)
        return a, c, s, None
    except ValueError as e:               # unparseable SMILES — record, don't crash the batch
        return None, None, None, str(e)


def _git_version():
    try:
        return subprocess.check_output(
            ["git", "describe", "--tags", "--always", "--dirty"],
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unknown"


def _provenance(receptor_pdb, source):
    digest = hashlib.sha256(protocol.PROTOCOL.read_bytes()).hexdigest()
    return {
        "protocol_sha256": digest,
        "protocol_unmodified": digest == protocol.PROTOCOL_SHA,
        "receptor": receptor_pdb,
        "openafr_version": _git_version(),
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "precedent_source": source,
    }


def build(screen_dir, library_smi, receptor_pdb, top, source, near_path):
    rows = _ranked_geometry(screen_dir, receptor_pdb)
    if top is not None:
        rows = rows[:top]
    name2smi = {name: smi for name, smi in read_smi(library_smi)}
    neighbours = precedent.load_neighbours(near_path) if near_path else None

    records = []
    for r in rows:
        name = r["name"]
        smi = name2smi.get(name)
        geometry = {"rank": r["rank"], "best_fe": r["best_fe"], "n_passing": r["n_passing"],
                    "n_poses": r["n_poses"], "best_score": r["best_score"]}
        if smi is None:                    # in the screen but not the library — no structure
            records.append({"name": name, "smiles": None, "geometry": geometry,
                            "error": "not-in-library"})
            continue

        a, c, s, err = _profiles(smi, name, catalogued=False)
        rec = {"name": name, "smiles": smi, "geometry": geometry,
               "admet": a, "cyp": c, "synth": s, "error": err}

        if err is None and source != "none":
            try:
                rec["prec"] = precedent.lookup(smi, sources=tuple(source.split("+")))
            except Exception as e:         # network hiccup — isolate per molecule
                rec["prec"] = {"verdict": "lookup-failed", "error": str(e)}
        if err is None and neighbours:
            fp = precedent.fingerprint(smi)
            best = precedent.nearest_tested_neighbour(fp, neighbours) if fp is not None else None
            if best:
                rec["near"] = {"tanimoto": round(best["similarity"], 3),
                               "name": best["id"], "verdict": best["verdict"]}
        records.append(rec)

    return dossier.build_all(records, _provenance(receptor_pdb, source))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("screen_dir", help="directory of docked *.pdbqt poses")
    ap.add_argument("library_smi", help="the screened library, `<SMILES>\\t<name>` per line")
    ap.add_argument("receptor_pdb", nargs="?", default="work/receptor_A.pdb",
                    help="receptor with the heme iron (default work/receptor_A.pdb)")
    ap.add_argument("-o", "--out", help="write the JSON array here (default: stdout)")
    ap.add_argument("--split", metavar="DIR", help="also write one <name>.json per molecule")
    ap.add_argument("-n", "--top", type=int, default=None,
                    help="cap to the top-N ranked molecules (default: all)")
    ap.add_argument("--precedent", default="both",
                    choices=["both", "chembl", "pubchem", "none"],
                    help="network precedent source; 'none' skips the network (default: both)")
    ap.add_argument("--near", metavar="TESTED.smi",
                    help="nearest previously-tested neighbour annotation (soft negatives)")
    args = ap.parse_args(argv)

    source = "chembl+pubchem" if args.precedent == "both" else args.precedent
    dossiers = build(args.screen_dir, args.library_smi, args.receptor_pdb,
                     args.top, source, args.near)

    text = json.dumps(dossiers, indent=2, sort_keys=True)
    if args.out:
        pathlib.Path(args.out).write_text(text + "\n")
        sys.stderr.write("wrote %d dossiers to %s\n" % (len(dossiers), args.out))
    else:
        print(text)

    if args.split:
        d = pathlib.Path(args.split)
        d.mkdir(parents=True, exist_ok=True)
        for one in dossiers:
            (d / ("%s.json" % one["name"])).write_text(
                json.dumps(one, indent=2, sort_keys=True) + "\n")
        sys.stderr.write("wrote %d per-molecule dossiers to %s/\n" % (len(dossiers), d))

    # a one-line tier tally to stderr so a human sees the shape without parsing JSON
    tally = {"A": 0, "B": 0, "C": 0}
    for one in dossiers:
        tally[one["priority"]] = tally.get(one["priority"], 0) + 1
    sys.stderr.write("priority: A=%d  B=%d  C=%d\n" % (tally["A"], tally["B"], tally["C"]))


if __name__ == "__main__":
    main()
