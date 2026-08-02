"""Turn a screened, ranked library into a human-readable candidate report.

Pipeline stage: report/  (runs after screen.sh, the last stage before a human reads it)

Why this exists
---------------
rank_candidates.py emits a bare TSV — good for a machine, thin for the person who
has to decide whether any of these molecules is worth wet-lab time. This stage is
the presentation layer: for each candidate it lays out the *evidence* behind its
rank (how close it brings a nitrogen to the heme iron, how many of its poses make
that coordination contact, and — clearly demoted — the docking score), together
with the provenance a skeptical reader needs (which receptor, which frozen
protocol, whether that protocol is still unmodified) and the caveats that keep the
list honest.

It computes nothing new. Ordering is the same validated criterion rank_candidates
uses (minimum N-Fe distance, smaller first, unrankable last); the per-candidate
evidence is analyze_poses.analyze(). This file only arranges what those already
produce into Markdown a collaborator can read.

It does NOT re-run the gate, and it says so, loudly and at the top: a rank here is
a hypothesis for a wet lab, never a hit, and it means nothing unless
validate_gate2.py has PASSED under the frozen protocol first.

Usage: python scripts/report_candidates.py SCREEN_DIR [RECEPTOR_PDB] [-o REPORT.md] [-n TOP]
  SCREEN_DIR    directory of docked *.pdbqt poses (e.g. work/screen)
  RECEPTOR_PDB  receptor with the heme iron       (default work/receptor_A.pdb)
  -n / --top    how many candidates to expand in full detail (default 10)
"""
import argparse
import hashlib
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from openafr import protocol
from openafr.pdbqt import iron_position
from scripts.analyze_poses import FE_CUTOFF, analyze


def collect(screen_dir, receptor_pdb):
    """Return (fe, rows) where rows are per-molecule evidence, best candidate first.

    Each row: {name, best_fe, n_passing, n_poses, best_score}. best_fe is the
    closest nitrogen-to-iron approach over all poses, or None when the molecule
    produced no usable pose. Ordering is the validated criterion — smaller best_fe
    first, unrankable (None) last — identical to rank_candidates.py.
    """
    fe = iron_position(receptor_pdb)
    rows = []
    for f in sorted(os.listdir(screen_dir)):
        if not f.endswith(".pdbqt"):
            continue
        r = analyze(os.path.join(screen_dir, f), fe)
        if r is None:                       # docking produced no poses at all
            rows.append({"name": f[:-6], "best_fe": None, "n_passing": 0,
                         "n_poses": 0, "best_score": None})
            continue
        rows.append({"name": f[:-6], "best_fe": r["best_fe"], "n_passing": r["n_passing"],
                     "n_poses": r["n_poses"], "best_score": r["best_score"]})
    rows.sort(key=lambda x: (x["best_fe"] is None, x["best_fe"] if x["best_fe"] is not None else 0.0))
    return fe, rows


def _fmt(v, spec, na="—"):
    return format(v, spec) if v is not None else na


def render(screen_dir, receptor_pdb, fe, rows, top):
    """Return the Markdown report as one string."""
    digest = hashlib.sha256(protocol.PROTOCOL.read_bytes()).hexdigest()
    verified = digest == protocol.PROTOCOL_SHA
    dock = protocol.load()["docking"]
    gate = protocol.load()["gate"]

    rankable = [r for r in rows if r["best_fe"] is not None]
    unrankable = [r for r in rows if r["best_fe"] is None]

    L = []
    L.append("# OpenAFR candidate report")
    L.append("")
    L.append("Molecules ranked by the **one** criterion the validation gate proved works on "
             "held-out actives: the closest approach of a ligand nitrogen to the heme iron "
             "(**smaller = better**). The docking score is carried along but does **not** "
             "decide order — on this system it ranked worse than random.")
    L.append("")
    L.append("> **This is a list of hypotheses, not hits.** A high rank means one thing: a "
             "wet lab should look. It is meaningful *only* if `validate_gate2.py` has PASSED "
             "under the frozen protocol below — this script does not re-run the gate. Ranking "
             "is valid only within the declared applicability domain (≤ 45 heavy atoms; see "
             "`filter_library.py`).")
    L.append("")

    L.append("## Provenance")
    L.append("")
    L.append(f"- **Receptor:** `{receptor_pdb}` — heme iron at "
             f"({fe[0]:.2f}, {fe[1]:.2f}, {fe[2]:.2f})")
    L.append(f"- **Iron-bound cutoff:** a pose counts as coordinating when N–Fe < "
             f"{FE_CUTOFF:.1f} Å")
    L.append(f"- **Docking protocol:** {dock['size']} Å box, exhaustiveness "
             f"{dock['exhaustiveness']}, {dock['num_modes']} modes, seed {dock['seed']}")
    L.append(f"- **Gate pass bar:** AUC ≥ {gate['min_auc']}, EF@1% ≥ {gate['min_ef1']} "
             "(the thresholds the ranking criterion had to clear to earn its use here)")
    if verified:
        L.append(f"- **Protocol integrity:** ✅ `protocol.yaml` unmodified since freezing "
                 f"(sha256 `{digest[:16]}…`)")
    else:
        L.append(f"- **Protocol integrity:** ⚠️ **`protocol.yaml` HAS CHANGED since it was "
                 f"frozen** (expected `{protocol.PROTOCOL_SHA[:16]}…`, got `{digest[:16]}…`) — "
                 "any ranking under this protocol is not credible until re-frozen.")
    L.append(f"- **Screened:** {len(rows)} molecules from `{screen_dir}` "
             f"({len(rankable)} rankable, {len(unrankable)} produced no usable pose)")
    L.append("")

    L.append("## Ranked candidates")
    L.append("")
    L.append("| rank | ligand | N–Fe (Å) | iron-bound poses | Vina (kcal/mol) |")
    L.append("|---:|---|---:|---:|---:|")
    for i, r in enumerate(rows, 1):
        d = _fmt(r["best_fe"], ".3f", na="no pose")
        poses = f"{r['n_passing']}/{r['n_poses']}" if r["n_poses"] else "0/0"
        sc = _fmt(r["best_score"], ".1f", na="NA")
        L.append(f"| {i} | {r['name']} | {d} | {poses} | {sc} |")
    L.append("")
    L.append("*N–Fe = closest nitrogen-to-iron approach (the deciding criterion). "
             "Iron-bound poses = how many docked poses reach coordination distance. "
             "Vina = affinity, secondary only.*")
    L.append("")

    detail = rankable[:top]
    if detail:
        L.append(f"## Top {len(detail)} in detail")
        L.append("")
        for i, r in enumerate(detail, 1):
            L.append(f"### {i}. {r['name']}")
            L.append("")
            L.append(f"- **N–Fe closest approach:** {r['best_fe']:.3f} Å "
                     "— the validated criterion; this is why it ranks here.")
            if r["n_passing"]:
                L.append(f"- **Coordination:** {r['n_passing']} of {r['n_poses']} poses "
                         f"bring a nitrogen within {FE_CUTOFF:.1f} Å of the iron — consistent "
                         "with the azole mechanism the gate tested.")
            else:
                L.append(f"- **Coordination:** none of its {r['n_poses']} poses reach "
                         f"{FE_CUTOFF:.1f} Å; its rank rests on the closest approach only, "
                         "which is weaker evidence.")
            L.append(f"- **Vina affinity:** {_fmt(r['best_score'], '.1f', na='NA')} kcal/mol "
                     "(secondary — does not decide rank).")
            L.append("")

    L.append("---")
    L.append("*Generated by `scripts/report_candidates.py`. A rank is a hypothesis for a wet "
             "lab, never a hit.*")
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("screen_dir", help="directory of docked *.pdbqt poses")
    ap.add_argument("receptor_pdb", nargs="?", default="work/receptor_A.pdb",
                    help="receptor PDB with the heme iron (default work/receptor_A.pdb)")
    ap.add_argument("-o", "--out", help="write the Markdown report here (default: stdout)")
    ap.add_argument("-n", "--top", type=int, default=10,
                    help="how many candidates to expand in full detail (default 10)")
    args = ap.parse_args()

    fe, rows = collect(args.screen_dir, args.receptor_pdb)
    report = render(args.screen_dir, args.receptor_pdb, fe, rows, args.top)

    if args.out:
        pathlib.Path(args.out).write_text(report)
        print(f"# wrote {args.out} ({len(rows)} molecules)", file=sys.stderr)
    else:
        sys.stdout.write(report)


if __name__ == "__main__":
    main()
