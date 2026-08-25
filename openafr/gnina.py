"""Parser for GNINA `--score_only` output — the metal-aware baseline scorer.

Pipeline stage: validate/ (baseline comparison, work/PREREGISTRATION_baselines.md)

Why this exists
---------------
Run 2 showed the geometric iron-approach criterion beats *AutoDock Vina's* score on
held-out azoles. Vina is known to mishandle metal coordination, so the honest next
question is whether the criterion also beats a *competent, metal-aware* scorer. GNINA's
CNN rescoring is that opponent. This module turns GNINA's `--score_only` stdout into
per-pose numbers the existing gate math can rank — nothing else.

Deliberately split from the runner: parsing is pure and stdlib-only so it is testable and
runnable anywhere (including hosts that cannot run the GNINA binary). The actual
`gnina` invocation lives in scripts/run_gnina_rescore.sh, which saves stdout per ligand;
this module reads those saved files. Poses are emitted in input-file order, so a GNINA
block list aligns 1:1 with openafr.pdbqt.read_poses() over the same .pdbqt.

GNINA `--score_only` prints one block per input pose, e.g.

    Affinity: -7.12 (kcal/mol)
    CNNscore: 0.9123
    CNNaffinity: 6.234
    CNNvariance: 0.451

Blank lines / banner text between blocks are ignored. A block is complete once a
CNNaffinity line is seen; Affinity/CNNscore default to None if GNINA omitted them.
"""
import re

_AFFINITY = re.compile(r"^Affinity:\s*(-?\d+\.?\d*)")
_CNNSCORE = re.compile(r"^CNNscore:\s*(-?\d+\.?\d*)")
_CNNAFFINITY = re.compile(r"^CNNaffinity:\s*(-?\d+\.?\d*)")


def parse_score_only(text):
    """Parse GNINA --score_only stdout into a list of per-pose dicts, in file order.

    Each dict has keys 'affinity' (Vina-style kcal/mol, lower better, or None),
    'cnnscore' (pose-quality prob, higher better, or None) and 'cnnaffinity'
    (predicted affinity, higher better, or None). One dict per input pose; the
    CNNaffinity line closes a pose. Returns [] if no pose block is present.
    """
    poses, cur = [], {}
    for line in text.splitlines():
        line = line.strip()
        m = _AFFINITY.match(line)
        if m:
            cur["affinity"] = float(m.group(1))
            continue
        m = _CNNSCORE.match(line)
        if m:
            cur["cnnscore"] = float(m.group(1))
            continue
        m = _CNNAFFINITY.match(line)
        if m:
            cur["cnnaffinity"] = float(m.group(1))
            poses.append({"affinity": cur.get("affinity"),
                          "cnnscore": cur.get("cnnscore"),
                          "cnnaffinity": cur["cnnaffinity"]})
            cur = {}
    return poses


def best_cnnaffinity(poses, keep=None):
    """Return the highest CNNaffinity over `poses`, or None.

    `keep`: optional list of bools, same length/order as `poses`, selecting which poses
    are eligible (used to restrict to iron-bound poses for mode GC). None = all poses.
    Higher CNNaffinity = better; callers rank by the NEGATED value (lower = better) so it
    plugs into openafr.scoring.report unchanged.
    """
    vals = [p["cnnaffinity"] for i, p in enumerate(poses)
            if p.get("cnnaffinity") is not None and (keep is None or keep[i])]
    return max(vals) if vals else None
