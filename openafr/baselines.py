"""Multi-baseline comparison scaffold (issue #83) — the tools a skeptic would reach for.

Pipeline stage: validate/ (baseline comparison, work/PREREGISTRATION_baselines_multi.md)

Run 2 showed the geometric iron-approach criterion beats *AutoDock Vina's* score on the same
poses, and the GNINA pre-registration (work/PREREGISTRATION_baselines.md, #63) set up the one
metal-aware learned opponent. Issue #83 widens that to a *table* of independent baselines run
on the IDENTICAL held-out poses, so the criterion is quantified against real alternatives, not
just vanilla Vina.

This module is the pure, stdlib-only half: it turns each opponent's saved score files into
per-pose ranking keys (lower = better, so it plugs into openafr.scoring.report unchanged) and
describes each opponent as a BaselineSpec. It is deliberately split from the runner — the
scorer binaries are Linux/x86 native (the same wall the ERG11 re-caller and GNINA hit), so
they run on a cloud host (scripts/run_*_rescore.sh) while parsing + grading stay reproducible
on any host (scripts/rescore_baselines.py). Parsing here runs anywhere, including where the
binaries cannot.

Every parser returns a list of per-pose ranking keys IN INPUT-FILE ORDER, so it aligns 1:1
with openafr.pdbqt.read_poses() over the same .pdbqt. A per-pose key of None means that pose
was unscored (never guessed); molecule_key() takes the best (minimum) non-None key, or None if
the molecule has none — which openafr.scoring.report ranks LAST, never dropping it.
"""
import collections

from openafr.gnina import parse_score_only

BaselineSpec = collections.namedtuple(
    "BaselineSpec", "key label subdir ext metal_aware parse citation")


def _gnina_keys(text):
    """GNINA --score_only stdout -> per-pose keys. CNNaffinity is higher-is-better, so the
    ranking key is its negation (lower = better)."""
    return [None if p.get("cnnaffinity") is None else -p["cnnaffinity"]
            for p in parse_score_only(text)]


def parse_plants_features(text):
    """PLANTS `--mode rescore` writes features.csv, one row per scored pose in input order,
    with a header naming the columns. TOTAL_SCORE (ChemPLP, lower = better binding) is the
    ranking column. Returns per-pose TOTAL_SCORE in row order; a row with a missing or
    non-numeric TOTAL_SCORE yields None for that pose rather than a guess.

    ChemPLP carries explicit metal-acceptor terms, which is why PLANTS is the second
    *metal-aware* opponent (alongside GNINA) named in the pre-registration.
    """
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return []
    header = [h.strip() for h in lines[0].split(",")]
    if "TOTAL_SCORE" not in header:
        return []
    col = header.index("TOTAL_SCORE")
    keys = []
    for row in lines[1:]:
        cells = row.split(",")
        if col >= len(cells):
            keys.append(None)
            continue
        try:
            keys.append(float(cells[col]))
        except ValueError:
            keys.append(None)
    return keys


# The pre-registered opponents (work/PREREGISTRATION_baselines_multi.md). Each is a drop-in
# rescorer of the frozen run-2 poses; only the scorer changes. Ordered as reported.
REGISTRY = [
    BaselineSpec(
        key="gnina", label="GNINA CNNaffinity", subdir="work/gnina_scores", ext=".gnina",
        metal_aware=True, parse=_gnina_keys,
        citation="McNutt et al. 2021, GNINA CNN rescoring (--score_only, default model)"),
    BaselineSpec(
        key="plants", label="PLANTS ChemPLP", subdir="work/plants_scores", ext=".csv",
        metal_aware=True, parse=parse_plants_features,
        citation="Korb et al. 2009, PLANTS/ChemPLP (--mode rescore, metal-acceptor terms)"),
]


def molecule_key(per_pose_keys):
    """Best (minimum) non-None ranking key over a molecule's poses, or None if it has none."""
    vals = [k for k in per_pose_keys if k is not None]
    return min(vals) if vals else None


def spec_by_key(key):
    for s in REGISTRY:
        if s.key == key:
            return s
    raise KeyError(key)
