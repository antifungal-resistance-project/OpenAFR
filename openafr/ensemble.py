"""Ensemble aggregation of the mode-C criterion over multiple receptor conformers.

The within-azole ensemble look (work/PREREGISTRATION_ensemble.md) does NOT invent a new
feature. It re-uses mode C — the minimum ligand-nitrogen-to-heme-iron distance, the project's
validated primary criterion — computed once per conformer, then aggregates the per-conformer
values across the exhaustive 3-structure C. albicans CYP51 ensemble by ONE pre-committed rule:

    ensemble-minimum (primary): a molecule's score = the smallest N-Fe distance it achieves in
                                ANY conformer. Mechanistic reading: "does there EXIST an
                                experimentally-observed channel conformation in which this
                                molecule coordinates the iron like a real azole?"
    ensemble-mean (diagnostic): reported, never promoted to the gate (prereg stopping rule).

A per-conformer value is None when the molecule had no nitrogen-bearing pose in that conformer.
None values are ignored by both aggregators; a molecule with None in EVERY conformer is
unrankable (aggregate None) and, per the prereg, ranked LAST by the grader — never dropped.

These are pure functions over already-computed distances; the pose scoring itself stays in
openafr.pdbqt (min_nitrogen_iron_distance). Locked by tests/test_ensemble.py.
"""


def aggregate(per_conformer_values, mode="min"):
    """Aggregate one molecule's per-conformer best N-Fe distances (Angstrom).

    per_conformer_values: iterable that may contain None (no coordinating pose in a conformer).
    Returns the aggregate distance, or None if every conformer was None (unrankable).
    """
    present = [v for v in per_conformer_values if v is not None]
    if not present:
        return None
    if mode == "min":
        return min(present)
    if mode == "mean":
        return sum(present) / len(present)
    raise ValueError(f"unknown aggregation mode: {mode!r} (use 'min' or 'mean')")


def ensemble_rows(names, per_conformer, mode="min", exclude=None):
    """Build (name, aggregate_distance) rows over an ensemble.

    names:          molecule names to score.
    per_conformer:  {conformer_id: {name: best_nfe_or_None}} — one inner dict per conformer.
    exclude:        optional {name: {conformer_id, ...}} of conformers to DROP for that
                    molecule before aggregating (the self-docking guard: a co-crystal azole
                    is not scored against the conformer it was solved in).
    Returns [(name, aggregate)], aggregate None when the molecule is unrankable.
    """
    exclude = exclude or {}
    rows = []
    for name in names:
        skip = exclude.get(name, set())
        vals = [d.get(name) for cid, d in per_conformer.items() if cid not in skip]
        rows.append((name, aggregate(vals, mode)))
    return rows
