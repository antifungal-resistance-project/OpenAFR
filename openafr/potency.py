"""Continuous-potency correlation core — does the N-Fe geometry track *how potent* a molecule is?

Why this module exists
----------------------
Every validated headline so far is a *binary* active/inactive claim: the geometry
criterion enriches actives over inactives at AUC ~0.79-0.85 (work/RESULTS_run2_holdout.md,
work/RESULTS_verified_inactives.md). A stronger claim a wet lab could bet on is *monotonic*:
that a tighter nitrogen-to-heme-iron approach tracks a *higher measured affinity*, not merely
membership in the "active" label (issue #81). This module holds the two things a test can pin
for that stronger claim:

  1. The **aggregation rule** (`median_pchembl`) that turns one compound's several ChEMBL
     bioactivity records into a single continuous potency value, so the correlation has one
     row per compound and no compound is weighted by how many times it was assayed.
  2. The **rank statistics** (`spearman`, `permutation_p`, `bootstrap_ci`) the grader uses to
     score geometry-vs-potency, implemented in numpy only (the pinned env has no scipy).

The network/I/O shell that pulls ChEMBL and applies the applicability triage lives in
`scripts/make_potency_set.py`; the grader that reads the frozen dataset + the docked geometry
and computes the correlation lives in `scripts/validate_gate_potency.py`. Both import from here
so the rules a reviewer must trust sit in one tested place, exactly as `openafr.inactives` holds
the verified-inactive evidence rule.

What pChEMBL is (and why it is the right continuous axis)
--------------------------------------------------------
ChEMBL's `pchembl_value` is -log10(molar potency) — IC50/Ki/Kd/EC50 folded onto one
comparable -log scale (pChEMBL 7 = 100 nM, 8 = 10 nM). Higher = more potent. It only exists
for records with an exact ('=') molar-standardised value, so a censored ('>'/'<') or
non-molar (raw ug/mL MIC) record carries none and is dropped by construction. That is the
honest filter: we correlate geometry only against potencies ChEMBL itself standardised.

The hypothesis direction, fixed here so the grader cannot flip it: a *smaller* closest N-Fe
distance (tighter iron approach) should go with a *larger* pChEMBL (higher potency), i.e. a
**negative** Spearman rho. The grader tests the one-sided 'less' alternative.
"""
import math

import numpy as np

# The standard_types whose pchembl_value we treat as one comparable binding/functional
# potency axis. Frozen here, mirrored into the pre-registration; adding a type is a
# construction change, not a tuning knob. MIC is deliberately absent: ChEMBL does not emit
# a pchembl_value for raw MIC records (they are not reliably molar-standardised), so a MIC
# never reaches this rule in the first place.
POTENCY_TYPES = ("IC50", "Ki", "Kd", "EC50", "Potency")


def median_pchembl(records, allowed_types=POTENCY_TYPES):
    """One compound's ChEMBL records -> (median_pchembl, n_used, sorted_types_used).

    Only records that carry a numeric `pchembl_value`, an allowed `standard_type`, and an
    exact relation (`standard_relation` '=' or absent) are used; everything else is ignored.
    Returns (None, 0, ()) if the compound has no usable quantitative potency record — such a
    compound is not placed on the continuous axis at all (it is not a zero).
    """
    vals, types = [], set()
    for rec in records:
        st = rec.get("standard_type")
        if st not in allowed_types:
            continue
        rel = rec.get("standard_relation")
        if rel not in (None, "=", "'='"):
            continue
        raw = rec.get("pchembl_value")
        if raw in (None, ""):
            continue
        try:
            v = float(raw)
        except (TypeError, ValueError):
            continue
        vals.append(v)
        types.add(st)
    if not vals:
        return None, 0, ()
    return float(np.median(vals)), len(vals), tuple(sorted(types))


def rankdata(a):
    """Average-rank transform (ties share their mean rank), matching scipy's default."""
    a = np.asarray(a, dtype=float)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a), dtype=float)
    ranks[order] = np.arange(1, len(a) + 1, dtype=float)
    # Resolve ties to the average rank within each tied group.
    _, inv, counts = np.unique(a, return_inverse=True, return_counts=True)
    sums = np.zeros(len(counts))
    np.add.at(sums, inv, ranks)
    return (sums / counts)[inv]


def _pearson(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xd, yd = x - x.mean(), y - y.mean()
    denom = math.sqrt(float(np.dot(xd, xd)) * float(np.dot(yd, yd)))
    if denom == 0.0:
        return float("nan")
    return float(np.dot(xd, yd) / denom)


def spearman(x, y):
    """Spearman rank correlation = Pearson correlation of the average ranks."""
    if len(x) != len(y):
        raise ValueError("x and y must be the same length")
    return _pearson(rankdata(x), rankdata(y))


def permutation_p(x, y, *, n_perm=100000, seed=42, alternative="less"):
    """One-sided permutation p-value for the Spearman rho of (x, y).

    The null shuffles y against x; the p-value is the fraction of permuted rho at least as
    extreme as observed in the pre-registered direction. 'less' tests rho <= observed (the
    tighter-N-Fe/higher-potency hypothesis, a NEGATIVE rho). A +1 numerator/denominator
    (Davison-Hinkley) keeps the estimate unbiased and never reports p=0.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    obs = spearman(x, y)
    rx = rankdata(x)
    ry = rankdata(y)
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(ry)
        rho = _pearson(rx, perm)
        if alternative == "less":
            count += rho <= obs
        elif alternative == "greater":
            count += rho >= obs
        else:  # two-sided
            count += abs(rho) >= abs(obs)
    return (count + 1) / (n_perm + 1)


def bootstrap_ci(x, y, *, n_boot=10000, seed=42, alpha=0.05):
    """Percentile bootstrap CI for the Spearman rho, resampling (x, y) pairs with replacement."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    rng = np.random.default_rng(seed)
    rhos = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        rhos[i] = spearman(x[idx], y[idx])
    rhos = rhos[~np.isnan(rhos)]
    lo = float(np.quantile(rhos, alpha / 2))
    hi = float(np.quantile(rhos, 1 - alpha / 2))
    return lo, hi
