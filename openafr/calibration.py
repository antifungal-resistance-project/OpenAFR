"""Turn resistance scores into probabilities with a MEASURED reliability (issue #136).

The concordance engine (openafr/concordance.py) answers "does the caller emit the right
panel token?" -- a categorical, right/wrong question. This module answers the next,
harder one the diagnostics contract (docs/DIAGNOSTIC_VERDICT_CONTRACT.md) defers to #136:
*when the engine attaches a probability of resistance, is that number reliable?* A
probability of 0.8 is only honest if, across the isolates we stamp 0.8, about 80% are in
fact resistant. This module fits that mapping on a training split and MEASURES whether it
holds on a held-out test split, so the "calibrated probability" claim is backed, not asserted.

Two families of calibrator, because the two regimes the contract names will NOT share a curve:

  * STRATIFIED empirical-rate calibration (`stratified_rates` / `apply_rates`). The map is a
    lookup: P(resistant | stratum) = resistant carriers / carriers in that stratum on the
    training split, each with a Wilson 95% CI (reused from openafr.backtest). This is the
    honest calibrator for BOTH regimes TODAY:
      - known-variant regime  -> stratum = the specific panel token (Y132F, K143R, ...);
      - novel-variant regime  -> stratum = the structural evidence tier
        (direction x confidence), because openafr/structural.py DELIBERATELY emits a
        *direction*, never a magnitude -- there is no continuous score to fit a sigmoid to,
        so a continuous curve would fabricate a precision the structural module refuses to
        claim. Calibrating the ordinal tier is the honest ceiling there.

  * CONTINUOUS-score calibration (`fit_isotonic` / `fit_platt`). Isotonic regression (PAVA,
    monotone, non-parametric) and Platt scaling (a 1-D logistic with Platt's pseudo-count
    targets so it does not overfit small n). These are the textbook "score -> probability"
    tools #136 names; they apply the day a genuinely continuous resistance score exists
    (e.g. a pooled-C.-albicans MIC regression, #132 option C), and are unit-tested now so
    that path is ready rather than asserted.

Reliability is measured the same way regardless of calibrator: on the held-out predictions,
`evaluate_calibration` reports the Brier score, a binned reliability curve (mean predicted
vs observed frequency per bin), the expected calibration error (ECE), and calibration-in-the-
large (mean predicted vs overall observed rate). A low Brier with a diagonal reliability curve
is the evidence behind the word "calibrated."

Design discipline, matching openafr/concordance.py:
  * PURE + OFFLINE, stdlib only. It consumes already-scored per-isolate records -- it does
    NOT fetch panels, run callers, or import numpy. The orchestration that produces scores
    from a real panel lives in the (data/tool-gated) scripts/validate_calibration.py, exactly
    as concordance's scoring lives in validate_fks1_concordance.py. This keeps the math in the
    unit-tested core (tests/test_calibration.py).
  * FIT AND EVALUATE ARE SEPARATE CALLS. A calibrator is fit on training records and then
    applied to a disjoint test split; nothing here hides the split or evaluates on the fit
    data. The caller (and the frozen pre-registration, work/PREREGISTRATION_calibration.md)
    owns the split, the per-drug / per-regime stratification, and the seed.
  * UNRESOLVED / UNSCORED ROWS ARE EXCLUDED, NEVER GUESSED. A record without a usable score or
    label (score None, label None) drops out of fit and evaluation and is surfaced as a count
    -- the same "excluded, not counted negative" rule as concordance.evaluate and
    backtest.panel_prevalence. A test isolate whose stratum was never seen in training has no
    honest probability to assign; it is reported as `n_unseen_stratum`, not forced to a guess.
"""
import math

from openafr.backtest import wilson_interval


def _clean_labels(labels):
    """Coerce a label to a 0/1 int; anything else (None, '') -> None (excluded)."""
    out = []
    for y in labels:
        if y is None:
            out.append(None)
        elif isinstance(y, bool):
            out.append(1 if y else 0)
        elif isinstance(y, (int, float)):
            out.append(1 if y >= 0.5 else 0)
        else:
            s = str(y).strip().upper()
            if s in ("R", "1", "TRUE", "RESISTANT"):
                out.append(1)
            elif s in ("S", "0", "FALSE", "SUSCEPTIBLE"):
                out.append(0)
            else:
                out.append(None)
    return out


def _paired(scores, labels):
    """Drop rows where the score or the label is missing; return (xs, ys, n_excluded)."""
    ys = _clean_labels(labels)
    xs, keep, n_excluded = [], [], 0
    for x, y in zip(scores, ys):
        if x is None or y is None:
            n_excluded += 1
            continue
        xs.append(float(x))
        keep.append(y)
    return xs, keep, n_excluded


# --- Continuous-score calibrators ---------------------------------------------------------

class _StepModel:
    """A monotone non-decreasing step function score -> probability (isotonic PAVA result).

    Prediction is piecewise-constant with linear interpolation between the fitted knot
    scores, clamped to the end knots outside the training range (never extrapolated past
    [0, 1])."""
    def __init__(self, xs, ys):
        self.xs = list(xs)       # knot scores, sorted ascending
        self.ys = list(ys)       # calibrated probability at each knot, non-decreasing

    def predict_one(self, x):
        xs, ys = self.xs, self.ys
        if not xs:
            return None
        if x <= xs[0]:
            return ys[0]
        if x >= xs[-1]:
            return ys[-1]
        # binary search for the bracketing knots, then linear interpolate
        lo, hi = 0, len(xs) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if xs[mid] <= x:
                lo = mid
            else:
                hi = mid
        x0, x1, y0, y1 = xs[lo], xs[hi], ys[lo], ys[hi]
        if x1 == x0:
            return max(y0, y1)
        frac = (x - x0) / (x1 - x0)
        return y0 + frac * (y1 - y0)

    def predict(self, scores):
        return [self.predict_one(float(x)) if x is not None else None for x in scores]


def fit_isotonic(scores, labels):
    """Isotonic (monotone non-decreasing) calibration via Pool-Adjacent-Violators.

    Assumes a higher score means a higher resistance probability. Fits the monotone step
    function that minimises squared error to the 0/1 labels -- non-parametric, so it does not
    impose a sigmoid shape, at the cost of needing more data than Platt to be stable. Rows
    with a missing score or label are excluded. Returns a model with `.predict(scores)`."""
    xs, ys, _ = _paired(scores, labels)
    if not xs:
        return _StepModel([], [])
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    sx = [xs[i] for i in order]
    sy = [float(ys[i]) for i in order]

    # PAVA: each block is [sum_y, weight, value]; merge while monotonicity is violated.
    blocks = [[y, 1.0, y] for y in sy]
    i = 0
    while i < len(blocks) - 1:
        if blocks[i][2] <= blocks[i + 1][2]:
            i += 1
            continue
        s = blocks[i][0] + blocks[i + 1][0]
        w = blocks[i][1] + blocks[i + 1][1]
        blocks[i:i + 2] = [[s, w, s / w]]
        if i > 0:
            i -= 1

    # expand block values back to per-point, then collapse to knots at unique scores
    fitted, bi, used = [], 0, 0.0
    for b in blocks:
        for _ in range(int(round(b[1]))):
            fitted.append(b[2])
    knot_x, knot_y = [], []
    j = 0
    while j < len(sx):
        k = j
        while k + 1 < len(sx) and sx[k + 1] == sx[j]:
            k += 1
        # average the fitted values across tied scores (they share a PAVA value anyway)
        knot_x.append(sx[j])
        knot_y.append(sum(fitted[j:k + 1]) / (k - j + 1))
        j = k + 1
    return _StepModel(knot_x, knot_y)


class _PlattModel:
    """P(resistant | score) = sigmoid(A * score + B)."""
    def __init__(self, A, B):
        self.A = A
        self.B = B

    def predict_one(self, x):
        z = self.A * x + self.B
        # numerically stable logistic
        if z >= 0:
            return 1.0 / (1.0 + math.exp(-z))
        e = math.exp(z)
        return e / (1.0 + e)

    def predict(self, scores):
        return [self.predict_one(float(x)) if x is not None else None for x in scores]


def fit_platt(scores, labels, max_iter=100, tol=1e-7):
    """Platt scaling: fit a 1-D logistic score -> probability by regularised ML.

    Uses Platt's (1999) pseudo-count targets -- t+ = (N+ + 1)/(N+ + 2), t- = 1/(N- + 2) --
    instead of hard 0/1, so a perfectly separable tiny sample does not drive the slope to
    infinity (the small-n over-confidence the structural module also guards against). Solved
    by Newton-Raphson on the cross-entropy; falls back to the base-rate intercept if the
    Hessian is singular. Rows with a missing score/label are excluded."""
    xs, ys, _ = _paired(scores, labels)
    n_pos = sum(1 for y in ys if y == 1)
    n_neg = len(ys) - n_pos
    if not xs or n_pos == 0 or n_neg == 0:
        # no signal to fit a slope: constant model at the observed base rate
        rate = (n_pos / len(ys)) if ys else 0.0
        rate = min(max(rate, 1e-6), 1 - 1e-6)
        return _PlattModel(0.0, math.log(rate / (1 - rate)))

    hi = (n_pos + 1.0) / (n_pos + 2.0)
    lo = 1.0 / (n_neg + 2.0)
    t = [hi if y == 1 else lo for y in ys]

    A, B = 0.0, math.log((n_neg + 1.0) / (n_pos + 1.0)) * -1.0  # sensible intercept start
    for _ in range(max_iter):
        g0 = g1 = h00 = h01 = h11 = 0.0
        for x, ti in zip(xs, t):
            z = A * x + B
            p = 1.0 / (1.0 + math.exp(-z)) if z >= 0 else math.exp(z) / (1.0 + math.exp(z))
            d = p - ti
            w = max(p * (1.0 - p), 1e-12)
            g0 += d * x
            g1 += d
            h00 += w * x * x
            h01 += w * x
            h11 += w
        det = h00 * h11 - h01 * h01
        if abs(det) < 1e-12:
            break
        dA = (h11 * g0 - h01 * g1) / det
        dB = (h00 * g1 - h01 * g0) / det
        A -= dA
        B -= dB
        if abs(dA) < tol and abs(dB) < tol:
            break
    return _PlattModel(A, B)


# --- Stratified empirical-rate calibrator (the honest calibrator for both regimes today) ---

def stratified_rates(strata, labels):
    """Calibration map by lookup: per stratum, P(resistant) = carriers resistant / carriers.

    `strata` is a per-isolate stratum key (a panel token for the known regime, a structural
    evidence tier for the novel regime); `labels` the matching R/S phenotype. Returns
    {stratum: {k, n, point, ci95}} with a Wilson 95% CI per stratum. Rows with a missing
    stratum or label are excluded. This is the calibrator the novel regime uses, because
    openafr/structural.py emits an ordinal direction, not a continuous score."""
    counts = {}
    ys = _clean_labels(labels)
    for stratum, y in zip(strata, ys):
        if stratum is None or y is None:
            continue
        slot = counts.setdefault(str(stratum), [0, 0])
        slot[1] += 1
        if y == 1:
            slot[0] += 1
    rates = {}
    for stratum, (k, n) in counts.items():
        rates[stratum] = {"k": k, "n": n,
                          "point": (k / n if n else None),
                          "ci95": wilson_interval(k, n)}
    return rates


def apply_rates(rates, strata):
    """Assign each isolate its stratum's training probability. An unseen stratum -> None.

    Returns (probs, n_unseen). A None probability means the training split never contained
    that stratum, so there is no honest calibrated number for it -- the caller excludes those
    isolates from reliability (the 'excluded, not guessed' rule), it does NOT default them to
    the base rate."""
    probs, n_unseen = [], 0
    for stratum in strata:
        rec = rates.get(str(stratum)) if stratum is not None else None
        if rec is None or rec["point"] is None:
            probs.append(None)
            n_unseen += 1
        else:
            probs.append(rec["point"])
    return probs, n_unseen


# --- Reliability metrics (measured on the held-out predictions) ---------------------------

def brier_score(probs, labels):
    """Mean squared error of predicted probability vs the 0/1 outcome (lower = better).

    Rows with a missing probability or label are excluded. Returns (score, n); score is None
    when nothing is scorable."""
    ys = _clean_labels(labels)
    se, n = 0.0, 0
    for p, y in zip(probs, ys):
        if p is None or y is None:
            continue
        se += (p - y) ** 2
        n += 1
    return (se / n if n else None), n


def reliability_curve(probs, labels, n_bins=10):
    """Binned reliability: for each of n_bins equal-width probability bins, the mean predicted
    probability, the observed resistant frequency (with a Wilson CI), and the bin count.

    Perfect calibration puts every bin on the diagonal (mean predicted == observed). Empty
    bins are omitted. Rows with a missing probability/label are excluded."""
    ys = _clean_labels(labels)
    bins = [[0.0, 0, 0] for _ in range(n_bins)]   # [sum_pred, n_pos, n]
    for p, y in zip(probs, ys):
        if p is None or y is None:
            continue
        idx = min(int(p * n_bins), n_bins - 1)
        bins[idx][0] += p
        bins[idx][1] += (1 if y == 1 else 0)
        bins[idx][2] += 1
    curve = []
    for i, (sp, npos, n) in enumerate(bins):
        if n == 0:
            continue
        curve.append({
            "bin": i,
            "range": (i / n_bins, (i + 1) / n_bins),
            "mean_predicted": sp / n,
            "observed": npos / n,
            "ci95": wilson_interval(npos, n),
            "n": n,
        })
    return curve


def expected_calibration_error(probs, labels, n_bins=10):
    """ECE: the count-weighted mean gap |mean_predicted - observed| across the bins. 0 = the
    reliability curve is exactly on the diagonal. Returns None when nothing is scorable."""
    curve = reliability_curve(probs, labels, n_bins=n_bins)
    total = sum(b["n"] for b in curve)
    if total == 0:
        return None
    return sum(b["n"] * abs(b["mean_predicted"] - b["observed"]) for b in curve) / total


def calibration_in_the_large(probs, labels):
    """(mean predicted probability, overall observed rate, n). A large gap means the model is
    biased high or low on average even if its ranking is fine. Rows missing p/label excluded."""
    ys = _clean_labels(labels)
    sp, npos, n = 0.0, 0, 0
    for p, y in zip(probs, ys):
        if p is None or y is None:
            continue
        sp += p
        npos += (1 if y == 1 else 0)
        n += 1
    if n == 0:
        return {"mean_predicted": None, "observed": None, "n": 0}
    return {"mean_predicted": sp / n, "observed": npos / n, "n": n}


def evaluate_calibration(probs, labels, n_bins=10):
    """Bundle the reliability evidence for one held-out split: Brier, reliability curve, ECE,
    calibration-in-the-large, and the scorable / excluded counts. This is the object a
    calibration run reports and the pre-registered bar (work/PREREGISTRATION_calibration.md)
    is applied to."""
    brier, n_scored = brier_score(probs, labels)
    n_excluded = sum(1 for p, y in zip(probs, _clean_labels(labels))
                     if p is None or y is None)
    return {
        "n_scored": n_scored,
        "n_excluded": n_excluded,
        "brier": brier,
        "ece": expected_calibration_error(probs, labels, n_bins=n_bins),
        "calibration_in_the_large": calibration_in_the_large(probs, labels),
        "reliability_curve": reliability_curve(probs, labels, n_bins=n_bins),
    }


def _fmt_ci(ci):
    return f"[{ci[0]:.1%}, {ci[1]:.1%}]" if ci else "[n/a]"


def format_report(summary):
    """Human-readable reliability report of an evaluate_calibration() result (no side effects)."""
    lines = [
        f"scored: {summary['n_scored']} isolate(s); {summary['n_excluded']} excluded "
        f"(no score/label/unseen stratum -- not guessed)",
    ]
    if summary["brier"] is None:
        lines.append("  Brier: n/a (nothing scorable)")
        return "\n".join(lines)
    cil = summary["calibration_in_the_large"]
    lines += [
        f"  Brier score:            {summary['brier']:.4f}  (lower = better)",
        f"  expected calib. error:  {summary['ece']:.1%}",
        f"  calibration-in-large:   predicted {cil['mean_predicted']:.1%} vs "
        f"observed {cil['observed']:.1%} (n={cil['n']})",
        "  reliability curve (predicted -> observed [95% CI], n):",
    ]
    for b in summary["reliability_curve"]:
        lines.append(
            f"    {b['range'][0]:.1f}-{b['range'][1]:.1f}: "
            f"pred {b['mean_predicted']:.1%} -> obs {b['observed']:.1%} "
            f"{_fmt_ci(b['ci95'])}  (n={b['n']})")
    return "\n".join(lines)
