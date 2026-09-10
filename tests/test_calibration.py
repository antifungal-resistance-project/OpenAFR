"""Tests for the probability-calibration engine (openafr/calibration.py, issue #136).

These lock the guarantees the "calibrated probability with known reliability" claim leans on:
  * the stratified empirical-rate calibrator is the honest per-variant / per-tier lookup with
    Wilson CIs, and an UNSEEN test stratum gets None (excluded, never defaulted to a guess);
  * unresolved / unscored rows (score None, label None) are EXCLUDED from fit and from every
    reliability denominator, and surfaced as a count -- never scored as a negative;
  * isotonic calibration is monotone non-decreasing and recovers the right shape;
  * Platt scaling recovers a known sigmoid and does not blow up on a separable tiny sample;
  * Brier / reliability curve / ECE / calibration-in-the-large are the exact arithmetic, so a
    well-calibrated prediction scores near-zero ECE and a biased one is caught.
No network, no reads, no panel: everything runs against synthetic scored isolates.
"""
import math

from openafr import calibration as cal


# --- label coercion + exclusion discipline -----------------------------------------------

def test_clean_labels_coerces_rs_bool_and_excludes_unknown():
    assert cal._clean_labels(["R", "S", True, False, 1, 0]) == [1, 0, 1, 0, 1, 0]
    assert cal._clean_labels([None, "", "maybe"]) == [None, None, None]


def test_paired_drops_missing_score_or_label():
    xs, ys, n_ex = cal._paired([0.1, None, 0.9, 0.5], ["S", "R", None, "R"])
    assert xs == [0.1, 0.5]
    assert ys == [0, 1]
    assert n_ex == 2


# --- stratified empirical-rate calibrator (the regime calibrator) -------------------------

def test_stratified_rates_are_per_stratum_fractions_with_wilson_ci():
    strata = ["Y132F", "Y132F", "Y132F", "K143R", "K143R"]
    labels = ["R", "R", "S", "S", "S"]
    rates = cal.stratified_rates(strata, labels)
    assert rates["Y132F"]["k"] == 2 and rates["Y132F"]["n"] == 3
    assert math.isclose(rates["Y132F"]["point"], 2 / 3)
    assert rates["K143R"]["point"] == 0.0
    lo, hi = rates["Y132F"]["ci95"]
    assert 0.0 <= lo <= rates["Y132F"]["point"] <= hi <= 1.0


def test_stratified_rates_exclude_missing_rows():
    rates = cal.stratified_rates(["Y132F", None, "Y132F"], ["R", "R", None])
    assert rates["Y132F"]["n"] == 1 and rates["Y132F"]["k"] == 1


def test_apply_rates_unseen_stratum_is_none_not_guessed():
    rates = cal.stratified_rates(["Y132F", "Y132F"], ["R", "S"])
    probs, n_unseen = cal.apply_rates(rates, ["Y132F", "NOVEL_T123I"])
    assert probs[0] == 0.5           # seen stratum -> training rate
    assert probs[1] is None          # unseen -> no honest number
    assert n_unseen == 1


# --- isotonic calibration -----------------------------------------------------------------

def test_isotonic_is_monotone_non_decreasing():
    # scores correlate with outcome but noisily; PAVA must still be non-decreasing
    scores = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    labels = [0, 0, 0, 1, 0, 1, 1, 1, 1]
    model = cal.fit_isotonic(scores, labels)
    preds = model.predict(scores)
    assert all(b >= a - 1e-9 for a, b in zip(preds, preds[1:]))
    assert 0.0 <= min(preds) and max(preds) <= 1.0


def test_isotonic_clamps_outside_training_range():
    model = cal.fit_isotonic([0.2, 0.8], [0, 1])
    assert model.predict([-5.0])[0] == model.predict([0.2])[0]
    assert model.predict([5.0])[0] == model.predict([0.8])[0]


def test_isotonic_recovers_a_clean_threshold():
    scores = [0.1, 0.2, 0.3, 0.7, 0.8, 0.9]
    labels = [0, 0, 0, 1, 1, 1]
    model = cal.fit_isotonic(scores, labels)
    assert model.predict([0.15])[0] < 0.5 < model.predict([0.85])[0]


# --- Platt scaling ------------------------------------------------------------------------

def test_platt_recovers_a_known_sigmoid():
    # generate labels from a steep logistic in score; Platt should rank/shape it back
    scores = [i / 50.0 - 1.0 for i in range(101)]          # -1.0 .. +1.0
    labels = [1 if s > 0 else 0 for s in scores]
    model = cal.fit_platt(scores, labels)
    assert model.predict([-0.8])[0] < 0.2
    assert model.predict([0.8])[0] > 0.8
    assert abs(model.predict([0.0])[0] - 0.5) < 0.2
    assert model.A > 0                                      # higher score -> higher P


def test_platt_does_not_diverge_on_separable_tiny_sample():
    # perfectly separable n=4: hard-label ML would push the slope to infinity; the
    # pseudo-count targets keep the probabilities strictly inside (0, 1).
    model = cal.fit_platt([0.0, 0.1, 0.9, 1.0], [0, 0, 1, 1])
    preds = model.predict([0.0, 1.0])
    assert all(0.0 < p < 1.0 for p in preds)
    assert preds[0] < preds[1]


def test_platt_all_one_class_is_constant_base_rate():
    model = cal.fit_platt([0.1, 0.5, 0.9], [0, 0, 0])
    assert model.A == 0.0
    assert model.predict([0.1])[0] < 0.5


# --- reliability metrics ------------------------------------------------------------------

def test_brier_score_is_mean_squared_error_and_excludes_missing():
    # (0.0-0)^2 + (1.0-1)^2 = 0 over 2 scorable rows; the None-prob row is excluded
    score, n = cal.brier_score([0.0, 1.0, None], [0, 1, 1])
    assert n == 2 and score == 0.0
    score2, _ = cal.brier_score([0.5, 0.5], [0, 1])
    assert math.isclose(score2, 0.25)


def test_perfect_calibration_has_near_zero_ece():
    # 100 isolates at p=0.3, exactly 30 resistant -> bin is on the diagonal
    probs = [0.3] * 100
    labels = [1] * 30 + [0] * 70
    ece = cal.expected_calibration_error(probs, labels, n_bins=10)
    assert ece < 1e-9


def test_biased_predictions_have_large_ece():
    # claim 0.9 for everyone but only 10% are resistant -> ECE ~ 0.8
    probs = [0.9] * 100
    labels = [1] * 10 + [0] * 90
    ece = cal.expected_calibration_error(probs, labels, n_bins=10)
    assert ece > 0.7


def test_reliability_curve_bins_and_omits_empty():
    probs = [0.05, 0.15, 0.95]
    labels = [0, 0, 1]
    curve = cal.reliability_curve(probs, labels, n_bins=10)
    assert [b["bin"] for b in curve] == [0, 1, 9]
    assert curve[0]["n"] == 1 and curve[0]["observed"] == 0.0
    assert curve[-1]["observed"] == 1.0


def test_calibration_in_the_large_compares_mean_pred_to_observed():
    cil = cal.calibration_in_the_large([0.8, 0.8, 0.8, 0.8], [1, 0, 0, 0])
    assert math.isclose(cil["mean_predicted"], 0.8)
    assert math.isclose(cil["observed"], 0.25)
    assert cil["n"] == 4


def test_evaluate_calibration_bundles_and_counts_excluded():
    probs = [0.3] * 10 + [None]
    labels = [1] * 3 + [0] * 7 + [1]
    summary = cal.evaluate_calibration(probs, labels, n_bins=10)
    assert summary["n_scored"] == 10
    assert summary["n_excluded"] == 1
    assert summary["ece"] < 1e-9
    assert "scored: 10" in cal.format_report(summary)


def test_format_report_handles_nothing_scorable():
    summary = cal.evaluate_calibration([None], [None])
    assert "n/a" in cal.format_report(summary)


# --- end-to-end: fit on train, MEASURE reliability on held-out test -----------------------

def test_fit_train_then_evaluate_held_out_test_is_well_calibrated():
    # a stratified fit on a balanced train set, applied to a disjoint test set drawn from the
    # same per-stratum rates, should land near the diagonal on the held-out split.
    train_strata = (["Y132F"] * 40 + ["K143R"] * 40)
    train_labels = ([1] * 30 + [0] * 10) + ([1] * 10 + [0] * 30)   # 0.75 and 0.25
    rates = cal.stratified_rates(train_strata, train_labels)
    assert math.isclose(rates["Y132F"]["point"], 0.75)
    assert math.isclose(rates["K143R"]["point"], 0.25)

    test_strata = (["Y132F"] * 20 + ["K143R"] * 20)
    test_labels = ([1] * 15 + [0] * 5) + ([1] * 5 + [0] * 15)      # matches train rates
    probs, n_unseen = cal.apply_rates(rates, test_strata)
    assert n_unseen == 0
    summary = cal.evaluate_calibration(probs, test_labels, n_bins=10)
    assert summary["ece"] < 1e-9
