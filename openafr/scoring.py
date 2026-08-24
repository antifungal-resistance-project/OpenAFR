"""Ranking metrics — the math the entire OpenAFR result rests on.

Two functions, locked by tests/test_gate_metrics.py:

    metrics(ordered_flags)  pure scoring on a list of activity flags ALREADY
                            sorted best-first -> (AUC, [EF@1/5/10%], BEDROC)

    report(label, ranked)   takes (name, key, is_active) rows, lower key = better,
                            key None = unrankable. Sorts, appends unrankable LAST
                            (never drops them), then calls metrics() and prints.

All three metrics come from rdkit.ML.Scoring, the canonical implementations:
    AUC    - overall ranking power. 0.5 = coin flip, 1.0 = perfect.
    EF@x%  - enrichment: how many times more actives in the top x% than random.
    BEDROC - AUC that weights the top of the list heavily (alpha=20), which is
             what matters when you can only afford to test the top few hundred.

EF_FRACTIONS and BEDROC_ALPHA are the metric *definitions*; the tunable pass bar
(min AUC / min EF) lives in the hash-frozen protocol.yaml, not here.
"""
import random

from rdkit.ML.Scoring import Scoring

EF_FRACTIONS = [0.01, 0.05, 0.10]
BEDROC_ALPHA = 20.0


def metrics(ordered_flags, ef_fractions=EF_FRACTIONS, bedroc_alpha=BEDROC_ALPHA):
    """Score a list of activity flags (1=active, 0=decoy) already ordered best-first."""
    r = [[1 if a else 0] for a in ordered_flags]
    return (Scoring.CalcAUC(r, 0),
            [Scoring.CalcEnrichment(r, 0, [f])[0] for f in ef_fractions],
            Scoring.CalcBEDROC(r, 0, bedroc_alpha))


def report(label, ranked, gating=False, ef_fractions=EF_FRACTIONS, bedroc_alpha=BEDROC_ALPHA):
    """ranked: list of (name, key, is_active), lower key = better. None key = ranked last."""
    ok = [r for r in ranked if r[1] is not None]
    bad = [r for r in ranked if r[1] is None]
    ordered = sorted(ok, key=lambda r: r[1]) + bad          # unrankable go last
    flags = [r[2] for r in ordered]
    auc, efs, bedroc = metrics(flags, ef_fractions, bedroc_alpha)
    n_act = sum(flags)
    print(f"\n{label}{'   [PRIMARY — gates the result]' if gating else '   (secondary, reported only)'}")
    print(f"  ranked           : {len(ordered)}  ({n_act} actives, {len(ordered)-n_act} decoys)")
    if bad:
        nb = sum(1 for r in bad if r[2])
        print(f"  unrankable       : {len(bad)} molecules ({nb} actives) — placed LAST per pre-registration")
    print(f"  AUC              : {auc:.3f}")
    for f, ef in zip(ef_fractions, efs):
        print(f"  EF@{int(f*100):>3}%          : {ef:.2f}x")
    print(f"  BEDROC (a=20)    : {bedroc:.3f}")
    top = [r[0] for r in ordered[:15] if r[2]]
    print(f"  actives in top 15: {len(top)}/{n_act}  {top if top else ''}")
    return auc, efs, bedroc


def _decoys_before(ordered_flags):
    """Per active, how many decoys are ranked ahead of it. (n_decoys, [count per active])"""
    seen_dec, per_active = 0, []
    for flag in ordered_flags:
        if flag:
            per_active.append(seen_dec)
        else:
            seen_dec += 1
    return seen_dec, per_active


def auc_from_positions(n_decoys, decoys_before):
    """Mann-Whitney AUC straight from active positions: the mean fraction of decoys each
    active outranks. Identical to Scoring.CalcAUC on untied data (locked by a test), but it
    takes positions rather than a flag list, so actives can be resampled for a bootstrap."""
    if not decoys_before or n_decoys == 0:
        return 0.0
    return sum(n_decoys - d for d in decoys_before) / (len(decoys_before) * n_decoys)


def permutation_p(ordered_flags, n_shuffles=20000, seed=42):
    """One-sided permutation test on AUC: how often does a random relabelling of the SAME
    ranking do at least as well? Reported as (1 + #{perm >= observed}) / (1 + n_shuffles) --
    the add-one form, so a p of exactly 0 is never claimed from a finite number of shuffles.
    Deterministic under the given seed."""
    rng = random.Random(seed)
    observed = metrics(ordered_flags)[0]
    flags = list(ordered_flags)
    at_least = 0
    for _ in range(n_shuffles):
        rng.shuffle(flags)
        n_dec, before = _decoys_before(flags)
        if auc_from_positions(n_dec, before) >= observed:
            at_least += 1
    return (1 + at_least) / (1 + n_shuffles)


def bootstrap_auc_ci(ordered_flags, n_boot=10000, seed=42, alpha=0.05):
    """Percentile bootstrap CI for AUC, resampling the ACTIVES with replacement and holding
    the decoy pool fixed -- the small-n set is the actives, and they are what the interval is
    about. Returns (lo, hi). Deterministic under the given seed."""
    rng = random.Random(seed)
    n_dec, before = _decoys_before(ordered_flags)
    if not before or n_dec == 0:
        return (0.0, 0.0)
    draws = sorted(auc_from_positions(n_dec, [rng.choice(before) for _ in before])
                   for _ in range(n_boot))
    lo = draws[int(alpha / 2 * n_boot)]
    hi = draws[min(n_boot - 1, int((1 - alpha / 2) * n_boot))]
    return lo, hi
