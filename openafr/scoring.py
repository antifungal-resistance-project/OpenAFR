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
