"""File-drawer bound — how much can a 'no-precedent' label be hiding? (issue #80)

Why this module exists
----------------------
openafr.precedent already says, loudly and in prose, that a `no-precedent` verdict means
"no deposited enzyme record was found", NOT "never tested" — because a failing measurement
is often never deposited (publication bias, the *file-drawer problem*, Rosenthal 1979). That
caveat is correct but qualitative. This module makes it *systematic*: it encodes a small,
assumption-explicit model of the residual risk so every `no-precedent` output can carry a
quantified caveat, and so the payoff of adding sources (issue #78) is expressed as a number
rather than a hope.

The full derivation, the literature grounding, and the honest limitations live in the
referenced note work/METHODS_file_drawer.md. This module is that note made executable and
testable, so the two cannot drift.

What is (and is NOT) being bounded — read this before using a number here
-------------------------------------------------------------------------
The thing a reader is tempted to want is P(this compound has quietly failed | no-precedent).
We deliberately DO NOT return that, because it factors as

    P(hidden prior failure | no-precedent)
        = prevalence  π   × leak  L
        = P(a prior negative CYP51 measurement exists for this compound)
          × P(that measurement is absent from every queried source | it exists)

and π is unknowable from our data — it depends entirely on how much prior attention the
compound has had, which varies molecule to molecule. Reporting a single P(failed) would be
false precision dressed as rigour. What the file-drawer literature *does* let us reason about
is the leak L: given a negative result was produced, how often does it evade the stores we
query? That — and only that — is what this module quantifies, and it returns it as a BAND,
never a point.

The model of the leak L
-----------------------
Let d = probability that a given negative CYP51 measurement is *discoverable* in at least one
queried store (deposited into ChEMBL's curated literature or PubChem BioAssay's HTS dumps).
For a single independent source with deposit probability d_i the miss probability is (1 - d_i);
independent sources compound, so

    L = Π_i (1 - d_i).

Two facts fix the plausible range of d, both from the note:
  * Negative/null results are substantially under-published. The file-drawer problem and its
    successors put the fraction of conducted work that reaches the literature well below 1 —
    commonly summarised as roughly half, with negatives lagging further. So for the
    *literature* channel d is modest.
  * PubChem BioAssay is a targeted mitigant: HTS depositors upload the whole plate, actives
    AND inactives, so for the HTS subset a negative IS discoverable. This is exactly why
    openafr.precedent queries PubChem at all. It lifts d, but only for compounds that happened
    to sit on a deposited screen.

Because the mix of "measured in a targeted (often unpublished) study" vs "swept up by a
deposited HTS campaign" is itself unknown, d cannot be pinned to a point. The note adopts an
explicit, clearly-labelled *illustrative* band d ∈ [0.40, 0.70] for the two queried sources
taken together, giving L ∈ [0.30, 0.60]: of the negatives that exist, an estimated 30–60%
would evade both ChEMBL and PubChem today. These numbers are assumption-driven, not measured
for CYP51 (that measurement is impossible — it is the file drawer), and are here to size the
effect and to be argued with, not to be quoted as fact.

The one honest lever: more sources
----------------------------------
L = Π(1 - d_i) shrinks geometrically as independent sources are added. Adding BindingDB /
patents / full-text (issue #78) at even a modest per-source deposit probability measurably
lowers the leak — which is precisely why #78 raises confidence. `leak_after_adding_sources`
makes that payoff computable so #78 can be justified with a number.
"""
from __future__ import annotations

# Illustrative deposit-probability band for the TWO currently-queried sources
# (ChEMBL + PubChem) taken together. Assumption-driven, grounded in the under-publication
# literature; see work/METHODS_file_drawer.md. NOT a measured CYP51 rate.
DEPOSIT_RATE_LO = 0.40
DEPOSIT_RATE_HI = 0.70

METHODS_NOTE = "work/METHODS_file_drawer.md"


def leak(deposit_rate: float) -> float:
    """Fraction of *existing* negative measurements that evade a source with this deposit rate.

    Returns 1 - deposit_rate. This is the conditional leak L for one channel; it is NOT the
    probability a compound has failed (that needs the unknowable prevalence π — see module
    docstring). `deposit_rate` must be a probability in [0, 1].
    """
    if not 0.0 <= deposit_rate <= 1.0:
        raise ValueError(f"deposit_rate must be a probability in [0, 1], got {deposit_rate}")
    return 1.0 - deposit_rate


def residual_hidden_band(lo: float = DEPOSIT_RATE_LO, hi: float = DEPOSIT_RATE_HI):
    """The illustrative leak band (L_lo, L_hi) for the currently-queried sources.

    Returned as a (low, high) pair — a BAND, by design. There is no `residual_hidden_point`:
    a single number here would be false precision, because the true deposit rate for
    undeposited negatives is, definitionally, unmeasurable. `lo`/`hi` are deposit rates with
    lo <= hi; the returned leaks are the reverse-ordered (1 - hi, 1 - lo).
    """
    if not 0.0 <= lo <= hi <= 1.0:
        raise ValueError(f"need 0 <= lo <= hi <= 1, got lo={lo}, hi={hi}")
    # higher deposit rate -> lower leak, so the band flips.
    return (leak(hi), leak(lo))


def leak_after_adding_sources(new_source_rates, base_band=(DEPOSIT_RATE_LO, DEPOSIT_RATE_HI)):
    """Leak band after adding independent sources with the given deposit rates (issue #78).

    Independent sources compound: L = Π_i (1 - d_i). `new_source_rates` is an iterable of
    per-source deposit probabilities for the sources being ADDED to the two already queried,
    whose combined rate is `base_band`. Returns the new (low, high) leak band. This is the
    executable form of #78's payoff: each source added drives the leak down geometrically.
    """
    base_lo, base_hi = base_band
    if not 0.0 <= base_lo <= base_hi <= 1.0:
        raise ValueError(f"base_band must be 0 <= lo <= hi <= 1, got {base_band}")
    lo_leak = leak(base_hi)   # best case for coverage -> smallest leak
    hi_leak = leak(base_lo)   # worst case for coverage -> largest leak
    for d in new_source_rates:
        m = leak(d)
        lo_leak *= m
        hi_leak *= m
    return (lo_leak, hi_leak)


def caveat_line(band=None) -> str:
    """A one-line, quantified caveat for a `no-precedent` label, safe to print inline.

    States the residual honestly: it is the fraction of *existing* negatives that would evade
    the queried sources, NOT a probability the compound failed, and it is illustrative. Cites
    the note so a reader can find the assumptions.
    """
    lo, hi = residual_hidden_band() if band is None else band
    return (f"'no-precedent' is not 'untested': of negative CYP51 results that exist, an "
            f"illustrative {round(lo * 100)}-{round(hi * 100)}% would evade the queried sources "
            f"(file-drawer leak; assumption-driven, not a P(failed) — see {METHODS_NOTE}).")
