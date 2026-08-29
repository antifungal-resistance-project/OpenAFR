# Methods note — bounding the file-drawer limit of `no-precedent` (issue #80)

**What this note is.** A short, referenced bound on how much a `no-precedent` precedent
verdict can be hiding. The precedent check (`openafr/precedent.py`,
`scripts/check_precedent.py`) already says in prose that `no-precedent` means *"no deposited
enzyme record was found"*, not *"never tested"*. This note makes that caveat systematic:
it states a model, grounds it in the publication-bias literature, computes the residual, and
— most importantly — is honest about what cannot be bounded. The model is executable and
tested in `openafr/filedrawer.py` / `tests/test_filedrawer.py`, so the prose and the number
cannot drift.

**Bottom line, first.** A `no-precedent` label does **not** bound the chance a compound has
quietly failed to anything near zero, and no amount of extra querying can drive that to zero —
only reduce it. What we *can* quantify is the *leak*: of the negative CYP51 measurements that
exist, an illustrative **30–60%** would evade the two sources we currently query. This is an
assumption-driven size estimate, not a measured CYP51 rate (that rate is, by definition,
unmeasurable — it is the file drawer). Every `no-precedent` output now carries this caveat
inline via `openafr.filedrawer.caveat_line()`.

## The failure mode

`rank_candidates.py` orders molecules by heme-iron geometry. A high rank on a molecule that a
wet lab *already measured inactive* against fungal CYP51 is a wasted bench experiment. The
precedent check catches such molecules **when the measurement was deposited** in ChEMBL or
PubChem BioAssay. The file-drawer risk is the residue: a real prior failure that was never
deposited, so the check returns `no-precedent` and the reader is tempted to read "novel".

## The model

The quantity a reader wants is `P(hidden prior failure | no-precedent)`. It factors into two
pieces, only one of which we can touch:

```
P(hidden prior failure | no-precedent) = π × L

  π  (prevalence)  = P(a prior negative CYP51 measurement exists for this compound)
  L  (leak)        = P(it is absent from every queried source | it exists)
```

**π is unknowable from our data.** It depends entirely on how much prior attention a specific
molecule has had — a repurposed clinical drug is far more likely to have been screened than an
obscure fragment. Emitting a single `P(failed)` would multiply in a number we simply do not
have and dress an assumption as a measurement. We therefore report **L only**, as a band.

**Modelling the leak L.** Let `d_i` be the probability that a negative CYP51 measurement is
*discoverable* in source *i* (deposited into ChEMBL's curated literature or PubChem BioAssay's
HTS records). A single source misses it with probability `1 - d_i`; independent sources
compound:

```
L = Π_i (1 - d_i)
```

Two facts from the literature fix the plausible range of `d` for the two sources we query:

1. **Negative results are substantially under-published.** The file-drawer problem
   (Rosenthal 1979) and the meta-research that followed it find that a large fraction of
   conducted work — and negatives especially — never reaches the literature; Fanelli (2012)
   documents negatives *disappearing* from the published record over time. So for the
   **literature** channel (the ChEMBL side), `d` is well below 1.
2. **PubChem BioAssay is a deliberate mitigant.** HTS depositors upload the entire plate —
   actives *and inactives* — so a negative that sat on a deposited screen **is** discoverable.
   This is exactly why `openafr.precedent` queries PubChem at all, not only ChEMBL. It lifts
   `d`, but only for the subset of compounds that happened to appear on a deposited screen.

Because the mix of "measured in a targeted (often unpublished) study" versus "swept up by a
deposited HTS campaign" is itself unknown, `d` cannot be pinned to a point. We adopt an
explicit, **illustrative** combined band for the two queried sources:

```
d ∈ [0.40, 0.70]   ⇒   L = 1 - d ∈ [0.30, 0.60]
```

i.e. of the negatives that exist, an estimated **30–60%** evade both ChEMBL and PubChem today.
These endpoints are assumption-driven, chosen to sit within the under-publication literature's
range while crediting PubChem's inactive capture; they are here to *size* the effect and to be
argued with, not to be quoted as a CYP51 fact.

## The one honest lever: more sources (issue #78)

Because `L = Π(1 - d_i)`, the leak falls **geometrically** as independent sources are added.
Adding BindingDB / patents / full-text search (issue #78) at even a modest per-source deposit
probability measurably lowers the leak. `filedrawer.leak_after_adding_sources` makes this
computable. For example, adding one source at `d = 0.5` multiplies the leak by 0.5:

```
current band          L ∈ [0.30, 0.60]
+ 1 source at d=0.5   L ∈ [0.15, 0.30]
+ 2 sources at d=0.5  L ∈ [0.075, 0.15]
```

This is the quantified justification for #78: each independent source roughly halves the
fraction of existing negatives that can hide.

## Limitations (the honest part)

- **π is not modelled.** This note bounds the *leak*, not the probability a compound failed.
  A low leak does not make a compound novel; it makes the *check* more trustworthy.
- **The band is illustrative, not measured.** The true deposit rate for undeposited negatives
  cannot be measured — that is the definition of the file drawer. The `[0.40, 0.70]` deposit
  band is a defensible assumption, not data. `filedrawer` refuses to return a point estimate
  for exactly this reason.
- **Source independence is assumed.** ChEMBL and PubChem overlap (some literature assays are
  also deposited), so `Π(1 - d_i)` slightly *understates* the true leak; the band's width is
  meant to absorb this.
- **Near-neighbour precedent is out of scope here.** `precedent.nearest_tested_neighbour`
  separately catches close analogues of tested-inactive compounds; that reduces π-side risk
  and is not double-counted in L.

## What changed in the outputs

Every `no-precedent` label now carries the quantified caveat inline, sourced from this note:

- `scripts/check_precedent.py` footer prints `filedrawer.caveat_line()`.
- `scripts/report_candidates.py` appends the quantified band to each `no-precedent` row and to
  the report's "Prior testing" preamble.

No output states or implies that `no-precedent` means novel or untested.

## References

- Rosenthal, R. (1979). *The file drawer problem and tolerance for null results.*
  Psychological Bulletin, 86(3), 638–641. — origin of the file-drawer problem.
- Fanelli, D. (2012). *Negative results are disappearing from most disciplines and countries.*
  Scientometrics, 90, 891–904. — under-publication of negative results, worsening over time.
- Kim, S. et al. *PubChem* (Nucleic Acids Research, annual database updates). — PubChem
  BioAssay stores depositor-adjudicated **inactives**, the basis for querying it as a
  file-drawer catch.
- Gaulton, A. et al. *The ChEMBL database.* (Nucleic Acids Research). — curated, quantitative
  literature bioactivity, the on-target verdict source.
