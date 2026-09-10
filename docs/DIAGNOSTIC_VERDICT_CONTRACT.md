# The diagnostic verdict contract (issue #133) — concordance-only, day one

The output spec for the resistance-interpretation engine (`.context/DIAGNOSTICS_CONCEPT.md`).
It defines *exactly* what the engine emits per isolate and — just as importantly — what it
refuses to claim. This contract **is** the product: the honest boundary is the differentiator.

**Day-one scope is fixed by #132.** The calibration-grade panel does not exist in public form
for *C. auris* (class-imbalanced, ERG11-hotspot-skewed; see
`work/RESULTS_diagnostics_panel_spike.md`). So this contract is **concordance-only**: it detects
*named resistance markers*, with no calibrated resistance probability. The probability field is
deferred to the pooled-*C. albicans* calibration track (#136), and re-enters this contract only
after that track produces a validated calibration — never before.

Built entirely from primitives the existing callers (`openafr/recaller.py` ERG11,
`openafr/fks1_caller.py`) already emit. Nothing here invents a new schema; it names and
constrains what the callers produce.

## The verdict object (per isolate, per drug-class)

An isolate yields one verdict per drug-class it was typed for: **azole** (from ERG11) and/or
**echinocandin** (from FKS1). Each verdict is:

```
{
  drug_class:      "azole" | "echinocandin",
  gene:            "ERG11" | "FKS1",
  resolution:      "resolved" | "unresolved",
  resolution_note: <reason string when unresolved: coverage_gap | in_window_indel_refused
                    | orchestration_failed>,          # excluded, NOT a negative
  called_tokens: [                                     # verbatim, non-synonymous only
    { token: "Y132F", class: "known_resistance" },
    { token: "T123I", class: "uncharacterized" },      # the honest "I don't know" slot (#135)
    ...
  ],
  verdict:         <one of the four enum values below>,
  structural:      <ERG11 pocket verdict where mappable (#135); null for FKS1 by design>,
  provenance:      { reference_accession, reference_sha256, caller_version,
                     known_panel_version },
  scope:           "RUO — research use only; not a clinical determination"
}
```

## The verdict enum (exactly four values — categorical, never a probability)

1. **`RESISTANCE_MARKER_DETECTED`** — resolved, and ≥1 `called_token` is a **known-resistance
   panel** substitution. ERG11 panel: V125A, F126L, Y132F, K143R (`RESISTANCE_PANEL` in
   `recaller.py`). FKS1 HS1 panel: S639F, S639P, S639Y (`fks1_caller.py`).
2. **`UNCHARACTERIZED_VARIANT`** — resolved, a non-synonymous change is present but **none** is a
   known-resistance panel token. This is the explicit **"I don't know"** verdict — surfaced, not
   hidden. #135 attaches a mechanism-based structural best-guess for ERG11 here; the verdict
   itself stays uncharacterized until a variant is promoted into the panel table.
3. **`NO_KNOWN_MARKER`** — resolved, and only wild-type / no panel token in the resistance
   windows. **This is NOT a susceptibility call** (see the load-bearing rule below).
4. **`UNRESOLVED`** — the window could not be called honestly (coverage gap, in-window indel
   refused, orchestration failure). **Excluded from every denominator, never scored as a
   negative** — the caller's "an uncalled codon is uncalled, never wild-type" rule
   (`recaller.py` honesty constraint #2) carried up to the verdict level.

## Four load-bearing honesty rules (the reason this contract exists)

1. **No probability field, day one.** The verdict is categorical. #132 found no data to
   calibrate a probability; emitting one now would be an unbacked number. Deferred to #136.
2. **`NO_KNOWN_MARKER` ≠ "susceptible."** A concordance tool detects *known markers in
   ERG11/FKS1*. It does not cover efflux (TAC1/MRR1/CDR1), ERG3, or non-target mechanisms, so
   absence of a marker cannot assert susceptibility. The verdict name says exactly what was and
   wasn't checked, and the RUO scope string forbids reading it as a clinical "S".
3. **`UNCHARACTERIZED_VARIANT` is a first-class verdict, not silence.** Lookup-table tools fall
   silent on a novel variant; this contract requires an explicit verdict for it — the moat
   (#135).
4. **`UNRESOLVED` is excluded, not negative.** Same "excluded not guessed" discipline as
   `concordance.evaluate` and `backtest.panel_prevalence`. A caller that honestly refuses to
   guess must never be scored as if it called wild-type.

## Per-drug-class asymmetry (declared, not glossed)

- **ERG11 / azole** carries a `structural` pocket verdict where the variant is mappable (the
  CYP51/heme-iron geometry moat, #135 / `openafr/structural.py`).
- **FKS1 / echinocandin is detection-only, no structural so-what.** Echinocandins don't
  coordinate a metal, so the geometry moat does not transfer (`fks1_caller.py` docstring). Its
  `structural` field is always `null`, and no verdict may borrow the azole pocket rationale.

## Out of scope for this contract (declared, not silently dropped)

- **Calibrated probability / MIC value.** Deferred to #136; not emitted day one.
- **Clinical actionability.** Every verdict is RUO; the go/no-go for a clinical claim is #139.
- **Species identification.** The contract assumes the organism/gene is already typed; species
  ID is a separate upstream concern, not part of the verdict.
- **Non-target resistance mechanisms** (efflux, ERG3, promoter/TR). Their existence is *why*
  rule 2 holds, but calling them is not in this gene-panel contract.

## What #133 delivers next

This spec, then: a `verdict` schema module + unit tests that convert the callers' existing
token/resolution output into the object above, with the four-value enum and the RUO scope
string enforced. #134 wires the callers into it; #135 fills the `UNCHARACTERIZED` structural
best-guess; #137 validates the whole thing as concordance against the C. auris benchmarks.
