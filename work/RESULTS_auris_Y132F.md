# Results — *C. auris* Y132F resistance-pocket transfer test

Date: 2026-08-06. Pre-registered in
[PREREGISTRATION_auris_Y132F.md](PREREGISTRATION_auris_Y132F.md) (sha256 689de2fd…,
verified unmodified at grading time). Protocol hash 6eca16b1… verified unmodified.

> **Outcome: the gate FAILED on the pre-registered bar.** Per the pre-registration this is
> reported as the outcome, unrevised. No new ranking metric is substituted to call it a
> pass. This is a publishable negative — and an informative one.

## What was tested

The run-2 geometric criterion (rank azoles by minimum nitrogen-to-heme-iron distance) was
validated in the **wild-type** *C. albicans* pocket. This run asked the narrower question
the mission actually needs: **does that criterion still work once the pocket carries the
dominant *C. auris* resistance mutation Y132F?**

Controlled single variable: the *same* 7 held-out azole actives and *same* property-matched
decoys as run 2, the *same* frozen docking protocol (26 Å box @ 70.61/66.28/4.18,
exhaustiveness 32, seed 42) — the **only** change is the receptor side chain at position
132 (Tyr → Phe, built by deterministic deletion of the para-hydroxyl;
`work/receptor_auris_Y132F.pdb`, sha256 924cf7b8…). The heme iron is unchanged.

Docking: 355/357 ligands docked; the 2 failures were both **decoys** (for fenticonazole),
so all 7 actives are present and none were dropped.

## The gate result

| mode | AUC | EF@1% | EF@5% | BEDROC | role |
|---|---:|---:|---:|---:|---|
| **C — iron-approach distance** | **0.805** | **0.00x** | 5.63x | 0.235 | **PRIMARY — FAIL** |
| A — Vina score | 0.411 | 0.00x | 0.00x | 0.002 | secondary (worse than random) |
| D — combined rank | 0.677 | 0.00x | 1.41x | 0.040 | secondary |

Pass bar (both required): **AUC ≥ 0.70 AND EF@1% ≥ 5.0x**. AUC clears it; **EF@1% = 0.00x
does not**. Gate = FAIL.

## What actually happened (the texture matters)

The failure is **not** "azoles stopped reaching the iron." Overall separation held — AUC
0.805 is essentially unchanged from wild-type (0.794), and the best drug still coordinates
tightly. The failure is confined to the **very top** of the ranking:

| held-out active | rank (of 355) | N–Fe (Å) | percentile |
|---|---:|---:|---:|
| efinaconazole | 7 | 2.722 | top 2.0% |
| luliconazole | 16 | 2.813 | top 4.5% |
| bifonazole | 42 | 3.085 | 11.8% |
| sertaconazole | 46 | 3.206 | 13.0% |
| oxiconazole | 50 | 3.268 | 14.1% |
| fenticonazole | 79 | 3.551 | 22.3% |
| butoconazole | 264 | 5.342 | 74.4% |

The top 6 slots are all **decoys**, reaching the iron 2.62–2.72 Å — a hair closer than the
best real drug (efinaconazole, 2.722 Å). Since EF@1% only looks at the top ~3–4 molecules,
one thin reshuffle at the tip drops it to zero.

Compared to wild-type run 2 (EF@1% 12.68x), the collapse is stark, but its cause is a
**sub-0.1 Å reshuffle among near-tied molecules** — i.e. it lands right in the
per-molecule noise floor the repurposing run already flagged (the clotrimazole outlier).
The Y132F change was just enough to let property-matched decoys tie and edge past the real
drugs at the extreme top.

## Interpretation (pre-committed, not post-hoc)

Per the pre-registration's FAIL branch: **the WT-validated criterion does not transfer to
the Y132F pocket at the pre-registered bar, and screening a novel library against Y132F is
NOT justified under this protocol.** The wild-type result must not be silently assumed to
hold on the resistant target — which is exactly the risk VISION.md named.

This is a useful negative for two reasons: (1) it is the kind of "does it still work when
the target changes?" check that open docking repos almost never run, and (2) it says
something concrete — a single resistance substitution is enough to erase the top-rank
enrichment that made the tool useful, even while broad separation survives.

## Honest limitations (declared in the pre-registration, all apply)

1. **A rigid single-atom deletion under-models real resistance.** Clinical Y132F also
   involves conformational change, H-bond-network rearrangement, and co-occurring
   efflux/expression changes. This captures the pocket-geometry change only, so the FAIL is
   evidence about *this specific perturbation*, not proof the drugs are clinically dead.
2. **The result lives near the tool's noise floor.** The top-rank reshuffle is sub-0.1 Å;
   at that resolution single poses are not reliable, as the repurposing run showed.
3. **Y132F only.** K143R (clade I) and F126L (clade III) were not tested — they need a
   side-chain repacker not yet in the pinned environment.
4. **Decoys presumed inactive** (unchanged from run 2); 2 decoys failed to dock and are
   simply absent from the pool.

## Reproduce

```
conda activate openafr
python scripts/mutate_receptor.py Y132F -o work/receptor_auris_Y132F.pdb
python scripts/prep_receptor.py --mutant work/receptor_auris_Y132F.pdb \
    work/receptor_auris_Y132F.pdbqt
cat data/ligands/actives_holdout_final.smi data/ligands/decoys_holdout.smi > /tmp/run2_all.smi
python scripts/prep_ligands.py /tmp/run2_all.smi work/run2_ligands
RECEPTOR=work/receptor_auris_Y132F.pdbqt scripts/screen.sh work/run2_ligands work/screen_auris_Y132F
python scripts/validate_gate_auris.py work/screen_auris_Y132F     # exits 0 on pass, 1 on fail
```

---

## Addendum — durable-metric reading (added 2026-08-21; the verdict above is UNCHANGED)

> This section is a *reading* of the numbers already graded above, added later. It does
> **not** re-grade the run and does **not** flip the pre-registered outcome: the gate FAILED
> and stays FAILED. Re-scoring a sha-frozen FAIL into a pass on a different metric is the
> exact p-hacking move pre-registration exists to prevent (refocus plan D5). What this
> addendum does instead is read the same table through the lens of the
> [reliability run](RESULTS_reliability.md) (2026-08-11), which settled *which metrics on
> this decoy set are trustworthy* — a determination made after this run's pre-registration
> was frozen.

**The metric that fails the bar is the one later shown to be unreliable.** The Y132F FAIL is
driven entirely by `EF@1% = 0.00x` (AUC 0.805 clears the `≥ 0.70` half of the bar; only the
`EF@1% ≥ 5.0x` half fails). The reliability run then demonstrated — on the *wild-type* poses,
pre-committed, two robust aggregators — that top-1% enrichment on this property-matched decoy
set is a **single-shot docking artifact, not a durable property of the criterion**: run five
searches and combine them any way and WT EF@1% collapses 12.68x → 0.00x, because with more
sampling the decoys reliably reveal genuine iron-coordinating poses. The project's committed
conclusion there was that **AUC and EF@5% are the defensible metrics and EF@1% is noise.**

**Read on the durable metrics, the criterion transfers to the Y132F pocket essentially intact:**

| durable metric | wild-type (run 2) | Y132F | reading |
|---|---:|---:|---|
| AUC | 0.794 | **0.805** | transfers (unchanged, even marginally higher) |
| EF@5% | 5.63x | **5.63x** | transfers (identical) |
| EF@1% *(discredited)* | 12.68x → 0.00x under resampling | 0.00x | not a durable signal — see reliability run |

So the honest technical reading is: **broad enrichment (a trustworthy top-~5% shortlist)
survives the dominant *C. auris* azole-resistance substitution; razor-top-1% ranking does
not — but that top-1% figure was already known not to survive even in wild-type.** The FAIL
is a FAIL against a bar that leans on a since-discredited metric; the pre-registration's
EF@1% clause is a recorded weakness of *the bar*, not new evidence the tool broke on Y132F.

**The catch that keeps this honest (refocus plan, outside-voice #1 — do NOT overclaim).**
That AUC is *unchanged* by the mutation (0.805 ≈ 0.794) cuts both ways. It means the geometric
criterion is **nearly insensitive to the Y132F change** — so "broad enrichment transfers to
Y132F" is close to a **tautology**, not evidence the tool can pick molecules that *defeat*
resistance. The held-out actives here are the 7 known azoles that Y132F clinically defeats;
ranking defeated drugs well on the mutant pocket demonstrates robustness of enrichment, **not**
resistance-breaker discovery. A genuine resistance-aware claim needs actives that are
*measured active against the resistant strain* (the T1 active-set expansion), not this panel.
This reading must never be pitched as "OpenAFR triages molecules that beat resistance."

**Net for the wet-lab package (T3):** lead with the pre-registered FAIL as evidence of
discipline ("I ran the does-it-still-work-on-the-mutant check almost no docking repo runs, and
reported the negative honestly"), then add this durable-metric reading as nuance — broad
enrichment is robust to Y132F; top-1% ranking is not durable in *either* pocket; and neither
result claims resistance-breaker discovery, which is what the expanded active set is for.
