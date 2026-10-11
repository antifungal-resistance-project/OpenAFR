# Quickstart: interpret a *C. auris* isolate (for lab users)

**You have a sequenced *Candida auris* isolate and want to know whether it carries a known
echinocandin- or azole-resistance marker.** This page is the five-minute path to one verdict per
drug-class — no developer setup assumed. For the full output specification see
[`DIAGNOSTIC_VERDICT_CONTRACT.md`](DIAGNOSTIC_VERDICT_CONTRACT.md).

> **Research use only (RUO) — not a clinical test.** This tool reports *concordance with named
> ERG11/FKS1 resistance markers*. It does **not** measure an MIC, does **not** emit a
> probability, and does **not** cover non-target mechanisms (efflux, ERG3, promoter/tandem
> repeats). It must not be used to guide treatment.

## 1. What it does

You give it the genotype of one *C. auris* isolate — at one of two genes — and it returns a
single **categorical verdict per drug-class**:

| Gene | Drug class | What a marker there means |
|---|---|---|
| **FKS1** | echinocandin (caspofungin, micafungin, anidulafungin) | hot-spot mutation → reduced echinocandin susceptibility |
| **ERG11** | azole (fluconazole, voriconazole, …) | pocket mutation → reduced azole susceptibility |

It checks **only the named marker panels** for those two genes. Its honesty boundary *is* the
product: it tells you exactly which markers it did and did not find, and says "I don't know"
out loud rather than guessing (see [§4](#4-how-to-read-the-verdict)).

## 2. Install

The fastest input — a list of variant tokens from your own variant caller — needs **only Python
3 from the standard library**, no external bioinformatics tools:

```bash
conda env create -f environment.yml     # or: mamba env create -f environment.yml
conda activate openafr
```

(The consensus-FASTA paths in [§3](#3-run-it-on-your-isolate) use the same environment.)

## 3. Run it on your isolate

Pick whichever input shape your sequencing pipeline already produces.

**A. Variant tokens — the easiest path.** If your variant caller already reports protein
substitutions, just name them:

```bash
python scripts/interpret_isolate.py --isolate MY_ISOLATE --gene FKS1 --variants S639P
```

Type the isolate for **both** drug-classes in one call by repeating `--gene` with its own
variants (an empty `--variants` is a valid "nothing found"):

```bash
python scripts/interpret_isolate.py --isolate MY_ISOLATE \
    --gene FKS1  --variants S639P \
    --gene ERG11 --variants Y132F
```

If a marker region was **not sequenced**, declare the gap so it is not read as wild-type:

```bash
python scripts/interpret_isolate.py --isolate MY_ISOLATE --gene FKS1 --variants --uncalled 639
```

**B. A consensus CDS FASTA.** A single-record, in-frame consensus of the gene:

```bash
python scripts/interpret_isolate.py --isolate MY_ISOLATE --gene ERG11 --cds erg11_consensus.fasta
```

**C. Per-window FKS1 consensus** (the robust path — what `scripts/recall_fks1.py recall`
produces from reads):

```bash
python scripts/interpret_isolate.py --isolate MY_ISOLATE \
    --gene FKS1 --window HS1=hs1.fasta --window HS2=hs2.fasta --window HS3=hs3.fasta
```

Supply **all three** hot-spot windows (HS1, HS2, HS3); a window you leave out is reported
`missing` and forces `UNRESOLVED` at its markers — the tool never assumes an un-read window is
wild-type.

Add `--json` for a machine-readable report (for a LIMS / pipeline) instead of the human brief.

**Starting from raw reads?** Build the consensus first with the FKS1 orchestrator
(`python scripts/recall_fks1.py recall …`, needs sra-tools/minimap2/samtools ≥1.13), then feed
its per-window output to path **C** above. This CLI starts at genotype/consensus by design.

### See it work with no data of your own

```bash
python scripts/demo_recallers.py --check
```

This runs both callers on a tiny committed sample set (a wild-type control + one known
resistance mutation per gene) and checks them against a committed oracle — the on-ramp before
you point it at your own isolate.

## 4. How to read the verdict

Every drug-class gets exactly one of four values. Three of them are *not* a susceptibility call.

| Verdict | Meaning | What it is **not** |
|---|---|---|
| **`RESISTANCE_MARKER_DETECTED`** | A known, high-PPV resistance mutation is present. | Still not an MIC — it is concordance with the named panel. |
| **`UNCHARACTERIZED_VARIANT`** | A real change is present but it is **not** a trusted resistance marker — the explicit "I don't know". (ERG11 adds a mechanism-based structural best-guess, flagged *calibrated-low*.) | Not a resistance call, and not a clearance. |
| **`NO_KNOWN_MARKER`** | Only wild-type across the checked windows. | **NOT "susceptible".** Only the named ERG11/FKS1 markers were checked — not efflux or other mechanisms. |
| **`UNRESOLVED`** | The window could not be read honestly (coverage gap or indel). | **NOT a negative.** It is excluded from every denominator, never scored as susceptible. |

The two traps to internalize: **`NO_KNOWN_MARKER` ≠ susceptible**, and **`UNRESOLVED` ≠
negative.** A low-PPV FKS1 token (e.g. D642Y) is recognized but *abstains* — it is never read as
resistance. Full enum and provenance fields: [`DIAGNOSTIC_VERDICT_CONTRACT.md`](DIAGNOSTIC_VERDICT_CONTRACT.md).

## 5. Limits (read before acting)

- **RUO**, not a clinical determination — no MIC, no calibrated probability.
- The **echinocandin (FKS1) arm is the validated one**: against measured phenotype it hits
  very-major error 0% and major error 2.2% at a ~10% honest-abstention price. The **azole
  (ERG11) arm is underpowered** and carried as a declared limit. The calibrated-probability
  track is blocked on non-public data. Full write-up:
  [`../work/PREPRINT_diagnostics_engine.md`](../work/PREPRINT_diagnostics_engine.md).
- Coverage is **US-centric** at the surveillance level; this per-isolate tool is unaffected, but
  the population picture on the [resistance-weather page](../README.md#openafr-resistance-weather--the-surveillance-product)
  reflects that bias.
