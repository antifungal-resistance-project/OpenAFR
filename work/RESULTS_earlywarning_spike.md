# Spike #20 — can NCBI Pathogen Detection feed a C. auris azole early-warning?

**Question (issue #20):** before building issues #21–27, confirm the killer premise —
that a public feed delivers fresh, dated, geolocated *C. auris* isolates whose
azole-resistance state we can read.

**Probe:** `scripts/probe_ncbi_auris.py` (stdlib only). Run against release
**PDG000000067.670**, last modified **2026-08-12** (the day before this write-up),
**29,111 isolates**. Reproduce: `python scripts/probe_ncbi_auris.py`.

Source of record:
`https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Candidozyma_auris/` — C. auris now
sits under its reassigned genus *Candidozyma*. Each `PDG000000067.<N>/` release
carries `Metadata/*.metadata.tsv` (one row per isolate, 67 columns) and
`AMR/*.amr.metadata.tsv`.

## Verdict: QUALIFIED PASS — but the premise is only half-true

The feed is an **excellent epidemiological backbone** (who / where / when, at scale,
refreshed within days) and a **non-existent resistance-genotype feed**. NCBI does not
call a single azole-resistance determinant for C. auris. The early-warning track is
still viable, but the AMR-calling layer is now unambiguously **ours to build**, not
something we read off NCBI. This does **not** trip the kill criterion (dates,
geography, and cadence are all strong, and mutations are derivable from linked reads),
but it materially reshapes issues #21–27.

## The four go/no-go metrics

| # | Metric | Result | Status |
|---|--------|--------|--------|
| 1 | Collection dates | `collection_date` 99.2%; `target_creation_date` 100%; `sra_release_date` 96.6% | ✅ strong |
| 2 | Geography | `geo_loc_name` (country/region) 95.8%; `lat_lon` 23.5% | ✅ country-level; ⚠️ coords sparse |
| 3 | Azole-resistance calls **exposed** | `AMR_genotypes` **0.0%**; `computed_types` **0.0%**; AMRFinderPlus applied to **0** isolates; ERG11 hot-spot mentions **0**; `AST_phenotypes` **0.1%** (23 isolates) | ❌ **not exposed** |
| 3′| Azole-resistance calls **derivable** | SRA raw reads (`Run`) 96.6%; assembled genome (`asm_acc`) 3.4% | ✅ derivable by us, from reads |
| 4 | Cadence | ~400–800 new DB targets/month through 2026; 670 releases; latest built 2026-08-12 | ✅ days-to-weeks |

## What this means for the track (the load-bearing finding)

**NCBI runs no AMR pipeline on C. auris.** `number_amr_genes = 0`,
`amrfinder_applied = 0`, and `amrfinder_version = NULL` for **all 29,111** isolates.
AMRFinderPlus is a bacteria-first tool; the fungal azole determinants (ERG11/CYP51
`Y132F`, `K143R`, `F126L`, the `VF125AL`/`F126L` clade-III haplotype, tandem
duplications, `TAC1B` GOF, `ERG3`) are simply not in its reference set here. So the
feed's `AMR_genotypes` / `computed_types` columns are permanently empty for our
purpose. **There is no NCBI resistance call to consume.**

Consequences, in priority order:

1. **The ERG11 re-caller (issue #23 / T1) is not optional — it is the product.**
   The one thing that would let us *read* resistance off the feed does not exist
   upstream. We must call ERG11 azole substitutions ourselves from the **SRA reads**
   (96.6% coverage), not from NCBI genotype calls and not from assemblies (only 3.4%
   are assembled). Plan issue #23 around a reads→ERG11 caller.

2. **The backtest (#27) cannot use NCBI labels as ground truth.** Only 23 isolates
   (0.1%) carry any `AST_phenotypes`, and those are phenotype (MIC-derived S/I/R),
   not genotype. That is too thin to validate an emergence signal against. The
   backtest needs an external truth set (published clade/mutation panels, CDC AR
   Isolate Bank strains with known ERG11 status) — pre-register that in #27 rather
   than assuming NCBI provides it.

3. **The feed is US-centric.** 25,705 / 29,111 (88%) are `USA` (largely a CDC/state
   surveillance pipeline). China (284), South Africa (260), Colombia (216), Pakistan
   (165), Qatar (122) trail far behind. "Global early warning" is honestly "US
   early warning + thin international." Say so in the alert copy (#25) and scope
   claims (#24) — do not imply worldwide coverage.

4. **Cadence is genuinely good.** New isolates enter at a few hundred to ~1,500 per
   month with no multi-year lag, and the release folder rebuilt the day before this
   probe. An "early warning" cadence claim is defensible on the *arrival* axis. The
   real lag to characterize (do this in #27's honesty section) is
   **collection_date → target_creation_date** (sample taken → visible in NCBI), which
   this probe does not yet measure; add it before promising a warning *lead time*.

## What the feed DOES give us for free (keep)

- Dated, geolocated isolate stream at scale (28.8k with both date and country).
- Linked SRA run accessions (96.6%) → the substrate for our own ERG11 calling.
- SNP clusters and phylogenetic trees (`Clusters/`, `SNP_trees/`) — clade context for
  distinguishing a *new* resistant lineage from clonal spread of a known one (useful
  for the "emergence vs. noise" logic in #22).

## Bottom line

Do not pivot the source — NCBI Pathogen Detection is the right backbone. But rewrite
the mental model: **it is a metadata + genome-pointer feed, not an AMR feed.** The
resistance signal is manufactured by *our* re-caller over SRA reads, and the structural
so-what (#24) is what differentiates this from a bare novelty stream. Proceed to #21
(snapshot store + baseline) with that framing; make #23 (ERG11 re-caller) the critical
path; and pre-register #27's ground truth externally, not from NCBI.
