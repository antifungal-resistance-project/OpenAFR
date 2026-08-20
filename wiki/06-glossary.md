# Glossary

One-line definitions for the terms scattered through the codebase and the rest of this wiki.
Grouped by area. For the fuller story behind the biology terms, see the
[primer](01-antifungal-resistance-primer.md).

## Biology & chemistry

- **Antifungal resistance** — a fungus surviving a drug that used to kill it, usually because a
  mutation changed the drug's target.
- **Azole** — the main class of antifungal drug (fluconazole, itraconazole, posaconazole,
  voriconazole…). Works by jamming the enzyme CYP51.
- ***Candida auris* / *Candidozyma auris*** — a lethal, drug-resistant, hospital-spreading
  yeast; a WHO critical-priority pathogen. The target organism. (Second name is a recent
  reclassification; same organism.)
- **CYP51** — the fungal enzyme azoles jam; makes ergosterol. The *protein*. Same thing as
  ERG11 (the gene).
- **ERG11** — the *gene* that encodes CYP51. In DNA/mutation contexts this project says ERG11;
  in structure contexts, CYP51. Same thing.
- **Ergosterol** — the essential sterol in fungal cell membranes (the fungal analog of
  cholesterol). CYP51 is needed to make it.
- **Heme** — an iron-containing chemical group at the core of CYP51 (same family as the heme in
  blood). Its **iron atom (Fe)** is where azoles bind.
- **N–Fe distance** — the distance from a ligand's closest **nitrogen (N)** to the heme
  **iron (Fe)**. Short = the molecule is coordinating the iron like a real azole. **The
  measurement track 1 ranks by.**
- **Residue** — one position (one amino acid) in a protein chain, e.g. "residue 132."
- **Wild-type** — the natural, unmutated form of a gene or protein. Mutations are described
  relative to it.
- **Mutation shorthand (e.g. Y132F)** — wild-type letter + position + new letter: "at 132,
  Tyrosine (Y) became Phenylalanine (F)." The four this project tracks: **Y132F, K143R, F126L,
  V125A**.
- **Panel / resistance panel** — the set of known azole-resistance mutations the code tags
  (V125A, F126L, Y132F, K143R).
- **Substitution** — a mutation swapping one amino acid for another (a "point substitution").
- **CDS (coding sequence)** — the stretch of DNA that codes for the protein; the re-caller works
  on the ERG11 CDS.
- **Codon** — a group of three DNA bases that codes for one amino acid.
- **Consensus (sequence)** — the single best-supported sequence assembled from many overlapping
  raw sequencing reads of one isolate.
- **Isolate** — one sample of the fungus taken from one patient/place/time.
- **Clade** — a genetic lineage/family of *C. auris*; different clades carry different
  signature mutations.
- **TR34 / TR-type** — a *tandem-repeat / promoter* resistance mechanism (extra copies making
  more enzyme). **Not** a point substitution — the code refuses to map it as one.
- **Efflux pump** — a resistance mechanism where the fungus pumps the drug back out (e.g. via
  *TAC1B*). Out of scope here.
- **Selectivity / cytotoxicity** — killing the fungus without harming human cells; the central
  difficulty of antifungals (fungi are eukaryotes, like us).
- **MIC (minimum inhibitory concentration)** — the standard wet-lab measure of how much drug it
  takes to stop a microbe; the first real yes/no a computational hit must survive.

## Structure & docking (track 1)

- **PDB** — a file format (and the public Protein Data Bank) storing 3D atomic coordinates of a
  molecule. `5TZ1` is the CYP51 structure used here.
- **5TZ1** — the specific *C. albicans* CYP51 crystal structure (with heme + a bound azole) this
  project docks against. Used as a stand-in because no *C. auris* structure existed.
- **Receptor** — the target protein in docking (CYP51). **Ligand** — the small candidate
  molecule being fitted into it.
- **Docking** — simulating how a ligand best fits into a receptor's pocket. Done here with
  **AutoDock Vina**.
- **Vina / AutoDock Vina** — the docking program. The `vina` CLI.
- **PDBQT** — Vina's input/output structure format (PDB + charge + atom-type). Parsed by
  `openafr/pdbqt.py`.
- **Pose** — one predicted position+orientation of a ligand in the pocket. Docking yields
  several per molecule.
- **(Affinity / docking) score** — Vina's own estimate of binding strength (kcal/mol, more
  negative = stronger). **This project deliberately distrusts it** — it ranked worse than random
  here.
- **SMILES** — a compact text encoding of a molecule's structure; the `.smi` ligand files.
- **Active** — a molecule known to work (a real azole). **Decoy** — a molecule presumed inactive,
  included to make the test hard.
- **Held-out** — actives locked away and never used to form the hypothesis, so the final test is
  blind.
- **Property-matched decoys** — decoys matched to actives on size / greasiness / heavy-atom count
  so the method can't "win" by a shortcut (see the circular-decoy and size-shortcut traps).
- **cLogP** — a computed measure of a molecule's greasiness (lipophilicity); one of the matched
  properties.
- **Applicability domain** — the region of chemical space the validation actually covers here,
  ≤ 45 heavy atoms; outside it the method isn't trusted.
- **Induced fit** — a receptor flexing to accommodate a ligand. This pipeline uses a *rigid*
  receptor (a stated limitation).
- **Repacker** — a tool that rebuilds a mutated side chain's atoms. The pinned environment has
  none, so only atom-*deletion* mutations can be built deterministically.

## Metrics & statistics

- **AUC** — Area Under the ROC Curve: overall ranking quality, 0–1. 0.5 = coin flip, 1.0 =
  perfect. **The durable headline number** (0.794).
- **EF@x% (enrichment factor)** — how many times more actives land in the top x% of the ranked
  list than random. **EF@1%** is the razor-top number that proved fragile under multi-seed
  docking.
- **BEDROC** — an AUC variant weighting the top of the list heavily (α=20 here); what matters
  when you can only test the top few.
- **Permutation test** — shuffle the labels many times to see how often chance alone beats your
  result; gave p ≈ 0.0028.
- **Bootstrap CI** — a confidence interval built by resampling the data; the AUC's was
  0.590–0.946.
- **Wilson (score) interval** — a confidence interval for a proportion that stays well-behaved at
  small n and near 0/1 (unlike the textbook Wald interval); used for the azole event frequency.
- **Azole event frequency / panel prevalence** — of the isolates the re-caller could *resolve*,
  the fraction carrying a known panel mutation. The eventual deliverable number of track 2;
  "NOT MEASURED YET" until the re-caller runs.

## Project-specific concepts

- **Track 1 / Track 2** — drug discovery / genomic early-warning surveillance.
- **The gate** — a program (`validate_gate*.py`) that grades a run against a frozen pass bar and
  **exits nonzero (fails)** if the method can't re-discover known drugs. The enforcement arm of
  the founding rule.
- **Pre-registration** — writing the rules + pass thresholds and SHA-256-freezing them *before*
  generating data, so they can't be tuned to pass. Files: `PREREGISTRATION_*.md`.
- **Frozen protocol** — `protocol.yaml`, the hash-pinned docking box + pass bar read by both the
  docking script and the gate (so what's docked can't drift from what's graded).
- **Hash-freezing** — pinning a file's SHA-256 in code and re-checking it at run time; a changed
  file makes any "pass" print as "not credible."
- **Uncalled / resolved / pending** — a `resistance_source` status. An **uncalled** isolate is
  never assumed susceptible; only **resolved** isolates enter a prevalence denominator; pending /
  partial / failed are reported, never counted negative.
- **Emergence signal** — the detector's output that a mutation is *first-appearing* or *rising*.
- **Backlog (`possible_backlog`)** — a fake-looking spike caused by a lab depositing old samples
  now; flagged via the gap between collection date and NCBI-visible date.
- **So-what priority** — ACT-NOW / WATCH / CONTEXT, a triage tier derived from emergence +
  structural verdicts (never invented).
- **No-op (delivery)** — a run that produces nothing because there's no genuine new news. The
  intended common case.
- **Evidence tier / confidence** — a structural verdict is `measured` (a real dock exists),
  `estimate` (a direction, no magnitude), or `none` (honest no-call).
- **ERG11 re-caller** — the module (`openafr/recaller.py` + `scripts/recall_erg11.py`) that turns
  raw reads into a resistance call. **The single blocker for track 2.**
- **NCBI Pathogen Detection** — the public feed track 2 reads; a metadata + genome-pointer feed
  with *no* resistance calls for *C. auris*.
- **SRA** — the Sequence Read Archive; where NCBI stores the raw sequencing reads the re-caller
  would download.
- **Snapshot / baseline / diff** — one normalized NCBI pull; the isolates seen up to a date; the
  change between two pulls.

## Tooling

- **conda / mamba** — the package/environment managers used to install the pinned toolchain
  (`environment.yml`).
- **rdkit** — the only third-party Python import in the code: molecule handling + the scoring
  math.
- **openbabel (`obabel`)** — CLI the docking scripts use to convert molecule formats.
- **fasterq-dump / minimap2 / samtools** — the bioinformatics binaries the re-caller's
  reads→consensus step shells out to (Linux/x86, via bioconda).
- **pytest / pytest-cov** — the test runner and the coverage gate CI enforces.
