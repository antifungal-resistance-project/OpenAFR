# Antifungal Resistance — a primer

*Audience: someone who has never studied biology past high school. No code on this page.
By the end you'll understand every biological term the rest of the wiki uses.*

---

## 1. Why fungal infections are a serious problem

Bacteria, viruses, and fungi can all make us sick. We hear constantly about antibiotic
resistance (bacteria) and antivirals (viruses), but **fungi** get far less attention — and
that neglect is now dangerous.

Fungi are a bigger problem than most people realize:

- Invasive fungal infections kill an estimated **~1.5 million people a year** worldwide.
- They mostly strike people whose immune systems are already weakened — cancer patients on
  chemotherapy, transplant recipients, people in intensive care, premature babies.
- We have **very few** antifungal drugs, and research into new ones is dramatically
  underfunded compared to antibacterials.

### The hard part: fungi are a lot like us

Here is the core difficulty that makes antifungals special. Bacteria are *cells without a
nucleus* (prokaryotes) — very different from human cells, so a drug can attack something
bacteria have and we don't. But **fungi are eukaryotes**, just like humans. Their cells are
built on the same basic plan as ours. That means a chemical that poisons a fungus has a
much higher chance of also poisoning the patient. Finding something that kills the fungus
**without killing the person** is genuinely hard — this is called achieving *selectivity*,
and it's the central challenge of antifungal drug design.

---

## 2. *Candida auris*: the pathogen this project targets

*Candida* is a genus (a group) of yeasts. One species, *Candida albicans*, is the common
one behind ordinary yeast infections. A newer species, ***Candida auris*** (often written
*C. auris*), is the villain of this project.

*C. auris* was first identified in 2009 (in the ear canal of a patient in Japan — *auris*
is Latin for "ear"). Since then it has spread through hospitals around the world and earned
a fearsome reputation:

- **The World Health Organization lists it as a "critical priority" fungal pathogen** — the
  top tier of concern.
- It is often **multidrug-resistant**: several of our few antifungal drugs simply don't
  work against many strains.
- It survives on surfaces (bed rails, equipment) for weeks and resists common disinfectants,
  so it spreads person-to-person in hospitals and causes outbreaks.
- Bloodstream infections with it have high mortality.

> **Naming note.** *C. auris* has recently been formally reclassified and renamed
> ***Candidozyma auris***. You'll see both names. The NCBI database this project reads from
> (track 2) uses the folder name `Candidozyma_auris`. They are the same organism.

So: a lethal, spreading, drug-resistant hospital fungus that we're running out of drugs to
treat. That's the problem OpenAFR is aimed at.

---

## 3. How the main antifungal drugs work: the azoles

The most widely used class of antifungal drugs is the **azoles** — names like
**fluconazole**, **itraconazole**, **posaconazole**, and **voriconazole**. (The "-azole"
ending is a hint about their chemistry: they all contain a small nitrogen-rich ring called
an azole ring. Hold onto that nitrogen detail — it becomes central.)

To understand how azoles work, you need to know what they attack.

### The target: an enzyme called CYP51 (a.k.a. ERG11)

Fungal cells are wrapped in a membrane, and to build that membrane they need a molecule
called **ergosterol** (the fungal equivalent of the cholesterol in our own cell membranes).
No ergosterol → no working membrane → the fungus can't grow or survive.

Ergosterol is made in a multi-step assembly line. One essential step is performed by an
enzyme called **CYP51**. (An *enzyme* is a protein that speeds up a specific chemical
reaction.) This same enzyme is also called **ERG11** — that's the name of the *gene* that
carries the instructions for building it. So:

- **ERG11** = the gene (the DNA blueprint).
- **CYP51** = the protein/enzyme that the gene builds.

Throughout this project the two names are used almost interchangeably. When we talk about
DNA and mutations we tend to say ERG11; when we talk about the 3D protein structure we tend
to say CYP51. **Same thing.**

### The mechanism: jamming the enzyme at its iron core

CYP51 belongs to a huge family of enzymes called **cytochrome P450s**. What makes them tick
is a special chemical group buried at the enzyme's center called a **heme** — the very same
kind of iron-containing group that makes your blood red. At the middle of the heme sits a
single **iron atom (Fe)**. That iron is where the enzyme's chemistry happens.

Here's the trick azoles use:

> An azole drug slides into the enzyme's pocket and reaches a **nitrogen atom** from its
> azole ring right up to the **heme iron**, forming a direct bond to it. That nitrogen
> "coordinates" the iron — it plugs the socket. With the iron plugged, CYP51 can't do its
> job. Ergosterol production stalls. The fungus can't build its membrane.

**This single fact is the intellectual heart of the entire project.** The way an azole kills
a fungus is *geometric and physical*: a nitrogen atom getting close enough to the iron atom
to bond with it. Remember this — the whole drug-discovery pipeline (track 1) is built around
measuring exactly that distance.

---

## 4. How resistance happens

A drug that jams an enzyme works only as long as it *fits* the enzyme. Resistance is,
mostly, the fungus changing the fit. There are a few ways it does this, but the one this
project focuses on is **mutations in the ERG11 gene**.

### Mutations and how to read them

DNA is a sequence of chemical "letters." A **mutation** is a change in that sequence. A
protein is a chain of building blocks called **amino acids**, and the DNA sequence dictates
which amino acid goes where. Change the DNA and you can change an amino acid in the protein —
which can subtly reshape the protein's 3D structure.

Scientists write these amino-acid changes in a compact shorthand you'll see constantly in
this project, like **Y132F**:

```
 Y      132      F
 │       │       │
 │       │       └─ the new amino acid now at that spot   (F = Phenylalanine)
 │       └───────── the position number in the protein chain (residue 132)
 └───────────────── the original ("wild-type") amino acid  (Y = Tyrosine)
```

So **Y132F** means: *"at position 132, the amino acid that used to be Tyrosine (Y) is now
Phenylalanine (F)."* Each amino acid has a standard one-letter code. A position in the chain
is also called a **residue**.

> The natural, unmutated version of a gene or protein is called **wild-type**. A mutation is
> described relative to wild-type.

### The resistance mutations this project watches

For CYP51/ERG11 in *C. auris*, a handful of specific mutations are known to cause azole
resistance. This project tracks four of them (its "resistance panel"):

| Mutation | Meaning | Note |
|----------|---------|------|
| **Y132F** | Tyrosine→Phenylalanine at 132 | The single most important, most common global resistance mutation |
| **K143R** | Lysine→Arginine at 143 | Common |
| **F126L** | Phenylalanine→Leucine at 126 | Seen in certain lineages |
| **V125A** | Valine→Alanine at 125 | Part of a specific lineage's signature |

These mutations sit **in or near the pocket where the azole binds**. By slightly changing the
shape or chemistry of that pocket, they make the drug fit worse — so it no longer jams the
iron effectively, and the fungus survives the drug. Notice **Y132F** in the table: it's the
same mutation used earlier as an example, and it comes up again and again in this project
because it's the dominant real-world resistance mutation.

> **Why these exact names matter for the code.** Track 2 of this project reads fungal genome
> data and tries to *call* (detect) mutations like these. Track 1 tries to design drugs that
> still work *despite* them. Both need to name residues precisely — an off-by-one in the
> numbering would turn a real "Y132F" into nonsense, which is why the code obsessively
> verifies that residue 132 really is a Tyrosine in its reference before it trusts any call.

### Other resistance routes (mentioned so you recognize them)

ERG11 point mutations aren't the only way *C. auris* resists azoles. You'll see these terms
in the code and docs, flagged as **out of scope** for now:

- **Efflux pumps** — the fungus pumps the drug back out before it can act (genes like
  *TAC1B* turning these up). Not a change to CYP51 itself.
- **TR34 / tandem repeats ("TR-type")** — extra copies of a control region that make the
  fungus produce *more* CYP51, so a fixed dose of drug can't jam all of it. Written like
  `TR34`. This is **not** a simple amino-acid substitution, which is why the code explicitly
  refuses to treat `TR34` as one (see track 2's mapping stage).

The project's scope is deliberately narrow: **azole resistance via ERG11 point mutations.**
That focus is a feature, not a gap — see the [Project Overview](02-project-overview.md).

---

## 5. Two ways to fight back — and this project does both

Given a lethal, azole-resistant fungus, there are two broad strategies:

1. **Find better drugs.** Design new molecules that jam CYP51 *even in resistant strains*
   where the pocket has been reshaped by mutation. → This is **Track 1** of OpenAFR
   ([the drug-discovery pipeline](03-drug-discovery-pipeline.md)).

2. **Watch the enemy evolve.** Monitor fungal genomes coming out of hospitals worldwide,
   catch new resistance mutations *as they emerge*, and warn people early — with an
   assessment of whether our existing drugs will still work. → This is **Track 2**
   ([the early-warning pipeline](04-early-warning-pipeline.md)).

Track 1 makes new weapons. Track 2 is the radar. They share the same enemy (*C. auris*),
the same enzyme (CYP51/ERG11), and — importantly — the same 3D structural understanding of
how azoles fit the pocket.

---

## 6. A few structural-biology terms you'll meet in the code

The drug-discovery pipeline works with the 3D structure of the CYP51 protein. Here are the
terms it uses, in plain language:

- **Protein structure / PDB file.** Scientists determine the 3D shape of proteins and store
  the atomic coordinates in files. A common format is **PDB** (from the Protein Data Bank, a
  public repository). This project uses a real experimentally-determined structure with the
  ID **5TZ1** — CYP51 from *Candida albicans* (a close relative of *C. auris*), captured with
  a heme and an azole-like molecule bound inside it.

  > *Wait, why the albicans structure and not auris?* Because at the time this was built,
  > **no experimental 3D structure of *C. auris* CYP51 existed.** The two enzymes are very
  > similar, so the project uses the albicans structure as a stand-in — and is upfront that
  > this is a limitation. Track 1 can also *graft* auris mutations onto this scaffold to
  > model resistant versions.

- **Docking.** A computer simulation that takes a small drug molecule (a "ligand") and a
  target protein (a "receptor") and tries to find how the drug best fits into a pocket on the
  protein — like a computer trying thousands of ways to fit a key into a lock. The software
  used here is **AutoDock Vina** (usually just "Vina").

- **Ligand vs. receptor.** The **ligand** is the small molecule (the candidate drug). The
  **receptor** is the big target protein (CYP51). Docking places ligand into receptor.

- **Pose.** One specific predicted arrangement of the ligand inside the pocket — a proposed
  position and orientation. Docking produces several poses per molecule.

- **Docking score / affinity.** Vina outputs a number (in kcal/mol, more negative = it
  *thinks* the fit is stronger) estimating how well each pose binds. **A crucial theme of
  this project is that this score is not very trustworthy** — see track 1.

- **Heme iron / N–Fe distance.** As you learned in section 3, azoles work by getting a
  nitrogen (N) close to the heme iron (Fe). The project measures the **N–Fe distance** — how
  many ångströms (Å, a tiny unit of length, about the width of an atom) separate the closest
  ligand nitrogen from the iron. A short N–Fe distance means the molecule is doing the thing
  a real azole does. This measured distance is track 1's secret weapon.

---

## 7. Recap — the vocabulary in one place

- **Fungal infections** kill ~1.5M/year; hard to treat because fungi are eukaryotes, like us.
- ***C. auris*** = a lethal, spreading, drug-resistant hospital yeast; WHO critical priority.
- **Azoles** (fluconazole, etc.) = the main antifungal drug class.
- **CYP51 = ERG11** = the enzyme azoles jam (the protein / the gene, same thing). It makes
  **ergosterol**, essential for the fungal membrane.
- **Heme iron (Fe)** = the reactive iron atom at CYP51's core. Azoles work by reaching a
  **nitrogen** to it — a geometric event.
- **Resistance** = mutations like **Y132F, K143R, F126L, V125A** reshape the pocket so the
  drug fits worse.
- **Mutation shorthand** = wild-type letter + position + new letter (e.g. Y132F). A position
  is a **residue**. Unmutated = **wild-type**.
- **Docking / ligand / receptor / pose / score** = the computational tools track 1 uses.
- **N–Fe distance** = the geometric measurement track 1 trusts more than the docking score.

You now know enough biology to read the rest of the wiki. Next:
**[Project Overview →](02-project-overview.md)**
