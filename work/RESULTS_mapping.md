# Mutation → pocket mapping (issue #23)

**What:** the structural half of the C. auris azole early-warning loop, and where the
moat starts. The emergence detector (#22) turns a snapshot diff into flagged ERG11
substitution tokens — "Y132F is rising", "F126L is new". This stage turns a *token*
into a *structural statement*: where the residue sits in the modeled azole pocket,
whether it touches the drug, how far it is from the catalytic heme iron, and whether
we can build the mutant and dock it. A genomic feed can flag the token; only this step
says what it means for the drug.

**Code:** `openafr/mapping.py` (library) · `scripts/map_mutation.py` (CLI:
`map <tokens>`, `from-emergence`, `--build`) · `tests/test_mapping.py` (11 offline
tests). Reproduce: `python scripts/map_mutation.py map Y132F K143R F126L TR34`.

It reuses the exact receptor the whole project docks against — the 5TZ1 chain-A
*C. albicans* CYP51 scaffold (`work/receptor_A.pdb`), same heme, same iron, same frame
as the frozen protocol and the *C. auris* Y132F port (#13) — so a mapping verdict lives
in the same structural universe as every gate run.

## How the pocket is defined — measured, not hand-picked

A residue **lines the pocket** if any of its heavy atoms comes within **4.5 Å** (a
standard vdW/contact distance) of any atom of the co-crystallized azole in
`work/ref_VT1.pdb` (the VT1 ligand from 5TZ1). This is the textbook definition of a
binding-site residue, computed from the structure at map time — so "Y132 contacts the
drug" is *measured*, never asserted. Each verdict also reports the residue's closest
approach to the heme iron, the atom azoles coordinate.

## The verdicts on the three canonical C. auris ERG11 mutations

```
  Y132F   POCKET-CONTACT drug 3.58 A  iron 8.45 A   [buildable(deletion)]
  F126L   POCKET-CONTACT drug 3.68 A  iron 10.11 A  [needs-repacker]
  K143R   second-shell   drug 7.39 A  iron 8.17 A   [needs-repacker]
```

The structural so-what the metadata feed cannot give you:

- **Y132F and F126L directly line the azole pocket** — 3.58 Å and 3.68 Å from the
  co-crystal drug. Mutating them directly reshapes the drug-contact surface.
- **K143R is second-shell.** At 7.39 Å it does **not** contact the drug, consistent
  with the literature reading that K143R perturbs the pocket more indirectly than the
  two direct contacts. This is exactly the distinction an alert needs to prioritise a
  signal, and it falls straight out of the structure.

## The three honesty constraints (inherited from #22 / mutate_receptor)

1. **We assert the wild type, we do not assume it.** A token's wt letter must match the
   residue actually present at that position. `A132F` is refused — *"position 132 is TYR
   in the model, token expects ALA — numbering or structure mismatch"* — a numbering or
   structure error is surfaced, never mapped through.

2. **Non-residue tokens are refused, not forced.** `TR34` (a promoter / tandem-repeat
   call) is not a point substitution, so it cannot map to one residue and says so. An
   absent residue (`Y9999F`, outside the model's 45–528) is likewise refused. We do not
   pick a residue for a token that does not name one.

3. **Buildability is honest about the pinned environment.** The env has no side-chain
   repacker, so only pure atom-**deletion** mutations are deterministically buildable
   (`Y132F`, `V125A`) — the *same* whitelist `mutate_receptor.py` enforces, locked in
   sync by a test. A mutation that **adds** atoms (`K143R`, `F126L`) is *mapped* (we
   know where it is and whether it touches the drug) but reported as **needs-repacker** —
   never a silently guessed side chain. Mapping a residue and building its mutant side
   chain are two different claims, kept separate.

## The loop closes, and the port is reproduced

`map_mutation.py from-emergence` runs the #22 demo detector and maps its own output:
the four flagged tokens (`F126L, Y132F, K143R, TR34`) come back as three mapped verdicts
(two pocket-contact) and one honest refusal — a complete detect → map loop with no
manual hand-off.

`--build Y132F` reuses the C. auris port from #13 end-to-end: it invokes
`mutate_receptor.py` to graft the mutation onto the 5TZ1 scaffold and emits the exact
`prep_receptor` / `make_pose_view` / `render_views` command chain to view the mutant
pocket. The rebuilt `work/receptor_auris_Y132F.pdb` is **byte-identical** to the #13
artifact (sha256 `924cf7b8…`), heme iron unchanged — the structural path from a flagged
token to a dockable mutant receptor is live and reproducible.

## What is deferred, honestly

- **Add-atom mutations (K143R clade I, F126L clade III) are mapped, not docked.** Their
  mutant side chains need a pinned repacker — the same documented follow-up as the auris
  port (Path A, step 2). The pocket verdict for them is real; the mutant-docking verdict
  waits on repacker tooling.
- **Mapping is against a single modeled chain** (5TZ1 chain A, *C. albicans* numbering).
  C. auris ERG11 mutations are reported in C. albicans-equivalent numbering, which the
  wild-type-assertion check enforces per token; a token in a different numbering scheme
  will trip the mismatch guard rather than map silently.
