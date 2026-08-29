# Results — heme/iron model audit (issue #72)

Date: 2026-08-29. Reproduce with `python scripts/audit_heme_model.py`
(deterministic, no network, no docking engine — it reads only the checked-in
crystallographic coordinates `work/receptor_A.pdb` and `work/ref_VT1.pdb`).
Numbers are locked by `tests/test_heme_audit.py`.

> **Why this exists.** The entire OpenAFR criterion is one measured number: the
> closest approach of a ligand nitrogen to the heme iron
> (`openafr.pdbqt.min_nitrogen_iron_distance`). If the iron's coordination
> environment is modeled loosely, every distance that ranks every candidate is
> systematically off. That measurement had never been audited head-on. This note
> establishes that the heme model is physically defensible and that the criterion
> reproduces known azole coordination to within the canonical band.

## 1. The heme parameterization we actually dock into

Receptor `work/receptor_A.pdb` = 5TZ1 (*C. albicans* CYP51) chain A + its heme,
prepared by `scripts/prep_receptor.py` (chain A protein + chain-A `HEM`; waters
and the bound ligand dropped; `obabel -xr -p 7.4` for AutoDock atom types). The
iron's coordination sphere, measured directly from the crystal coordinates:

| coordination position | atom | distance to Fe | reference band |
|---|---|---:|---|
| equatorial | HEM ND | 1.869 Å | 1.8–2.2 Å (porphyrin N) |
| equatorial | HEM NA | 1.983 Å | 1.8–2.2 Å |
| equatorial | HEM NC | 2.055 Å | 1.8–2.2 Å |
| equatorial | HEM NB | 2.096 Å | 1.8–2.2 Å |
| **proximal axial (5th)** | **Cys470 SG (thiolate)** | **2.347 Å** | 2.2–2.5 Å |
| **distal axial (6th)** | **— open —** | — | this is where the azole binds |

This is the canonical **cysteine-thiolate–ligated (CYP450) heme**: four
equatorial porphyrin nitrogens plus the conserved proximal Cys470 thiolate as
the fifth ligand, at the expected ~2.35 Å Fe–S distance. The **sixth (distal)
axial position is deliberately open** — the model carries no explicit water on
the iron (0 `HOH` in the receptor; the prep step drops waters). That is the
correct state for docking: a resting ferric CYP has a distal aqua ligand, but the
whole point of an azole is that its N1 *displaces* the distal water to coordinate
the iron directly, so leaving the sixth site open is what lets the incoming
ligand reach the metal. Modeling an explicit distal water would occlude the very
site the criterion measures into.

## 2. The criterion reproduces known azole coordination

`work/ref_VT1.pdb` is the crystallographic VT-1161 (oteseconazole) pose from
5TZ1. Run through the **production measurement code**
(`min_nitrogen_iron_distance`) against the receptor iron:

- **crystal N–Fe = 2.254 Å**, inside the canonical 2.0–2.3 Å azole-coordination
  band. (Matches `work/RESULTS_redock_VT1.md`, which independently cites 2.25 Å
  for the crystal reference.)

So the number the whole product rests on lands on the physically correct
coordination distance when handed a real, experimentally determined azole–CYP51
complex. This is the direct "does the criterion reproduce known coordination?"
check #72 asked for, and it passes.

## 3. The min-nitrogen heuristic latches the *coordinating* nitrogen

`min_nitrogen_iron_distance` takes the closest of **any** ligand nitrogen, not a
hand-annotated coordinating atom. A fair audit has to check the heuristic isn't
silently measuring a non-coordinating amine. In the VT-1161 crystal:

| ligand N | distance to Fe |
|---|---:|
| **NAF (coordinating azole N1)** | **2.254 Å** |
| NAG | 3.109 Å |
| NAH | 4.155 Å |
| NAD | 4.280 Å |
| NBK | 8.368 Å |

The coordinating triazole nitrogen wins by a **0.855 Å margin** over the next
nitrogen. The heuristic is not ambiguous here: the atom it selects is the atom
that actually forms the coordination bond.

## 4. Charge / oxidation / spin state — settled, not re-run

#72 flags the iron charge/spin state as a possible source of systematic error.
It is not one for this pipeline, and the reason is already on record:
`work/RESULTS_redock_VT1.md` finding 3 showed that **Vina's scoring ignores
receptor partial charges** — patching the heme iron to a formal +3.000 produced
byte-identical poses and scores. AutoDock Vina scores on atom *types* and
geometry, not electrostatics, so the docked pose — and therefore the measured
N–Fe geometry — is invariant to whatever iron charge/oxidation/spin state we
nominally assign. `prep_receptor.py` accordingly writes obabel's all-zero charges
and documents this (lines 29–30). The audit measures geometry precisely because
geometry, not charge, is what the criterion and the engine both act on.

## 5. Sensitivity of the ranking to heme-model variation

Two model choices could in principle move the ranking; both are bounded:

- **Iron charge/oxidation/spin** → *no effect on ranking* (§4: Vina is charge-blind;
  the poses are byte-identical, so the ranking is identical).
- **Distal water present vs. open** → the model keeps the site **open** (§1). An
  explicit distal aqua ligand would sterically block the coordination site and
  push every ligand's N–Fe outward — but that is the wrong physical state for a
  ligand-displacement binder, so it is excluded by construction rather than
  swept over. This is a deliberate, documented choice, not an untested knob.

The remaining foundational geometry (porphyrin N₄ + proximal thiolate) is fixed
by the experimental crystal structure and is not a free parameter.

## Verdict

All five audit checks pass (`scripts/audit_heme_model.py`, exit 0):
four porphyrin N in band, proximal thiolate in band, distal site open, crystal
reproduces coordination (2.254 Å ∈ 2.0–2.3 Å), heuristic latches the coordinating
N (0.855 Å margin). **The heme/iron model the criterion rests on is defensible,
and the criterion reproduces known azole coordination.** The one honest caveat is
that ground-truth reproduction is anchored on a single high-quality co-crystal
(VT-1161/5TZ1, the same complex the redock validation and the docking box derive
from); adding a second independent azole–CYP51 co-crystal would widen this from
n=1 and is the natural follow-up.

**Source:** issue #72. Depends on nothing beyond the checked-in receptor and
crystal-ligand PDBs; reproduces offline via `scripts/audit_heme_model.py`.
