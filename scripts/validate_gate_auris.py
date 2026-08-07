"""C. auris Y132F resistance-pocket gate — implements work/PREREGISTRATION_auris_Y132F.md.

Same grading engine, pass bar and rank-unrankable-actives-last rule as the run-2 gate
(scripts/validate_gate2.py); the only differences are the mutant receptor and this run's
own frozen pre-registration. See that pre-registration for the pre-committed interpretation
of BOTH a pass (criterion transfers to the resistant pocket) and a fail (it does not —
a publishable negative).

Usage:
    python scripts/validate_gate_auris.py [SCREEN_DIR]
        SCREEN_DIR  docked poses of the held-out set in the Y132F pocket
                    (default work/screen_auris_Y132F)

Produce SCREEN_DIR first by docking the run-2 held-out ligands into the mutant receptor:
    python scripts/mutate_receptor.py Y132F -o work/receptor_auris_Y132F.pdb
    python scripts/prep_receptor.py --mutant work/receptor_auris_Y132F.pdb \
        work/receptor_auris_Y132F.pdbqt
    RECEPTOR=work/receptor_auris_Y132F.pdbqt \
        scripts/screen.sh work/run2_ligands work/screen_auris_Y132F
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from scripts.validate_gate2 import run_gate

PREREG = "work/PREREGISTRATION_auris_Y132F.md"
PREREG_SHA = "689de2fda439bee49e40fef1a8571ebd65a40e1069dae523a82c5c6759acaff1"
RECEPTOR = "work/receptor_auris_Y132F.pdb"  # iron read here; identical to WT iron


def main():
    screen = sys.argv[1] if len(sys.argv) > 1 else "work/screen_auris_Y132F"
    return run_gate(
        screen,
        receptor=RECEPTOR,
        prereg=PREREG,
        prereg_sha=PREREG_SHA,
        title="C. AURIS Y132F — PRE-REGISTERED RESISTANCE-POCKET TRANSFER TEST",
    )


if __name__ == "__main__":
    sys.exit(main())
