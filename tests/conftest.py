"""Make the scripts/ directory importable so tests can reach the gate code.

The gate currently lives at scripts/validate_gate2.py (not yet a package).
When eng-review task T1 moves metrics()/report() into an openafr/ package,
change the single import in the test module; this path shim can then go away.
"""
import pathlib
import sys

_SCRIPTS = pathlib.Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(_SCRIPTS))
