"""Make the repository root importable so tests can reach the openafr package.

metrics()/report() now live in openafr/scoring.py (eng-review task T1), so the
tests import them straight from the package. This shim only ensures the repo root
is on sys.path when pytest is invoked from outside the project root.
"""
import pathlib
import sys

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))
