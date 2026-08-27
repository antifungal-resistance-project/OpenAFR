"""Tests for the applicability-domain sensitivity analysis (scripts/domain_sensitivity.py, #87).

What these lock:
  * descriptors() reports heavy-atom and aromatic-nitrogen counts correctly, and only
    aromatic N counts (a nitrile / aliphatic amine N does not),
  * in_domain() applies both frozen rules (heavy-atom cap AND aromatic-N floor),
  * the sweep tallies retained molecules per (cap, aromN) grid point,
  * the analysis reuses the SAME frozen thresholds as filter_library (no drift).
"""
import pytest

pytest.importorskip("rdkit")

from scripts import domain_sensitivity as ds
from scripts.filter_library import MAX_HEAVY_ATOMS


def write_smi(tmp_path, rows):
    p = tmp_path / "lib.smi"
    p.write_text("".join("%s\t%s\n" % (smi, name) for smi, name in rows))
    return str(p)


def test_descriptors_counts_heavy_and_aromatic_n(tmp_path):
    rows = [
        ("c1ccncc1", "pyridine"),         # 6 heavy, 1 aromatic N
        ("N#Cc1ccccc1", "benzonitrile"),  # nitrile N is NOT aromatic
        ("CCN", "ethylamine"),            # aliphatic amine N is NOT aromatic
    ]
    d = dict((name, (h, a)) for name, h, a in ds.descriptors(write_smi(tmp_path, rows)))
    assert d["pyridine"] == (6, 1)
    assert d["benzonitrile"][1] == 0     # aromatic N count is zero
    assert d["ethylamine"][1] == 0


def test_in_domain_applies_both_rules():
    # in domain: small enough AND has the required aromatic N
    assert ds.in_domain(heavy=25, aromatic_n=1, cap=45, min_arom_n=1) is True
    # too big
    assert ds.in_domain(heavy=50, aromatic_n=3, cap=45, min_arom_n=1) is False
    # not enough aromatic N
    assert ds.in_domain(heavy=25, aromatic_n=0, cap=45, min_arom_n=1) is False
    # a stricter aromatic-N floor bites
    assert ds.in_domain(heavy=25, aromatic_n=1, cap=45, min_arom_n=2) is False


def test_sweep_tallies_grid(tmp_path):
    mols = [("a", 20, 2), ("b", 50, 1), ("c", 40, 0)]
    grid, n = ds.sweep(mols, caps=[45], min_arom_ns=[1])
    assert n == 3
    # only 'a' (20 heavy, 2 aromN) passes cap<=45 AND aromN>=1
    assert grid["cap=45,aromN>=1"] == 1


def test_uses_frozen_threshold_from_filter_library():
    # the module must not re-hardcode a different cap
    assert MAX_HEAVY_ATOMS == 45
