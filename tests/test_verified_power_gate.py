"""Guards for the verified-inactive POWER gate (scripts/validate_gate_verified_power.py).

The gate's ranking/subgroup math is the same code look #6 uses and is covered by
tests/test_verified_gate.py. What is new here and worth locking:

  * the three input hashes the grader refuses to score without must match the real files on
    disk, so a future accidental edit to the frozen actives/inactives/protocol is caught by
    the test suite, not silently graded;
  * the heavy-atom split (S4 domain-edge) reads the count RDKit reports, not a stale cache.
"""
import hashlib

import scripts.validate_gate_verified_power as g


def _sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def test_pinned_input_hashes_match_the_files_on_disk():
    """The grader certifies a pass only if these three match; keep them honest."""
    assert _sha(g.ACTIVES) == g.ACTIVES_SHA
    assert _sha(g.INACTIVES) == g.INACTIVES_SHA


def test_pinned_actives_and_inactives_have_the_expected_sizes():
    """N=300 measured-active azoles, 279 measured inactives — the single-variable design."""
    n_act = sum(1 for l in open(g.ACTIVES) if "\t" in l)
    n_ina = sum(1 for l in open(g.INACTIVES) if "\t" in l)
    assert (n_act, n_ina) == (300, 279)


def test_heavy_atoms_counts_from_structure(tmp_path):
    """Benzene has 6 heavy atoms; the domain-edge split reads that, not a cached number."""
    p = tmp_path / "set.smi"
    p.write_text("c1ccccc1\tbenzene\nCCO\tethanol\nnot-a-smiles\tbad\n")
    got = g._heavy_atoms(str(p))
    assert got == {"benzene": 6, "ethanol": 3}


def test_no_name_collision_between_actives_and_inactives():
    """A shared name would be graded as whichever set is read second — must be disjoint."""
    act = {l.split("\t")[1].strip() for l in open(g.ACTIVES) if "\t" in l}
    ina = {l.split("\t")[1].strip() for l in open(g.INACTIVES) if "\t" in l}
    assert act.isdisjoint(ina)
