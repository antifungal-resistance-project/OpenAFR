"""Tests for the synthesizability / availability profiler (openafr/synth.py, #86).

What these lock:
  * SAscore is on the Ertl 1..10 scale and its band (easy/moderate/hard) tracks the
    thresholds, with a trivial molecule easy and a fused/complex one harder,
  * availability is honest: 'catalogued-drug' only when told so, 'not-checked' by default
    (no network), and the PubChem presence probe is a WEAK proxy parsed offline via an
    injected opener — a hit -> 'pubchem-record', a miss/None -> 'no-pubchem-record',
  * a network failure in the probe never raises (returns None -> 'no-pubchem-record'),
  * an unparseable SMILES is refused (ValueError), never scored as empty.
"""
import io
import json
import urllib.error

import pytest

from openafr import synth

FLUCONAZOLE = "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F"
ETHANOL = "CCO"
# a strained, polycyclic, stereo-dense natural-product-like fragment -> higher SAscore
COMPLEX = "O=C1[C@@H]2CC[C@H]3[C@@H](CC[C@]4(C)[C@H]3CC[C@]4(O)C#C)[C@@]2(C)CC1"


def fake_opener(payload):
    """Return an opener that yields `payload` (dict -> JSON, or an exception to raise)."""
    def _open(url, timeout=None):
        if isinstance(payload, Exception):
            raise payload
        return io.BytesIO(json.dumps(payload).encode("utf-8"))
    return _open


def test_sascore_scale_and_bands():
    easy = synth.profile(ETHANOL)
    assert 1.0 <= easy["sa_score"] <= 10.0
    assert easy["sa_band"] == "easy"
    hard = synth.profile(COMPLEX)
    assert hard["sa_score"] > easy["sa_score"]


def test_band_thresholds():
    assert synth.sa_band(synth.SA_EASY - 0.01) == "easy"
    assert synth.sa_band(synth.SA_EASY) == "moderate"
    assert synth.sa_band(synth.SA_HARD) == "hard"


def test_availability_is_honest_by_default():
    p = synth.profile(FLUCONAZOLE, name="fluconazole")
    assert p["availability"] == "not-checked"   # no network unless asked
    assert p["pubchem_cid"] is None


def test_catalogued_flag():
    p = synth.profile(FLUCONAZOLE, name="fluconazole", catalogued=True)
    assert p["availability"] == "catalogued-drug"


def test_pubchem_probe_hit_offline():
    opener = fake_opener({"IdentifierList": {"CID": [3365]}})
    p = synth.profile(FLUCONAZOLE, name="fluconazole", opener=opener)
    assert p["availability"] == "pubchem-record"
    assert p["pubchem_cid"] == 3365


def test_pubchem_probe_miss_and_network_error_never_raise():
    # a 404-style miss
    p = synth.profile(FLUCONAZOLE, name="zzznope", opener=fake_opener({"Fault": {}}))
    assert p["availability"] == "no-pubchem-record"
    # a network failure must be swallowed, not raised
    p2 = synth.profile(FLUCONAZOLE, name="fluconazole",
                       opener=fake_opener(urllib.error.URLError("down")))
    assert p2["availability"] == "no-pubchem-record"


def test_unparseable_smiles_is_refused():
    with pytest.raises(ValueError):
        synth.profile("not_a_molecule)))")
