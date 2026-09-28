"""Tests for the #137 v3 accuracy-fixture DERIVATION (scripts/harvest_diagnostic_accuracy_v3.py).

v3 re-tiers the FKS1 panel (high-PPV core detects; low-PPV D642Y/M690I abstain) and re-derives the
accuracy fixture from the frozen v2 caller output -- no reads, deterministic. These lock the
row-level re-derivation logic and the invariant that a verdict may flip ONLY from
RESISTANCE_MARKER_DETECTED -> UNCHARACTERIZED_VARIANT, and only for a low-PPV-only call.
"""
import importlib.util
import pathlib

from openafr import verdict as V

_SCRIPT = (pathlib.Path(__file__).resolve().parent.parent
           / "scripts" / "harvest_diagnostic_accuracy_v3.py")


def _load():
    spec = importlib.util.spec_from_file_location("harvest_diagnostic_accuracy_v3", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


H = _load()


def _row(called_verdict, panel_hits, susceptibility="R"):
    note = f"harvested=concordance-fks1;resolved=True;panel_hits={panel_hits}"
    return {"isolate": "Bxxxxx", "species": "Candida auris", "run_acc": "SRRx",
            "gene": "FKS1", "drug_class": "echinocandin", "called_verdict": called_verdict,
            "susceptibility": susceptibility, "note": note, "citation": "PMC12323592"}


def test_panel_hits_parsing():
    assert H._panel_hits_of("...;panel_hits=-") == []
    assert H._panel_hits_of("...;panel_hits=S639F") == ["S639F"]
    assert H._panel_hits_of("...;panel_hits=S639F,D642Y") == ["S639F", "D642Y"]
    assert H._panel_hits_of("no marker column") == []


def test_core_hit_stays_detected():
    v3v, tier = H._v3_verdict(_row(V.RESISTANCE_MARKER_DETECTED, "S639F"))
    assert v3v == V.RESISTANCE_MARKER_DETECTED
    assert tier == "core"


def test_low_ppv_only_hit_flips_to_abstain():
    for mut in ("D642Y", "M690I"):
        v3v, tier = H._v3_verdict(_row(V.RESISTANCE_MARKER_DETECTED, mut))
        assert v3v == V.UNCHARACTERIZED_VARIANT, mut
        assert tier == "low_ppv_abstain", mut


def test_core_plus_low_ppv_stays_detected():
    v3v, tier = H._v3_verdict(_row(V.RESISTANCE_MARKER_DETECTED, "S639F,D642Y"))
    assert v3v == V.RESISTANCE_MARKER_DETECTED
    assert tier == "core"


def test_panel_negative_verdict_carried_verbatim():
    for v2v in (V.NO_KNOWN_MARKER, V.UNCHARACTERIZED_VARIANT, V.UNRESOLVED):
        v3v, tier = H._v3_verdict(_row(v2v, "-"))
        assert v3v == v2v
        assert tier == "-"


def test_frozen_v2_fixture_flips_exactly_nine_low_ppv_calls():
    # Data-backed: the derivation over the real frozen v2 fixture flips exactly the 7 D642Y + 2
    # M690I isolates, all detected -> abstain, and nothing else.
    rows = H._read_v2()
    flips = []
    for r in rows:
        v3v, _ = H._v3_verdict(r)
        if v3v != r["called_verdict"].strip():
            flips.append((v3v, r["called_verdict"].strip(), H._panel_hits_of(r.get("note", ""))))
    assert len(flips) == 9
    for v3v, v2v, hits in flips:
        assert v2v == V.RESISTANCE_MARKER_DETECTED and v3v == V.UNCHARACTERIZED_VARIANT
        assert hits and all(h in ("D642Y", "M690I") for h in hits)
