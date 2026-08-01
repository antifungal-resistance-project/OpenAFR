"""OpenAFR — shared library for the self-validating antifungal screening pipeline.

Historically every script under scripts/ re-implemented the same PDBQT/Vina
parsers and the AUC/EF/BEDROC scoring by hand (eng-review finding: four copies of
the pose reader, two copies of the metrics). This package is the single home for
that load-bearing code so a fix or a bug lives in exactly one place.

    openafr.pdbqt     read docked poses, Vina scores, and the heme iron position
    openafr.scoring   metrics() / report() — the ranking math the whole result rests on
    openafr.protocol  the frozen, hash-verified run protocol (protocol.yaml)
"""
