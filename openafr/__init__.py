"""OpenAFR — shared library behind two products that share one structural core.

Two products live in this repo (refocus plan, 2026-08-21). They are deliberately
NOT merged — different users, different deliverables — but they share a load-bearing
structural core, and the coupling between them runs one direction only. Read this
docstring to tell the two apart without opening a refactor.

  1. CYP51 drug-triage (the original, validated product)
     ------------------------------------------------------
     A self-validating structural screen: rank azoles against the CYP51/ERG11 pocket
     by heme-iron-approach geometry, gated by a blind-holdout + permutation test that
     must re-discover known drugs before any run is trusted. User: a med-chemist /
     wet-lab PI who wants a credible shortlist. Most of it lives in scripts/
     (validate_gate*.py, screen.sh, mutate_receptor.py); in-package it is the shared
     core below plus the pocket-fit engine.

  2. Genomic surveillance (the walled-off second product)
     ------------------------------------------------------
     A C. auris resistance early-warning track (issues #20-27, + FKS1 v2): take NCBI
     isolates, re-call ERG11/FKS1 resistance substitutions from raw reads, detect
     NEW/RISING mutations, and emit an alert with a structural so-what. User: a
     hospital epidemiologist / public-health analyst — not a med-chemist. This is a
     candidate to spin out; keep it legibly separate.

       earlywarning  snapshot store + baseline (#21)
       emergence     NEW/RISING detection, archival-backlog guard (#22)
       mapping       mutation -> pocket-residue mapping (#23)
       structural    the so-what: azole-fit consequence of a flagged mutation (#24)
       alert         assemble the human/machine alert (#25)
       delivery      deliver / no-op decision (#26)
       backtest      honest validation vs an external truth set (#27)
       recaller      ERG11 reads -> azole-resistance call
       fks1_caller   FKS1 reads -> echinocandin-resistance call (v2)
       runlog        append-only provenance for the un-reproducible re-caller runs

  Shared structural core (load-bearing, one home per eng-review de-duplication)
  -----------------------------------------------------------------------------
     openafr.pdbqt     read docked poses, Vina scores, and the heme iron position
     openafr.scoring   metrics() / report() — the ranking math the whole result rests on
     openafr.protocol  the frozen, hash-verified run protocol (protocol.yaml)

  The coupling, and why it is benign
  ----------------------------------
  Surveillance depends on triage, never the reverse: `structural.py`'s so-what verdict
  (#24) consumes the CYP51 pocket tooling (pdbqt/protocol/geometry) that the drug-triage
  product owns. That is a consumer depending on a core, not a tangle — so the refocus
  plan (D7) documents the boundary here rather than forcing an `openafr/core/` extraction.
  Extract only when the code actually needs it.
"""
