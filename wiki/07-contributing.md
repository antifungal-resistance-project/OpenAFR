# Contributing

*Audience: you've read enough of the wiki to know what OpenAFR is, and you want to set up,
run the tests, and make a change. This page gets you from a fresh clone to a first
contribution.*

---

## 1. Set up the environment

The whole toolchain is version-pinned in [`environment.yml`](../environment.yml) (rdkit,
openbabel, AutoDock Vina, pytest, pytest-cov). Pinning is part of the reproducibility contract —
a result you can't rebuild is not a result.

```bash
conda env create -f environment.yml     # or: mamba env create -f environment.yml
conda activate openafr
```

Works on Apple Silicon and Linux. **One important exception:** the ERG11 re-caller's
reads→consensus step needs bioinformatics binaries (`fasterq-dump`, `minimap2`, `samtools`) that
bioconda only ships for Linux/x86 — see [track 2](04-early-warning-pipeline.md) and
`scripts/cloud_recaller_setup.sh`. Everything else, including the entire deterministic core and
all the tests, runs on either platform.

## 2. Run the tests first

Before changing anything, confirm a green baseline:

```bash
pytest
```

Every test is a **known-answer** test — it locks a specific numeric or structural output, so if
you accidentally change behavior, a test fails loudly with the before/after. CI also enforces a
coverage gate on the `openafr/` core (~97%). If you add code to `openafr/`, add a test for it.

## 3. See a pipeline actually run

You don't need NCBI or a wet lab to watch either track work.

**Track 2 (no external anything — synthetic demo modes):**

```bash
python scripts/detect_emergence.py demo
python scripts/structural_sowhat.py estimate Y132F K143R F126L V125A TR34
python scripts/deliver_alert.py demo --markdown
python scripts/backtest_earlywarning.py mechanism
```

**Track 1 (needs the docking toolchain from step 1):** follow the
[README's reproduce recipe](../README.md#reproducing) — prep the receptor, prep the validation
ligands, run the screen, and grade it with `validate_gate2.py`. A `0` exit is a pass.

## 4. The one non-negotiable: the honesty bar

This is the most important section on the page. The entire value of OpenAFR is that its claims
are trustworthy. A contribution that adds a feature but quietly weakens that trust is a net
negative. Before you open a PR, check your change against these — they're the patterns you'll
have seen repeated across the codebase (and collected in the
[overview](02-project-overview.md#the-recurring-design-patterns-learn-these-once)):

- **Never assume the convenient default.** If the pipeline doesn't know something, it must record
  that it doesn't know — an uncalled isolate is not "susceptible," a molecule with no usable pose
  is ranked *last*, not dropped. Don't add a code path that silently fills in the flattering
  answer.
- **Don't move a frozen goalpost.** If you need to change `protocol.yaml`, a `PREREGISTRATION_*.md`,
  or the pinned ERG11 reference, that's a deliberate act: change the file, **re-freeze its hash**
  (e.g. `python -m openafr.protocol --freeze`, then update `PROTOCOL_SHA`), and **record why** in
  the commit and the relevant write-up. Never edit a frozen file to make a run pass.
- **Separate "measured" from "estimated" from "can't say."** If you add a structural or
  predictive claim, tag its confidence explicitly. "No confident call" is a valid, expected
  output — prefer it to a guess.
- **State the limitation next to the result.** Follow the existing docs: when you report a number,
  report what would undermine it in the same place. The README keeping the inflated EF@1% figure
  *with its caveat attached* — rather than deleting it — is the model.
- **Keep the tested core pure.** Logic that must be correct goes in `openafr/` (stdlib or rdkit
  only) with a test. Code that shells out to external tools or the network stays in `scripts/`
  and is honestly labelled as untested-in-CI. Don't drag a heavy dependency (pandas, PyYAML,
  etc.) into the core — the small hand-rolled parsers exist on purpose.
- **Commit the small durable memory, not the big regenerable blob.** Provenance manifests
  (`INDEX.tsv`, run-log `.jsonl`, `STATE.json`) are committed; multi-MB snapshots and receptors
  are gitignored and rebuilt. Match that split.

If your change can't satisfy one of these, that's worth a conversation in the PR rather than a
workaround.

## 5. Where to start — good first contributions

- **This wiki.** If a page here overclaims, is unclear, or has drifted from the code, fix it.
  Accuracy here is a real contribution.
- **Tests.** Find an `openafr/` branch that isn't covered and lock it with a known-answer test.
- **Track 2 demo polish / new demo modes** — low-risk, self-contained, no external deps.
- **The big open items** are in [TODOS.md](../TODOS.md). Current headliners: the
  **coordinator-identity** follow-through on the shortlist (track 1), and a real **FKS1 `fill`**
  to measure the echinocandin event frequency (track 2). The FKS1 (and the earlier ERG11) run
  needs a Linux/x86 host and large downloads, so read the runbooks (`work/RUNBOOK_fks1_run.md`,
  `work/RUNBOOK_recaller_run.md`) and the TODO's run order first, and coordinate before starting.
  (The ERG11 re-caller has already been run — see `work/RESULTS_prevalence.md`.)

Before picking anything up, skim the relevant `work/RESULTS_*.md` write-up — most stages have one,
and it tells you what's been tried, what worked, and what was deliberately deferred.

## 6. Scope discipline

The project's power comes from being **narrow**: CYP51/ERG11 azole resistance in *C. auris*, plus
the one data-justified adjacency it has taken on — **FKS1/echinocandin *detection*** (v2), opened
*because* the measured 80.4% azole baseline showed azole emergence is a weak signal, and kept
honest as detection-only (no structural verdict). Everything past that — TAC1B efflux, ERG3,
diagnostics, other pathogens, or a *structural* FKS1 verdict — remains **out of scope until the
current work is validated** (a wet-lab hit for track 1; a wet-lab-anchored backtest for track 2).
[VISION.md](../VISION.md) explains the "broad mission, focused execution" reasoning. New scope is
a strategy conversation, not a quiet PR.

## 7. License & contributions

OpenAFR is under the **PolyForm Noncommercial License 1.0.0** ([LICENSE](../LICENSE)) — free for
noncommercial use, commercial use requires a separate license from the foundation. Copyright is
held by the foundation so it can offer those commercial licenses; outside contributions are
accepted under the same terms (a contributor license agreement may be requested) to keep that
possible. If that matters for your situation, raise it before contributing.

---

That's the whole loop: set up → test → run a demo → make a change that respects the honesty bar →
open a PR. Welcome aboard.

← Back to the [Wiki home](README.md).
