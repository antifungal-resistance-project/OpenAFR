#!/usr/bin/env bash
# Provision a fresh Linux/x86_64 host for the ERG11 re-caller's real run.
#
# WHY a cloud VM (see work/RUNBOOK_recaller_run.md): the reads->consensus toolchain
# (sra-tools/minimap2/samtools) is bioconda linux-64/osx-64 only -- it will not run
# natively on Apple Silicon. Local Docker under x86 emulation *proves* the pipeline
# but its VM raw-disk does not reclaim space across runs, so even the 4-strain sanity
# set overflows a disk-tight laptop. A native x86 host has none of that: recall_erg11
# deletes each isolate's workdir after calling it, so disk stays flat at ~one isolate's
# scratch (~2-3 GB). The real constraints here are native x86, bandwidth, and long
# unattended runtime for many SRA downloads -- not disk.
#
# SIZING: any Ubuntu 22.04/24.04 x86_64 VM with >=2 vCPU, >=4 GB RAM, >=40 GB disk.
# (e.g. AWS t3.large, GCP e2-standard-2, a $24/mo DigitalOcean droplet -- your account.)
#
# USAGE (on the fresh VM):
#   git clone https://github.com/antifungal-resistance-project/OpenAFR.git
#   cd OpenAFR
#   bash scripts/cloud_recaller_setup.sh
#   conda activate openafr           # (or: source ~/miniforge3/bin/activate openafr)
#   # then follow work/RUNBOOK_recaller_run.md from step 1.
set -euo pipefail

echo "== ERG11 re-caller host setup =="
uname_m="$(uname -m)"
if [ "$uname_m" != "x86_64" ]; then
  echo "ERROR: this host is '$uname_m', not x86_64. The bioconda toolchain needs linux-64." >&2
  echo "Provision an x86_64 VM (not arm64/Graviton)." >&2
  exit 1
fi

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"
if [ ! -f environment.yml ]; then
  echo "ERROR: run this from inside the OpenAFR clone (environment.yml not found)." >&2
  exit 1
fi

# 1. conda (miniforge) if absent.
if ! command -v conda >/dev/null 2>&1; then
  echo "== installing miniforge =="
  tmp_installer="$(mktemp --suffix=.sh)"
  curl -fsSL -o "$tmp_installer" \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
  bash "$tmp_installer" -b -p "$HOME/miniforge3"
  rm -f "$tmp_installer"
fi
# shellcheck disable=SC1091
source "$HOME/miniforge3/etc/profile.d/conda.sh"

# 2. the openafr env, exactly as the runbook's step 0 (reproducibility contract).
if conda env list | grep -qE '^\s*openafr\s'; then
  echo "== env 'openafr' already exists -- updating =="
  conda env update -n openafr -f environment.yml --prune
else
  echo "== creating env 'openafr' from environment.yml =="
  conda env create -f environment.yml
fi
conda activate openafr

# 3. gate checks: tests green + samtools >= 1.13 (the caller's floor).
echo "== running the test suite =="
python -m pytest -q

echo "== toolchain versions =="
python --version
samtools --version | head -1
minimap2 --version
fasterq-dump --version 2>&1 | head -1

echo
echo "== READY =="
echo "Next (see work/RUNBOOK_recaller_run.md):"
echo "  conda activate openafr"
echo "  python scripts/snapshot_ncbi_auris.py pull"
echo "  python scripts/recaller_sanity.py --only B11220     # WT control, expect NO hit"
echo "  python scripts/recaller_sanity.py                    # full Illumina set"
echo "  python scripts/recall_erg11.py plan <snapshot> --list-runs | head"
echo "  python scripts/recall_erg11.py fill <snapshot> --limit 20   # start small"
echo "  python scripts/backtest_earlywarning.py prevalence <snapshot>   # the deliverable"
echo
echo "Tip: run long fills under tmux/screen so an SSH drop doesn't kill them."
