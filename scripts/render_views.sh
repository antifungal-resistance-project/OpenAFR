#!/bin/bash
# Render pose scenes (*.pml from scripts/make_pose_view.py) to PNG with PyMOL.
#
# Visualization is an OPTIONAL step, deliberately kept out of the pinned analysis env
# (environment.yml). PyMOL lives in a separate scratch env so it never perturbs the
# reproducible docking/scoring toolchain. One-time setup:
#
#   mamba create -y -n pymolview -c conda-forge python=3.9 pymol-open-source
#
# Usage: scripts/render_views.sh work/views/*.pml
#   Each .pml already contains its own `png ...` output path (written by make_pose_view.py).
set -eu
ENV="${PYMOL_ENV:-pymolview}"
for pml in "$@"; do
  echo "rendering $pml"
  mamba run -n "$ENV" pymol -cq "$pml"
done
