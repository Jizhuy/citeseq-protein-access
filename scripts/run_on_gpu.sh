#!/usr/bin/env bash
# Activate the dedicated environment and prefer the WSL NVIDIA libraries.
# Usage (from the project root on 4060-server):
#   bash scripts/run_on_gpu.sh python scripts/server_env_check.py
set -euo pipefail

export PATH="/usr/lib/wsl/lib:${PATH}"
export LD_LIBRARY_PATH="/usr/lib/wsl/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

if [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
  # shellcheck source=/dev/null
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "${HOME}/mambaforge/etc/profile.d/conda.sh" ]]; then
  # shellcheck source=/dev/null
  source "${HOME}/mambaforge/etc/profile.d/conda.sh"
elif [[ -f "${HOME}/micromamba/etc/profile.d/micromamba.sh" ]]; then
  # shellcheck source=/dev/null
  source "${HOME}/micromamba/etc/profile.d/micromamba.sh"
  micromamba activate multiomics_robustness
fi

if command -v conda >/dev/null 2>&1; then
  conda activate multiomics_robustness
fi

exec "$@"
