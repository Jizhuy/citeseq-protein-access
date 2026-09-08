#!/usr/bin/env bash
# PHASE 10 on the RTX 4060 CUDA server.
# Gate: seed 0, PBMC10k→PBMC5k, then PBMC5k→PBMC10k.
# Full grid: seeds 1–4, both directions, only after the gate succeeds.
set -euo pipefail

cd /home/jizhu/research/multiomics_robustness
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate multiomics_robustness

export PATH="/usr/lib/wsl/lib:$PATH"
export LD_LIBRARY_PATH="/usr/lib/wsl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
mkdir -p "$HOME/tmp"
export TMPDIR="$HOME/tmp"

python -c "import torch; assert torch.cuda.is_available(), 'CUDA unavailable'"

echo "=== PHASE 10 gate: PBMC10k → PBMC5k, seed 0 ==="
python scripts/10_run_cross_dataset.py --direction pbmc10k_to_pbmc5k

echo "=== PHASE 10 gate: PBMC5k → PBMC10k, seed 0 ==="
python scripts/10_run_cross_dataset.py --direction pbmc5k_to_pbmc10k

echo "=== PHASE 10 full: seeds 1-4, both directions ==="
python scripts/10_run_cross_dataset.py --full

echo "PHASE 10 complete."
