#!/bin/bash
set -o pipefail
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate multiomics_robustness
cd /home/jizhu/research/multiomics_robustness
export TMPDIR="$HOME/tmp"
export MPLBACKEND=Agg
export PATH="/usr/lib/wsl/lib:$PATH"
export PYTHONPATH="/home/jizhu/research/multiomics_robustness"
LOG=results/revision_round2_methodological_fixes/controls/test_time_corruption/logs/orchestrator.log
echo "RELAUNCH_START=$(date -Is)" | tee -a "$LOG"
python -u results/revision_round2_methodological_fixes/scripts/revision2_testtime_corruption.py --run --gate C --resume 2>&1 | tee -a "$LOG"
echo "GATEC_EXIT=$?" | tee -a "$LOG"
echo "RELAUNCH_END=$(date -Is)" | tee -a "$LOG"
