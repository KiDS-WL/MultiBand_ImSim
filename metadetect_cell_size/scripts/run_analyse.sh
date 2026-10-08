#!/bin/bash
# @Author: lshuns
# @Date:   1969-12-31 16:00:00
# @Last Modified by:   lshuns
# @Last Modified time: 2026-10-06 07:11:35
#SBATCH --partition=torino
#SBATCH --account=rubin:commissioning
#SBATCH --job-name=mdet_cell_analyse
#SBATCH --output=logs/%x-%j.txt
#SBATCH --error=logs/%x-%j.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --time=04:00:00

### R, m and c for every cell set-up, the baseline 250/150 included.
###    Set-ups with no catalogue yet are skipped, incomplete ones are flagged,
###    so it can be run before every job has finished.
###    Reads the catalogues only, nothing is written back to them.
###    Summary saved as results_cell_size.csv
###
###    Usage: sbatch run_analyse.sh

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

echo "host      : $(hostname)"
echo "started   : $(date)"
python -u analyse_cell_size.py
echo "finished  : $(date)"
