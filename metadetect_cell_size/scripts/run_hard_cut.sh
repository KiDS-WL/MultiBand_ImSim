#!/bin/bash
#SBATCH --partition=torino
#SBATCH --account=rubin:commissioning
#SBATCH --job-name=mdet_cell_hard_cut
#SBATCH --output=logs/%x-%j.txt
#SBATCH --error=logs/%x-%j.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128g
#SBATCH --time=02:00:00

### Part b of test_hard_cut.py: m of galaxies near the cell centre, for grid stamps
###    of 44 px (9-arcsec grid), 74 px (15-arcsec grid) and 256 px (uncut), in cells
###    of 250 and 500 px, on synthetic noise-free scenes drawn with ImSim's code.
###    About 85 core-seconds per galaxy.
###    Saves hard_cut_response.npz and hard_cut_response.csv.
###
###    Usage: sbatch run_hard_cut.sh
###    To summarise the saved results again: python test_hard_cut.py b --summary-only

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

echo "host      : $(hostname)"
echo "started   : $(date)"
N_TARGETS=480 python -u test_hard_cut.py b
echo "finished  : $(date)"
