#!/bin/bash
#SBATCH --partition=torino
#SBATCH --account=rubin:commissioning
#SBATCH --job-name=mdet_cell_layout
#SBATCH --output=logs/%x-%j.txt
#SBATCH --error=logs/%x-%j.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --mem=240g
#SBATCH --time=02:00:00

### test_layout.py: m of galaxies near the cell centre for grids of 9, 15 and 18
###    arcsec, random positions at two densities, and no neighbours, in cells of
###    250 and 500 px. About 100 core-seconds per galaxy, 480 galaxies.
###    Saves layout_response.npz and layout_response.csv.
###
###    Usage: sbatch run_layout.sh
###    To summarise the saved results again: python test_layout.py --summary-only

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

echo "host      : $(hostname)"
echo "started   : $(date)"
N_TARGETS=480 python -u test_layout.py
echo "finished  : $(date)"
