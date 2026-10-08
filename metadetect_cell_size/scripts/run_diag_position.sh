#!/bin/bash
#SBATCH --partition=torino
#SBATCH --account=rubin:commissioning
#SBATCH --job-name=mdet_cell_diag_position
#SBATCH --output=logs/%x-%j.txt
#SBATCH --error=logs/%x-%j.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --time=02:00:00

### Where in the cell does the shear bias come from? (see diag_position_in_cell.py)
###    m1, m2, R11, R22 in bins of distance from the cell centre, and maps over the
###    cell, for every set-up, the baseline 250/150 included.
###    Reads the catalogues only, nothing is written back to them.
###    Saves diag_position_in_cell.csv and diag_position_in_cell_maps.npz, and the
###    plots diag_position_radial.png and diag_position_maps_zero_shear.png.
###    The log ends with a check that the bin 'all' reproduces results_cell_size.csv.
###
###    Usage: sbatch run_diag_position.sh
###    To redo the plots without reading the catalogues again:
###       python diag_position_in_cell.py --plot-only

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

echo "host      : $(hostname)"
echo "started   : $(date)"
python -u diag_position_in_cell.py
echo "finished  : $(date)"
