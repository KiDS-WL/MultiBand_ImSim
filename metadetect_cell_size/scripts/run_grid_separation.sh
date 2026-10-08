#!/bin/bash
#SBATCH --partition=torino
#SBATCH --account=rubin:commissioning
#SBATCH --job-name=mdet_cell_grid_separation
#SBATCH --output=logs/%x-%j.txt
#SBATCH --error=logs/%x-%j.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --time=02:00:00

### m with every second grid galaxy (18 arcsec apart), see test_grid_separation.py
###    Reads the catalogues only, nothing is written back to them.
###    Saves results_grid_separation.csv and m_grid_separation.png.
###    The log checks that the selection 'all' reproduces results_cell_size.csv.
###
###    Usage: sbatch run_grid_separation.sh
###    To redo the plot without reading the catalogues again:
###       python test_grid_separation.py --plot-only

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

echo "host      : $(hostname)"
echo "started   : $(date)"
python -u test_grid_separation.py
echo "finished  : $(date)"
