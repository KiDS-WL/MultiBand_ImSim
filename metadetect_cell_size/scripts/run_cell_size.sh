#!/bin/bash
#SBATCH --partition=torino
#SBATCH --account=rubin:commissioning
#SBATCH --job-name=mdet_cell
#SBATCH --output=logs/%x-%j.txt
#SBATCH --error=logs/%x-%j.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --mem-per-cpu=4g
#SBATCH --time=1-00:00:00
#SBATCH --gpus 0

### One shear tag of one cell set-up: task 6_2 (metadetect) only.
###
###    The images are the collaborator's, linked into <out_dir>/<runTag>/images,
###    so no image is simulated here. The resources and rng_seed are those of the
###    collaborator's cell_size = 250, central_size = 150 run, which took
###    1.2-1.4 h per tag; the cost grows with the number of cells (see
###    submit_cell_size.sh).
###
###    A tile whose catalogue already exists is skipped, so a job that ran out
###    of time can simply be resubmitted. Do not run two jobs on the same
###    set-up and tag at once.
###
###    Usage: sbatch run_cell_size.sh <runTag> <config.ini>
###       sbatch run_cell_size.sh p000p000 ./config_cell250_central50.ini
###    Submit everything with submit_cell_size.sh.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

TAG=${1:?need a runTag, e.g. p000p000}
CONFIG=${2:?need a config file, e.g. ./config_cell250_central50.ini}

echo "host      : $(hostname)"
echo "started   : $(date)"
echo "runTag    : ${TAG}"
echo "config    : ${CONFIG}"
grep -E "^(cell_size|central_size) " "${CONFIG}"

## --cosmic_shear is only used to simulate images (task 1), so it is not passed
python ../../modules/Run.py 6_2 --runTag "${TAG}" --threads "${SLURM_CPUS_PER_TASK:-120}" \
    --rng_seed 940120 -c "${CONFIG}" --sep_running_log

echo "finished  : $(date)"
