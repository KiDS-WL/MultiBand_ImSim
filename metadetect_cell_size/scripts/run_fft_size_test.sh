#!/bin/bash
#SBATCH --partition=torino
#SBATCH --account=rubin:commissioning
#SBATCH --job-name=mdet_cell_fft_size
#SBATCH --output=logs/%x-%j.txt
#SBATCH --error=logs/%x-%j.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --mem=240g
#SBATCH --time=03:00:00

### Which GalSim FFT inside metacal sets the bias? test_layout.py on grids of 43.5, 44 and
###    45 px with GalSim's sizes inside metacal changed (test_hard_cut.py MCAL_FFT_SIZE,
###    MCAL_PAD_FACTOR; scenes unchanged):
###       W    size of the drawImage FFT (default: GalSim's choice, 384 px for most draws of a
###            250-px cell, 512 px for a 500-px one)
###       pad  pad_factor of the cell's InterpolatedImage (default 4: a table P of 1024 px for
###            250-px cells, 2048 px for 500-px ones)
###    Variants are name:cells:W:pad (empty: GalSim's default). Each writes
###    layout_response_fft_<name>.npz/.csv. Plots: python plot_layout_scans.py
###
###    Usage: sbatch run_fft_size_test.sh
###           ./run_fft_size_test.sh --summary-only        # summarise the saved results again
###           VARIANTS="W512:250:512:" sbatch run_fft_size_test.sh

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"

VARIANTS=${VARIANTS:-"default:250:: W256:250:256: W320:250:320: W384:250:384: W448:250:448: W512:250:512: \
W640:250:640: W768:250:768: W1024:250:1024: W2048:250:2048: pad6:250::6 pad8:250::8 \
pad8_W512:250:512:8 pad6_W768:250:768:6 cell500_W768:500:768: cell500_W1024:500:1024:"}
export LAYOUTS=grid8.7,grid8.8,grid9
export N_TARGETS=480

echo "host      : $(hostname)"
echo "started   : $(date)"
echo "layouts   : ${LAYOUTS}"
for v in ${VARIANTS}; do
    IFS=: read -r name cells W pad <<< "${v}"
    echo
    echo "=========== ${name}: cells ${cells}, drawImage FFT ${W:-default}, pad_factor ${pad:-default}  ($(date +%T))"
    LAYOUT_TAG=_fft_${name} CELL_SIZES=${cells} MCAL_FFT_SIZE=${W} MCAL_PAD_FACTOR=${pad} \
        python -u test_layout.py "$@"
done
echo "finished  : $(date)"
