#!/bin/bash
#SBATCH --partition=torino
#SBATCH --account=rubin:commissioning
#SBATCH --job-name=mdet_cell_layout_tests
#SBATCH --output=logs/%x-%j.txt
#SBATCH --error=logs/%x-%j.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --mem=240g
#SBATCH --time=03:00:00

### Layout tests with test_layout.py (galaxies within 20 px of the cell centre,
###    noise-free, 480 galaxies), one per call:
###       spacing_scan  grids of 6 to 15 arcsec, cells of 250 and 500 px
###                     -> layout_response_spacing_scan.npz/.csv
###       fine_scan     grids of 8.5 to 9.5 arcsec in steps of 0.1, cells of 250 and 500 px
###                     -> layout_response_fine_scan.npz/.csv
###       jitter        9-arcsec grid with neighbours jittered by up to +-0.5, 1, 2, 5, 10 px,
###                     and with unsheared positions, cells of 250 and 500 px
###                     -> layout_response_jitter.npz/.csv
###       cell_scan     9-arcsec grid in cells of 240 to 500 px
###                     -> layout_response_cell_scan.npz/.csv
###       map           grids of 43 to 47 px (8.6 to 9.4 arcsec) in steps of 0.5 px, in cells
###                     of 244 to 256 px in steps of 2 px, and of 500 px for reference
###                     -> layout_response_map.npz/.csv
###       fft_margin    cells of 498 to 512 px (GalSim FFT 1024 px, so 1024 - 2 x cell = 28 to 0 px)
###                     and 250 px, grids of 43.5, 44, 45, 88 and 90 px; compared with 500 px
###                     -> layout_response_fft_margin.npz/.csv
###       map500        as map, with cells of 500 to 512 px (and 498 px as the reference)
###                     -> layout_response_map500.npz/.csv
###    Every test includes grid9, which repeats the scenes of run_layout.sh.
###    Plots: python plot_layout_scans.py
###
###    Usage: sbatch --job-name=mdet_cell_<test> run_layout_tests.sh <test>
###    To summarise the saved results again, with the same environment:
###       ./run_layout_tests.sh <test> --summary-only

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"

TEST=${1:?usage: run_layout_tests.sh spacing_scan|fine_scan|jitter|cell_scan|map|fft_margin|map500 [--summary-only]}
shift
export CELL_SIZES=250,500
case "${TEST}" in
    spacing_scan)
        export LAYOUTS=grid6,grid7,grid7.5,grid8,grid8.5,grid9,grid9.5,grid10,grid10.5,grid11,grid12,grid13.5,grid15 ;;
    fine_scan)
        export LAYOUTS=grid8.5,grid8.6,grid8.7,grid8.8,grid8.9,grid9,grid9.1,grid9.2,grid9.3,grid9.4,grid9.5 ;;
    jitter)
        export LAYOUTS=grid9,jitter9_0.5,jitter9_1,jitter9_2,jitter9_5,jitter9_10,grid9_unsheared ;;
    cell_scan)
        export LAYOUTS=grid9
        export CELL_SIZES=240,250,260,270,300,350,400,450,500 ;;
    map)
        export LAYOUTS=grid8.6,grid8.7,grid8.8,grid8.9,grid9,grid9.1,grid9.2,grid9.3,grid9.4
        export CELL_SIZES=244,246,248,250,252,254,256,500 ;;
    fft_margin)
        export LAYOUTS=grid8.7,grid8.8,grid9,grid17.6,grid18
        export CELL_SIZES=250,498,500,502,504,508,512
        export REF_CELL=500 ;;
    map500)
        export LAYOUTS=grid8.6,grid8.7,grid8.8,grid8.9,grid9,grid9.1,grid9.2,grid9.3,grid9.4
        export CELL_SIZES=498,500,502,504,506,508,510,512
        export REF_CELL=498 ;;
    *)
        echo "unknown test ${TEST}" >&2; exit 1 ;;
esac
export LAYOUT_TAG=_${TEST}
export N_TARGETS=480

echo "host      : $(hostname)"
echo "started   : $(date)"
echo "test      : ${TEST}"
echo "layouts   : ${LAYOUTS}"
echo "cells     : ${CELL_SIZES}"
python -u test_layout.py "$@"
echo "finished  : $(date)"
