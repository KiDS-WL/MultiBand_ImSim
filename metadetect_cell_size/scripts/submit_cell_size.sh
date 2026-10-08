#!/bin/bash
### Submit the metadetect cell-size test: one job per cell set-up and shear tag.
###
###    first round:
###       cell_size 250 with central_size  50, 100, 200
###       cell_size 500 with central_size 250, 300, 400
###    second round, the same central sizes as cell 250 to separate the two:
###       cell_size 500 with central_size  50, 100, 200
###    third round, a cell size whose GalSim FFT sizes inside metacal are commensurate
###    (image table / drawImage FFT = 4), with the baseline's central size:
###       cell_size 260 with central_size 150
###    (250/150 is the collaborator's run, used as the baseline)
###
###    5 jobs per set-up, one per shear tag. The set-ups are named explicitly, so
###       that finished ones are not resubmitted by accident.
###
###    Usage:
###       ./submit_cell_size.sh cell500_central200 cell500_central100 cell500_central50
###
###    Walltime: the cost grows with the number of cells, ~ 1/central_size^2, and
###       with the cell area. Measured per shear tag with 120 cores:
###          250/50 8.1-9.8 h   250/100 2.3 h   250/200 0.95 h
###          500/250 1.4 h      500/300 1.1 h   500/400 0.85 h
###       i.e. about 0.24 h + 2.6e-4 h per cell at cell 250, and 0.49 h + 6.8e-4 h
###       per cell at cell 500. The second round is then expected to take
###          500/200 ~1.9 h   500/100 ~6 h   500/50 ~22 h
###       so the central size 50 runs are given 2 days and the others 1 day.
###
###    After they finish, measure R, m and c with:
###       sbatch run_analyse.sh

cd "$(dirname "$0")"
mkdir -p logs

TAGS="m020m020 m020p020 p000p000 p020m020 p020p020"
SETUPS="$*"
if [ -z "${SETUPS}" ]; then
    echo "usage: $0 <setup> [<setup> ...]" >&2
    echo "available: $(ls config_cell*_central*.ini | sed 's/^config_//; s/\.ini$//' | tr '\n' ' ')" >&2
    exit 1
fi

## check everything before submitting anything
for setup in ${SETUPS}; do
    config=./config_${setup}.ini
    if [ ! -f "${config}" ]; then
        echo "missing config: ${config}" >&2; exit 1
    fi
    ## without the images link, task 6_2 would create an empty images folder
    out_dir=$(awk '/^out_dir *=/{print $3}' "${config}")
    for tag in ${TAGS}; do
        if [ ! -d "${out_dir}/${tag}/images/original" ]; then
            echo "no images for ${tag}: ${out_dir}/${tag}/images/original" >&2; exit 1
        fi
    done
done

for setup in ${SETUPS}; do
    case "${setup}" in
        *_central50) walltime=2-00:00:00 ;;
        *)           walltime=1-00:00:00 ;;
    esac
    for tag in ${TAGS}; do
        sbatch --job-name="mdet_${setup}_${tag}" --time="${walltime}" \
            run_cell_size.sh "${tag}" "./config_${setup}.ini"
    done
    echo "submitted 5 tags for ${setup} (walltime ${walltime})"
done
