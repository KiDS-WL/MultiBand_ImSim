# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-07-03 15:00:52
# @Last Modified by:   lshuns
# @Last Modified time: 2021-07-14 20:21:44

### a simple script to combine catalogues produced by the main pipeline
###     it can be used to combine catalogues from different running_tags
###     it includes a tile_label column for noise info
###     it includes a gal_rot column for rotation info
###     it only preserves specified columns

########################
# usage: combine_SimCata.py [-h] [--in_dir IN_DIR] [--outfile_path OUTFILE_PATH]
#                           [--wanted_cols WANTED_COLS] [--preserve_stars]
#                           [--preserve_false_detections]
# combine_SimCata.py: combine all ImSim outputs to one file.
# optional arguments:
#   -h, --help            show this help message and exit
#   --in_dir IN_DIR       the top directory containing all ImSim outputs.
#   --outfile_path OUTFILE_PATH
#                         the path for the final catalogue.
#   --wanted_cols WANTED_COLS
#                         path to the list of wanted columns.
#                             If not provided, save all columns.
#   --preserve_stars      preserve stars in the final catalogue.
#   --preserve_false_detections
#                         preserve false detections (those not matched with input catalogue) in the final catalogue.
########################

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ parser for command-line interfaces
parser = argparse.ArgumentParser(
    description=f"combine_SimCata.py: combine all ImSim outputs to one file.",
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--in_dir", type=str,
    help="the top directory containing all ImSim outputs.")
parser.add_argument(
    "--outfile_path", type=str,
    help="the path for the final catalogue.")
parser.add_argument(
    "--wanted_cols", type=str,
    help="path to the list of wanted columns.\n\
    If not provided, save all columns.")
parser.add_argument(
    "--preserve_stars", action="store_true",
    help="preserve stars in the final catalogue.")
parser.add_argument(
    "--preserve_false_detections", action="store_true",
    help="preserve false detections (those not matched with input catalogue) in the final catalogue.")

## arg parser
args = parser.parse_args()
in_dir = args.in_dir
outfile_path = args.outfile_path
wanted_cols = args.wanted_cols
if wanted_cols is not None:
    with open(wanted_cols) as f:
        wanted_cols = f.readlines()
    wanted_cols = [x.strip() for x in wanted_cols if x[0]!='#'] 
preserve_stars = args.preserve_stars
preserve_false_detections = args.preserve_false_detections

# +++++++++++++++++++++++++++++ workhorse
print('Combine all ImSim outputs to one file...')
# get subdirs
subdirs = [x for x in glob.glob(os.path.join(in_dir, '*')) if (os.path.isdir(x) and (not 'tmp' in os.path.basename(x)))]
print('Number of subdirs:', len(subdirs))

cata_final = []
for subdir in subdirs:

    run_tag = os.path.basename(subdir)

    # basic info for input shear values
    with open(os.path.join(subdir, 'basic_info.txt'), 'r') as opened_file:
        all_lines = opened_file.readlines()
    useful_line = [line for line in all_lines if 'g_cosmic' in line][0]
    g1g2 = useful_line.split('=')[1]
    g1 = float(g1g2.split()[0])
    g2 = float(g1g2.split()[1])

    # out catalogues from simulation
    file_list = glob.glob(os.path.join(subdir, 'catalogues', '*_combined.*'))
    if not file_list:
        raise Exception(f'Cannot find any combined catalogues in tag {run_tag}!\n\
        make sure taskID=7 is performed in the main pipeline!')
    print(f'Number of files in {run_tag}: {len(file_list)}' )
    for file in file_list:
        file_type = file[-3:]
        if file_type == 'csv':
            cata = pd.read_csv(file)
        elif file_type == 'her':
            cata = pd.read_feather(file)
        elif file_type == 'its':
            with fits.open(file) as hdul:
                cata = Table(hdul[1].data).to_pandas()
        else:
            raise Exception(f'Not supported input file type! {file}')

        # mask for stars & galaxies
        mask_stars = cata['perfect_flag_star']==1
        mask_not_stars = cata['perfect_flag_star']==0
        mask_gals = cata['id_input']>-999

        # select columns
        if wanted_cols is not None:
            cata = cata[wanted_cols]

        # assign cosmic shear
        cata.loc[:, 'g1_in'] = g1
        cata.loc[:, 'g2_in'] = g2

        # assign noise info
        cata.loc[:, 'tile_label'] = re.search(r'tile(.*)_rot', file).group(1)
        # assign rotation info
        cata.loc[:, 'gal_rot'] = float(re.search(r'_rot(\d+)', file).group(1))

        # discard detections as required
        mask_selected = mask_gals
        if preserve_stars:
            print('stars are preserved.')
            mask_selected = mask_selected | mask_stars
        if preserve_false_detections:
            print('false detections are preserved.')
            mask_selected = mask_selected | mask_not_stars

        cata = cata[mask_selected]
        cata_final.append(cata)

cata_final = pd.concat(cata_final)
cata_final.reset_index(drop=True, inplace=True)
print(f'Total number of source {len(cata_final)}')

## check existence
if os.path.isfile(outfile_path):
    os.remove(outfile_path)
## save
if outfile_path[-3:] == 'her':
    cata_final.to_feather(outfile_path)
elif outfile_path[-3:] == 'csv':
    cata_final.to_csv(outfile_path, index=False)
elif outfile_path[-3:] == 'its':
    Table.from_pandas(cata_final).write(outfile_path, format='fits')
else:
    raise Exception(f'Not supported output file type! {outfile_path}')
print(f'Combined catalogue saved as {outfile_path}')
