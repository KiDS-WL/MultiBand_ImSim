# @Author: lshuns
# @Date:   2021-03-03, 18:21:04
# @Last modified by:   lshuns
# @Last modified time: 2021-04-26, 20:44:24

### a simple script to combine catalogues produced by the main pipeline
###     it can be used to combine catalogues from different running_tags
###             and assign photo-z based on input id of galaxies
###                 in which case at least one running_tag should contain photo-z
###                     and same galaxies are simulated in different running_tags
###                     and same noise, same rotations
###                     that is the only difference is the shear input
###     it removes false detections (those not matched with the input catalogue)
###     it includes a tile_label column for noise info
###     it includes a gal_rot column for rotation info
###     it only preserves specified columns
###     it can generate mask (MASK_gold) for desired selections

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
    description=f"combine_SimCata.py: combine catalogues from different cosmic shear realization.",
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--main_dir", type=str,
### NOTE: different sub-directories have identical galaxies&noise but different shear inputs
    help="the top directory containing all simulation outputs.")
parser.add_argument(
    "--outfile_name", type=str, default='final_combined.fits',
    help="The file name for the combined file. \n\
    File format inferred from the suffix. Supported formats: fits, feather, csv")
parser.add_argument(
    "--photoz_tag", type=str, default=None,
    help="The tag name of which contains photoz info. \n\
    Not mandatory. \n\
    If not provided, no assignment of photoz will be performed.")
parser.add_argument(
    "--save_all_cols", action="store_true",
    help="preserve all the columns.")
parser.add_argument(
    "--ori_cata_file", type=str, default=None,
    help="The original catalogue used by ImSim. \n\
    Not mandatory. \n\
    If not provided, no cross-match will be performed.")
parser.add_argument(
    "--ori_cata_id", type=str, default='index',
    help="The unique galaxy ID in the original catalogue. \n\
    Not used if ori_cata_file is not provided.")

## arg parser
args = parser.parse_args()
main_dir = args.main_dir
outfile_name = args.outfile_name
photoz_tag = args.photoz_tag
save_all_cols = args.save_all_cols
ori_cata_file = args.ori_cata_file
ori_cata_id = args.ori_cata_id

# +++++++++++++++++++++++++++++ the only variables you may want to modify
# desired columns
###  its too large and unnecessary to preserve all the columns
###  'id_input' is the reference, therefore, should always be provided
if not save_all_cols:
    cols_final = ['NUMBER', 'X_WORLD', 'Y_WORLD', 'MAG_AUTO', 'FLUX_RADIUS', 'perfect_flag_star', 'id_input', 'Z_B', 'mask_meas_9bands',
                'e1_LF_r', 'e2_LF_r', 'SNR_LF_r', 'scalelength_LF_r', 'oldweight_LF_r', 'weight_global_LF_r', 'psf_Q11_LF_r', 'psf_Q22_LF_r', 'psf_Q12_LF_r', 'class_LF_r', 'contamination_radius_LF_r']

# +++++++++++++++++++++++++++++ workhorse

# === assign photoz
### should match catalogue in tile level
if photoz_tag is not None:
    print('Assign photoz to catalogues')
    print(f'     Use photoz from {photoz_tag}')
    print('     NOTE: catalogues from different run_tag should have same input galaxies and same rotations and same noise')
    print('             the only difference is the shear input')

    # the tag with photoz is the reference
    file_list = glob.glob(os.path.join(main_dir, photoz_tag, 'catalogues', '*_combined.*'))
    if not file_list:
        raise Exception(f'Cannot find any combined catalogues in tag {photoz_tag}!\n\
        make sure taskID=7 is performed in the main pipeline!')
    # get rotations
    rotations = np.unique([re.search(r'_rot(\d+)', file).group(1) for file in file_list])
    # get tiles (noise realization)
    tiles = np.unique([re.search(r'tile(.*)_rot', file).group(1) for file in file_list])
    print(f'Number of files in each run_tag: {len(file_list)} (={len(tiles)} tiles * {len(rotations)} rotations)' )

    # loop over all tags
    subdirs = [x for x in glob.glob(os.path.join(main_dir, '*')) if (os.path.isdir(x) and os.path.basename(x)!=photoz_tag)]
    subdirs = [os.path.join(main_dir, photoz_tag)] + subdirs
    print('Number of subdirs:', len(subdirs))
    print('     NOTE: this should match with the number of input shear realizations.')
    cata_final = []
    for subdir in subdirs:

        run_tag = os.path.basename(subdir)

        # basic info for input shear values
        with open(os.path.join(subdir, 'basic_info.txt'), 'r') as opened_file:
            all_lines = opened_file.readlines()
        useful_line = [line for line in all_lines if 'g_cosmic' in line][0]
        try:
            g1 = float(useful_line.split()[2])
        except ValueError:
            g1 = float(useful_line.split()[2][:-1])
        g2 = float(useful_line.split()[3])

        # main catalogues within the tag
        if not 'cata_z_tag' in locals():
            cata_z_tag = [] # to save z info
        i_file = 0
        for tile in tiles:
            for rot in rotations:

                file = glob.glob(os.path.join(subdir, 'catalogues', f'tile{tile}_rot{rot}_combined.*'))[0]

                # load catalogue
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
                if 'cols_final' in locals():
                    try:
                        cata = cata[cols_final]
                    except KeyError:
                        cols_final.remove('Z_B')
                        cols_final.remove('mask_meas_9bands')
                        cata = cata[cols_final]

                # discard false-detections
                mask_true = (cata['id_input']>-999)
                cata = cata[mask_true]

                # assign cosmic shear
                cata.loc[:, 'g1_in'] = g1
                cata.loc[:, 'g2_in'] = g2

                # assign noise info
                cata.loc[:, 'tile_label'] = tile
                # assign rotation info
                cata.loc[:, 'gal_rot'] = float(rot)

                # for fast join
                cata.set_index('id_input', drop=False, inplace=True)

                # assign photoz
                if run_tag == photoz_tag:
                    cata_withz = cata[['Z_B', 'mask_meas_9bands']]
                    cata_z_tag.append(cata_withz)
                else:
                    cata = cata.join(cata_z_tag[i_file], how='left')
                # iterating
                i_file += 1

                cata_final.append(cata)
        print(f'All catalogues collected for {run_tag}')

    cata_final = pd.concat(cata_final)
    print(f'Total number of source {len(cata_final)}')

else:
    print('easy combine without assigning photoz')
    subdirs = [x for x in glob.glob(os.path.join(main_dir, '*')) if os.path.isdir(x)]
    print('Number of subdirs:', len(subdirs))
    print('     NOTE: this should match with the number of input shear realizations.')
    cata_final = []
    for subdir in subdirs:

        run_tag = os.path.basename(subdir)

        # basic info for input shear values
        with open(os.path.join(subdir, 'basic_info.txt'), 'r') as opened_file:
            all_lines = opened_file.readlines()
        useful_line = [line for line in all_lines if 'g_cosmic' in line][0]
        try:
            g1 = float(useful_line.split()[2])
        except ValueError:
            g1 = float(useful_line.split()[2][:-1])
        g2 = float(useful_line.split()[3])

        # main catalogues from simulation
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

            if 'cols_final' in locals():
                cata = cata[cols_final]

            # assign cosmic shear
            cata.loc[:, 'g1_in'] = g1
            cata.loc[:, 'g2_in'] = g2

            # assign noise info
            cata.loc[:, 'tile_label'] = re.search(r'tile(.*)_rot', file).group(1)
            # assign rotation info
            cata.loc[:, 'gal_rot'] = float(re.search(r'_rot(\d+)', file).group(1))
            
            # discard false-detections
            mask_true = (cata['id_input']>-999)
            cata = cata[mask_true]

            cata_final.append(cata)

    cata_final = pd.concat(cata_final)
    print(f'Total number of source {len(cata_final)}')

cata_final.reset_index(drop=True, inplace=True)
# dummy values for nan
cata_final.fillna(-999, inplace=True)

# cross-match with the original cata
if ori_cata_file is not None:
    # load catalogue
    file_type = ori_cata_file[-3:]
    if file_type == 'csv':
        ori_cata = pd.read_csv(ori_cata_file)
    elif file_type == 'her':
        ori_cata = pd.read_feather(ori_cata_file)
    elif file_type == 'its':
        with fits.open(ori_cata_file) as hdul:
            ori_cata = Table(hdul[1].data).to_pandas()
    else:
        raise Exception(f'Not supported input file type! {ori_cata_file}')

    # cross match based on the index
    ## for fast join
    ori_cata.set_index(ori_cata_id, drop=False, inplace=True)
    cata_final.set_index('id_input', drop=False, inplace=True)
    ## join
    cata_final = cata_final.join(ori_cata, how='left')
    ## reset index
    cata_final.reset_index(drop=True, inplace=True)

# save
out_path = os.path.join(main_dir, outfile_name)
## check existence
if os.path.isfile(out_path):
    os.remove(out_path)

if out_path[-3:] == 'her':
    cata_final.to_feather(out_path)
elif out_path[-3:] == 'csv':
    cata_final.to_csv(out_path, index=False)
elif out_path[-3:] == 'its':
    Table.from_pandas(cata_final).write(out_path, format='fits')
else:
    raise Exception(f'Not supported output file type! {out_path}')

print(f'Combined catalogue saved as {out_path}')
