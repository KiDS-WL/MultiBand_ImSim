# @Author: lshuns
# @Date:   2021-03-03, 18:21:04
# @Last modified by:   ssli
# @Last modified time: 2021-03-05, 16:09:16

### a simple script to combine catalogues produced by the main pipeline
###     it can be used to combine catalogues from different running_tags
###             and assign photo-z based on input id of galaxies
###                 in which case at least one running_tag should contain photo-z
###                     and same galaxies are simulated in different running_tags
###                     and same noise, same rotations
###                     that is the only difference is the shear input
###     it removes false detections (those not matched with the input catalogue)
###     it only preserves specified columns
###     it can generate mask (MASK_gold) for desired selections

import os
import re
import glob

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ the only variables you want to modify
# where to save the final catalogue
### supported formats: fits, csv, feather
###             guess from the file suffix
out_path = '/disks/shear15/ssli/ImSim/output/test_surfs_onePFSmean_combined.fits'

# where to find the input catalogues
in_dir = '/disks/shear15/ssli/ImSim/output/'

# a list of running tags being combined
run_tag_list = ['test_surfs_onePFSmean', 'test_surfs_onePFSmean1', 'test_surfs_onePFSmean2', 'test_surfs_onePFSmean3',
                    'test_surfs_onePFSmean4', 'test_surfs_onePFSmean5', 'test_surfs_onePFSmean6', 'test_surfs_onePFSmean7']

# assign photo-z ?
### NOTE: if that is set to True, the first run_tag should contain the photo-z info
assign_photoz = True

# apply selections ?
### supported choices: none (without selection), KiDS-lensfit (KiDS-lensfit-like selection)
### it produce a MASK_gold column with 0 for selected sample
mask_gold = 'KiDS-lensfit'

# desired columns
###  its too large and unnecessary to preserve all the columns
###  'id_input' is the reference, therefore, should always be provided
cols_final = ['NUMBER', 'X_WORLD', 'Y_WORLD', 'MAG_AUTO', 'perfect_flag_star', 'id_input', 'Z_B',
                'e1_LF_r', 'e2_LF_r', 'SNR_LF_r', 'scalelength_LF_r', 'weight_LF_r', 'psf_Q11_LF_r', 'psf_Q22_LF_r', 'psf_Q12_LF_r', 'class_LF_r', 'contamination_radius_LF_r']


# +++++++++++++++++++++++++++++ workhorse

# === assign photoz
### should match catalogue in tile level
if assign_photoz:
    print('Assign photoz to catalogues')
    print('     NOTE: catalogues from different run_tag should have same input galaxies and same rotations and same noise')
    print('             the only difference is the shear input')

    # the fist tag is the reference with photoz
    run_tag0 = run_tag_list[0]
    file_list = glob.glob(os.path.join(in_dir, run_tag0, 'catalogues', '*_combined.*'))
    if not file_list:
        raise Exception('Cannot find any catalogues from the main pipeline!\n\
        make sure the main pipeline ended properly to the end (taskID=7 is required)!')
    # get rotations
    rotations = np.unique([re.search(r'_rot(\d+)', file).group(1) for file in file_list])
    # get tiles (noise realization)
    tiles = np.unique([re.search(r'tile(.*)_rot', file).group(1) for file in file_list])
    print(f'Number of files in each run_tag: {len(file_list)} (={len(tiles)} tiles * {len(rotations)} rotations)' )

    # loop over all tags
    cata_final = []
    for run_tag in run_tag_list:

        # basic info for input shear values
        opened_file = open(os.path.join(in_dir, run_tag, 'basic_info.txt'), 'r')
        all_lines = opened_file.readlines()
        opened_file.close()
        useful_line = [line for line in all_lines if 'g_cosmic' in line][0]
        try:
            g1 = float(useful_line.split()[2])
        except ValueError:
            g1 = float(useful_line.split()[2][:-1])
        g2 = float(useful_line.split()[3])

        # main catalogues within the tag
        if not 'cata_z_tag' in locals():
            cata_z_tag = [] # to save z info, not used if assign_photoz=False
        i_file = 0
        for tile in tiles:
            for rot in rotations:

                file = glob.glob(os.path.join(in_dir, run_tag, 'catalogues', f'tile{tile}_rot{rot}_combined.*'))[0]

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
                try:
                    cata = cata[cols_final]
                except KeyError:
                    cols_final.remove('Z_B')
                    cata = cata[cols_final]

                # discard false-detections
                mask_true = (cata['id_input']>-999)
                cata = cata[mask_true]

                # assign cosmic shear
                cata.loc[:, 'g1_in'] = g1
                cata.loc[:, 'g2_in'] = g2

                # for fast join
                cata.set_index('id_input', inplace=True)

                # assign photoz
                if run_tag == run_tag0:
                    cata_withz = cata[['Z_B']]
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
    # easy combine without assign photoz
    cata_final = []
    for run_tag in run_tag_list:

        # basic info for input shear values
        opened_file = open(os.path.join(in_dir, run_tag, 'basic_info.txt'), 'r')
        all_lines = opened_file.readlines()
        opened_file.close()
        useful_line = [line for line in all_lines if 'g_cosmic' in line][0]
        try:
            g1 = float(useful_line.split()[2])
        except ValueError:
            g1 = float(useful_line.split()[2][:-1])
        g2 = float(useful_line.split()[3])

        # main catalogues from simulation
        file_list = glob.glob(os.path.join(in_dir, run_tag, 'catalogues', '*_combined.*'))
        if not file_list:
            raise Exception('Cannot find any catalogues from the main pipeline!\n\
            make sure the main pipeline ended properly to the end (taskID=7 is required)!')
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

            cata = cata[cols_final]

            # assign cosmic shear
            cata.loc[:, 'g1_in'] = g1
            cata.loc[:, 'g2_in'] = g2

            # discard false-detections
            mask_true = (cata['id_input']>-999)
            cata = cata[mask_true]

            cata_final.append(cata)

    cata_final = pd.concat(cata_final)
    print(f'Total number of source {len(cata_final)}')

# apply selections
if mask_gold == 'none':
    pass
elif mask_gold == 'KiDS-lensfit':

    ## mask
    mask_gal = (cata_final['perfect_flag_star']==0)
    # lensfit cuts with 0 means no issue and -9 means large galaxies
    fitcuts = ((cata_final['class_LF_r']==0) | (cata_final['class_LF_r']==-9))
    # remove invalid measurements
    weight_cuts = (cata_final['weight_LF_r']>0)
    # remove potentially blended sources
    blend_cuts = (cata_final['contamination_radius_LF_r']>4.25)
    # remove unresolved binary stars
    binary_star_cuts = ((np.hypot(cata_final['e1_LF_r'], cata_final['e2_LF_r'])<=0.8) | (cata_final['scalelength_LF_r']>=0.5*np.exp(0.65788*(24.2-cata_final['MAG_AUTO']))))

    # apply
    cata_final.loc[(mask_gal & fitcuts & weight_cuts & blend_cuts & binary_star_cuts), 'MASK_gold'] = 0
    print(f'Number after lensfit selection {np.sum(mask_gal & fitcuts & weight_cuts & blend_cuts & binary_star_cuts)}')
else:
    raise Exception(f'Unsupported mask_gold type! {mask_gold}')

cata_final.reset_index(drop=False, inplace=True)
# dummy values for nan
cata_final.fillna(-999, inplace=True)
cata_final = cata_final.astype({'MASK_gold': int})

# save
file_type = out_path[-3:]
if file_type == 'her':
    cata_final.to_feather(out_path)
elif file_type == 'csv':
    cata_final.to_csv(out_path, index=False)
elif file_type == 'its':
    Table.from_pandas(cata_final).write(out_path, format='fits')
else:
    raise Exception(f'Not supported output file type! {out_path}')

print(f'Combined catalogue saved as {out_path}')
