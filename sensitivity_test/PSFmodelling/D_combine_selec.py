# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-07-15 14:11:55
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-05 09:07:29

### a combined catalogue with selections

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd

# +++++++++++++++++++++++++++++ general info

# combined_tag = 'combined_PSFmodelling00_321'
combined_tag = 'combined_PSFmodelling41_321'

# main directory contains all the outputs
main_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7_PSFmodelling/'

# shear and rotation info
unique_shear_tags = ['m283m283', 'm283p283', 'p283m283', 'p283p283']
unique_rots = [0., 90.]

# the fiducial catalogue
dir_fiducial = '/disks/shear16/ssli/ImSim/output/skills_v07D7/'
parts = ['part0', 'part1', 'part2', 'part3', 'part4', 'part5']
# txt file includes wanted columns
wanted_cols = './columns_photometry.txt'
## wanted info
with open(wanted_cols) as f:
    cols_final = f.readlines()
cols_final = [x.strip() for x in cols_final if x[0]!='#'] 

# out info
outpath = os.path.join(main_dir, f'skills_v07D7_LF_321_{combined_tag}_noSG_noWeiCut_newCut.feather')

# +++++++++++++++++++++++++++++ workhorse

## loop over shear tags
cata_final = [] 
for shear_tag in unique_shear_tags:
    print('Running for', shear_tag)

    subdir = os.path.join(main_dir, shear_tag)

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
    file_list = glob.glob(os.path.join(subdir, 'catalogues', f'*_{combined_tag}.feather'))
    if not file_list:
        raise Exception(f'Cannot find any combined catalogues in tag {shear_tag}!\n\
        make sure taskID=7 is performed in the main pipeline!')
    print(f'Number of files: {len(file_list)}' )
    for file in file_list:
        cata = pd.read_feather(file)
        ## MAG_AUTO from fiducial
        cata.drop(columns=['MAG_AUTO'], inplace=True)

        # useful info
        print('working on file', os.path.basename(file))
        print('>>> cata NUMBER', cata['NUMBER'].values)
        tile_label = re.search(r'tile(.*)_rot', file).group(1)
        gal_rot = float(re.search(r'_rot(\d+)', file).group(1))

        # get photometry information from fiducial
        for part in parts:
            dir_fiducial_part = os.path.join(dir_fiducial, part)
            ## get all tile info in this part
            files_tmp = glob.glob(os.path.join(dir_fiducial_part, shear_tag, 'catalogues', 'tile*_rot0_combined_kids_filters.feather'))
            tile_labels_tmp = [re.search(r'tile(.*)_rot', os.path.basename(file)).group(1) for file in files_tmp]
            del files_tmp
            if (tile_label in tile_labels_tmp):
                print('finding in ', part)
                # load the catalogue
                file = os.path.join(dir_fiducial_part, shear_tag, 'catalogues', f'tile{tile_label}_rot{gal_rot:.0f}_combined_kids_filters.feather')
                cata0 = pd.read_feather(file)
                print('>>> cata0 NUMBER', cata0['NUMBER'].values)
                ## used columns
                cata0 = cata0[cols_final]
                ## combine
                cata = pd.concat([cata, cata0], axis=1)
                del cata0
                break

        # some columns for tile info
        ## input cosmic shear
        cata.loc[:, 'g1_in'] = g1
        cata.loc[:, 'g2_in'] = g2
        ## run tag
        cata.loc[:, 'run_tag'] = shear_tag        
        ## tile noise info
        cata.loc[:, 'tile_label'] = tile_label
        ## rotation info
        cata.loc[:, 'gal_rot'] = gal_rot

        # >>>>>>>>>> photometry-related
        ## 1. 9-band photometry cut
        flag_9 = np.zeros(len(cata)) 
        for band in ['u', 'g', 'r', 'i', 'Z', 'Y', 'J', 'H', 'Ks']:
            flag_9 += cata[f'FLAG_GAAP_{band}'].values
        mask_gaap = (flag_9==0)
        del flag_9

        ## 2. remove asteroids
        gmr = np.array(cata['MAG_GAAP_g']) - np.array(cata['MAG_GAAP_r'])
        imr = np.array(cata['MAG_GAAP_i']) - np.array(cata['MAG_GAAP_r'])
        mask_ast = (gmr <= 1.5) | (imr <= 1.5)
        del gmr, imr

        # >>>>>>>>>> LF-related
        #### a) remove unmeasured
        mask_psf = (cata['psf_Q11_LF_r'] != 0.0) & (cata['psf_Q22_LF_r'] != 0.0)

        #### b) remove binaries
        mask_binary = (np.hypot(cata['e1_LF_r'].values, cata['e2_LF_r'].values) <= 0.8) \
                    | (cata['scalelength_LF_r'] >= \
                        (0.5 * np.exp(0.65788*(24.2 - cata['MAG_GAAP_r']))))

        #### c) fitclass cut
        mask_class = (cata['class_LF_r']!=-1) \
                    & (cata['class_LF_r']!=-10) \
                    & (cata['class_LF_r']!=-4) \
                    & (cata['class_LF_r']!=1) \
                    & (cata['class_LF_r']!=2) \
                    & (cata['class_LF_r']!=-7) \
                    & (cata['class_LF_r']!=-3)

        #### d) magnitude cut
        mask_mag = (cata['MAG_AUTO']>20.0)

        #### e) blending cut
        mask_blending = (cata['contamination_radius_LF_r']>4.25)

        ## apply
        cata = cata[mask_gaap & mask_ast &
                    mask_psf & mask_binary & mask_class & mask_mag & mask_blending]
        cata.reset_index(drop=True, inplace=True)
        del mask_gaap, mask_ast
        del mask_psf, mask_binary, mask_class, mask_mag, mask_blending
        print('>>> number after fiducial selection', len(cata))

        # >>>>>>>>>> new selection
        ## 1. avoid negative variance or snr
        cata = cata[(cata['LS_variance_LF_r'].values>0)&(cata['SNR_LF_r'].values>0)]
        cata.reset_index(drop=True, inplace=True)
        ## 2. resolution cut
        ### circularised galaxy size
        emod = np.hypot(cata['e1_LF_r'].values, cata['e2_LF_r'].values)
        cata.loc[:, 'r_ab'] = np.array(cata['scalelength_LF_r'].values) * np.sqrt((1.-emod)/(1.+emod))
        del emod
        ### PSF size
        cata.loc[:, 'PSFsize'] = np.array(
                        (cata['psf_Q11_LF_r'].values*cata['psf_Q22_LF_r'].values \
                            - cata['psf_Q12_LF_r'].values**2.)**0.5)
        ### resolution parameter
        cata.loc[:, 'R'] = np.array(
                        cata['PSFsize'].values\
                            / (cata['r_ab'].values**2 + cata['PSFsize'].values))
        ### cut
        cata = cata[cata['R']<0.9]
        cata.reset_index(drop=True, inplace=True)
        ## 3. size cut
        cata = cata[((cata['scalelength_LF_r'].values 
                                - cata['scalelength_corr_LF_r'].values)>=0.5)]
        cata.reset_index(drop=True, inplace=True)
        print('>>> number after new selection', len(cata))

        # collect
        cata_final.append(cata)
        del cata

cata_final = pd.concat(cata_final, ignore_index=True)
print(f'Total number of source {len(cata_final)}')
cata_final.to_feather(outpath)
print('combined cata saved to', outpath)

# Total number of source 13784445
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_PSFmodelling/skills_v07D7_LF_321_combined_PSFmodelling41_321_noSG_noWeiCut_newCut.feather
# Elapsed:3:51.34,User=245.977,System=599.627,CPU=365.5%.

# >>> number after fiducial selection 83091
# >>> number after new selection 69039
# Total number of source 3253065
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_PSFmodelling/skills_v07D7_LF_321_kids_filters_PSFmodelling00_noSG_noWeiCut_newCut.feather