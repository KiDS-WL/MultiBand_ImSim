# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-08-10 13:49:16
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-28 15:32:35

### combine and select the results
####### photometry based on the fiducial full run, LF based on new run

import os
import re
import sys
import glob
import argparse

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ general info

# the fiducial catalogues
ori_dir = f'/disks/shear16/ssli/ImSim/output/skills_v07D7'

# main directory contains all the outputs
# main_dir = f'/disks/shear16/ssli/ImSim/output/skills_v07D7v'
main_dir = f'/disks/shear16/ssli/ImSim/output/skills_v07D7v_extra'

# does main dir include sub-dir from different input
split_ran_tags = ['part0', 'part1', 'part2', 'part3', 'part4', 'part5']

# # shear and rotation info
# ### constant shear
# shear_type = 'constant'
# unique_shear_tags = ['m283m283', 'm283p283', 'p283m283', 'p283p283']
# unique_shear_tags_ori = ['m283m283', 'm283p283', 'p283m283', 'p283p283']
# ### constant shear

#### variable shear
shear_type = 'variable'
# unique_shear_tags = ['CS0_rot0',  'CS0_rot90', 'CS0_rot180', 'CS0_rot270',
#                     'CS1_rot0',  'CS1_rot90', 'CS1_rot180', 'CS1_rot270',
#                     'CS2_rot0',  'CS2_rot90', 'CS2_rot180', 'CS2_rot270',
#                     'CS3_rot0',  'CS3_rot90', 'CS3_rot180', 'CS3_rot270']
unique_shear_tags = ['CS0_rot45',  'CS0_rot135', 'CS0_rot225', 'CS0_rot315',
                    'CS1_rot45',  'CS1_rot135', 'CS1_rot225', 'CS1_rot315',
                    'CS2_rot45',  'CS2_rot135', 'CS2_rot225', 'CS2_rot315',
                    'CS3_rot45',  'CS3_rot135', 'CS3_rot225', 'CS3_rot315']
unique_shear_tags_ori = ['m283p283', 'p283p283', 'p283m283', 'm283m283',
                        'm283p283', 'p283p283', 'p283m283', 'm283m283',
                        'm283p283', 'p283p283', 'p283m283', 'm283m283',
                        'm283p283', 'p283p283', 'p283m283', 'm283m283']
#### variable shear

# out info
# outpath = os.path.join(main_dir, f'skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_const.feather')
outpath = os.path.join(main_dir, f'skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_var.feather')

# +++++++++++++++++++++++++++++ workhorse

## loop over
cata_final = [] 
for part_tag in split_ran_tags:
    for i_shear, shear_tag in enumerate(unique_shear_tags):

        # the original dir
        subdir_ori = os.path.join(ori_dir, part_tag, unique_shear_tags_ori[i_shear])
        print('ori dir', subdir_ori)

        # new dir
        subdir = os.path.join(main_dir, part_tag, shear_tag)
        print('new dir', subdir)

        # basic info for input shear values
        if shear_type == 'constant':
            with open(os.path.join(subdir, 'basic_info.txt'), 'r') as opened_file:
                all_lines = opened_file.readlines()
            useful_line = [line for line in all_lines if 'g_cosmic' in line][0]
            try:
                g1 = float(useful_line.split()[2])
            except ValueError:
                g1 = float(useful_line.split()[2][:-1])
            g2 = float(useful_line.split()[3])

        # main catalogues from simulation
        file_list = glob.glob(os.path.join(subdir, 'catalogues', f'*_combined.*'))
        if not file_list:
            raise Exception(f'Cannot find any combined catalogues in tag {shear_tag}!\n\
            make sure taskID=7 is performed in the main pipeline!')
        print(f'Number of files: {len(file_list)}' )
        for file in file_list:
            file_name = os.path.basename(file)
            print(">>>>>>>>> file", file)

            # >>>>>>>>>>>> the target catalogue
            file_type = file[-3:]
            if file_type == 'csv':
                cata0 = pd.read_csv(file)
            elif file_type == 'her':
                cata0 = pd.read_feather(file)
            elif file_type == 'its':
                with fits.open(file) as hdul:
                    cata0 = Table(hdul[1].data).to_pandas()
            else:
                raise Exception(f'Not supported input file type! {file}')
            # some columns for tile info
            if shear_type == 'constant':
                ## input cosmic shear
                cata0.loc[:, 'g1_in'] = g1
                cata0.loc[:, 'g2_in'] = g2
            else:
                cata0.rename(columns={"gamma1_input": "g1_in", "gamma2_input": "g2_in"}, inplace=True)
            ## run tag
            cata0.loc[:, 'run_tag'] = shear_tag        
            ## tile noise info
            cata0.loc[:, 'tile_label'] = re.search(r'tile(.*)_rot', file_name).group(1)
            ## rotation info
            cata0.loc[:, 'gal_rot'] = float(re.search(r'_rot(\d+)', file_name).group(1))
            print('+++++++++ number ori blending', len(cata0))
            ## only galaxies
            cata0 = cata0[cata0['id_input']>0]
            cata0.reset_index(drop=True, inplace=True)
            print('+++++++++ number gal blending', len(cata0))

            # >>>>>>>>>>>> information from original catalogue
            file = os.path.join(subdir_ori, 'catalogues', file_name.replace('_combined.feather', '_combined_kids_filters.feather'))
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
            print('>>> number ori whole', len(cata))

            # >>>>>>>>>> photometry-related selections
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
            ## apply
            cata = cata[mask_gaap & mask_ast]
            cata.reset_index(drop=True, inplace=True)
            del mask_gaap, mask_ast
            print('>>> number after selection', len(cata))

            # >>>>>>>>>>>>> used columns
            cata = cata[['id_input', 'Z_B', 'MAG_GAAP_r']]

            # >>>>>>>>>>>>> combine with the blending cata
            cata0 = cata0.merge(cata, on='id_input')
            del cata
            print('+++++++++ number merged', len(cata0))
            # sys.exit()

            # >>>>>>>>>> LF-related selections
            #### a) remove unmeasured
            mask_psf = (cata0['psf_Q11_LF_r'] != 0.0) & (cata0['psf_Q22_LF_r'] != 0.0)

            #### b) remove binaries
            mask_binary = (np.hypot(cata0['e1_LF_r'].values, cata0['e2_LF_r'].values) <= 0.8) \
                        | (cata0['scalelength_LF_r'] >= \
                            (0.5 * np.exp(0.65788*(24.2 - cata0['MAG_GAAP_r']))))

            #### c) fitclass cut
            mask_class = (cata0['class_LF_r']!=-1) \
                        & (cata0['class_LF_r']!=-10) \
                        & (cata0['class_LF_r']!=-4) \
                        & (cata0['class_LF_r']!=1) \
                        & (cata0['class_LF_r']!=2) \
                        & (cata0['class_LF_r']!=-7) \
                        & (cata0['class_LF_r']!=-3)

            #### d) magnitude cut
            mask_mag = (cata0['MAG_AUTO']>20.0)

            #### e) blending cut
            mask_blending = (cata0['contamination_radius_LF_r']>4.25)

            ## apply
            cata0 = cata0[mask_psf & mask_binary & mask_class & mask_mag & mask_blending]
            cata0.reset_index(drop=True, inplace=True)
            del mask_psf, mask_binary, mask_class, mask_mag, mask_blending
            print('>>> number after fiducial selection', len(cata0))

            # >>>>>>>>>> new selection
            ## 1. avoid negative variance or snr
            cata0 = cata0[(cata0['LS_variance_LF_r'].values>0)&(cata0['SNR_LF_r'].values>0)]
            cata0.reset_index(drop=True, inplace=True)
            ## 2. resolution cut
            ### circularised galaxy size
            emod = np.hypot(cata0['e1_LF_r'].values, cata0['e2_LF_r'].values)
            cata0.loc[:, 'r_ab'] = np.array(cata0['scalelength_LF_r'].values) * np.sqrt((1.-emod)/(1.+emod))
            del emod
            ### PSF size
            cata0.loc[:, 'PSFsize'] = np.array(
                            (cata0['psf_Q11_LF_r'].values*cata0['psf_Q22_LF_r'].values \
                                - cata0['psf_Q12_LF_r'].values**2.)**0.5)
            ### resolution parameter
            cata0.loc[:, 'R'] = np.array(
                            cata0['PSFsize'].values\
                                / (cata0['r_ab'].values**2 + cata0['PSFsize'].values))
            ### cut
            cata0 = cata0[cata0['R']<0.9]
            cata0.reset_index(drop=True, inplace=True)
            ## 3. size cut
            cata0 = cata0[((cata0['scalelength_LF_r'].values 
                                    - cata0['scalelength_corr_LF_r'].values)>=0.5)]
            cata0.reset_index(drop=True, inplace=True)
            print('>>> number after new selection', len(cata0))

            # collect
            cata_final.append(cata0)
            del cata0

cata_final = pd.concat(cata_final, ignore_index=True)
print(f'Total number of source {len(cata_final)}')
cata_final.to_feather(outpath)
print('combined cata saved to', outpath)

# ###### varShear extra
# Total number of source 64110761
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7v_extra/skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_var.feather
# Elapsed:29:34.07,User=2083.549,System=5307.801,CPU=416.6%.

# ###### varShear
# Total number of source 64116903
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_var.feather
# Elapsed:38:33.98,User=2171.809,System=3954.969,CPU=264.7%.

# ###### constShear
# Total number of source 16024951
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_const.feather
# Elapsed:10:14.86,User=547.687,System=910.834,CPU=237.2%.
