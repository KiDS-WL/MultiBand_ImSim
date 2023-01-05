# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-08-22 13:03:59
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-01 10:55:35

### a combined catalogue with LF selection

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ general info

# ######### 1sigma
# # main_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7_nU10'
# # main_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7_qU5'
# main_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7_sizeU5'

# ######### 2sigma
# # main_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7_nU20'
# # main_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7_qU10'
# main_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7_sizeU10'

######### 3sigma
# main_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7_nU30'
# main_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7_qU15'
main_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7_sizeU15'


# main directory contains all the outputs
outpath_test = os.path.join(main_dir, f'skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_test.feather')

# the fiducial catalogue
dir_fiducial = '/disks/shear16/ssli/ImSim/output/skills_v07D7/'
outpath_fiducial = os.path.join(main_dir, f'skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_fiducial.feather')

# shear and rotation info
unique_shear_tags = ['m283m283', 'm283p283', 'p283m283', 'p283p283']
unique_rots = [0., 90.]

# parts
parts = ['part0']
# parts = ['part0', 'part1', 'part2', 'part3', 'part4', 'part5']

# +++++++++++++++++++++++++++++ workhorse

## loop over shear tags
cata_final = [] 
cata_final_fiducial = [] 
for part in parts:
    for shear_tag in unique_shear_tags:
        print('Running for', shear_tag)

        subdir = os.path.join(main_dir, part, shear_tag)
        subdir_f = os.path.join(dir_fiducial, part, shear_tag)

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
        file_list = glob.glob(os.path.join(subdir, 'catalogues', f'*_combined.feather'))
        if not file_list:
            raise Exception(f'Cannot find any combined catalogues in tag {shear_tag}!\n\
            make sure taskID=7 is performed in the main pipeline!')
        print(f'Number of files: {len(file_list)}' )
        for file in file_list:

            file_name = os.path.basename(file)
            print('working on file', file_name)

            # useful info
            tile_label = re.search(r'tile(.*)_rot', file_name).group(1)
            gal_rot = float(re.search(r'_rot(\d+)', file_name).group(1))

            for i_file, infile in enumerate([file, os.path.join(subdir_f, 'catalogues', file_name)]):

                cata = pd.read_feather(infile)
        
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

                # >>>>>>>>>> LF-related
                #### a) remove unmeasured
                mask_psf = (cata['psf_Q11_LF_r'] != 0.0) & (cata['psf_Q22_LF_r'] != 0.0)

                #### b) remove binaries
                mask_binary = (np.hypot(cata['e1_LF_r'].values, cata['e2_LF_r'].values) <= 0.8) \
                            | (cata['scalelength_LF_r'] >= \
                                (0.5 * np.exp(0.65788*(24.2 - cata['MAG_AUTO']))))

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
                cata = cata[mask_psf & mask_binary & mask_class & mask_mag & mask_blending]
                cata.reset_index(drop=True, inplace=True)
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
                if i_file == 0:
                    cata_final.append(cata)
                else:
                    cata_final_fiducial.append(cata)
                del cata

cata_final = pd.concat(cata_final, ignore_index=True)
print(f'Total number of source in test {len(cata_final)}')
cata_final.to_feather(outpath_test)
del cata_final
print('combined cata saved to', outpath_test)

cata_final_fiducial = pd.concat(cata_final_fiducial, ignore_index=True)
print(f'Total number of source in fiducial {len(cata_final_fiducial)}')
cata_final_fiducial.to_feather(outpath_fiducial)
print('combined cata saved to', outpath_fiducial)

###### nU10
# Total number of source in test 4180728
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_nU10/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_test.feather
# Total number of source in fiducial 4188088
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_nU10/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_fiducial.feather
# Elapsed:2:18.50,User=151.356,System=306.633,CPU=330.6%.

# ###### nU20
# Total number of source in test 4170830
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_nU20/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_test.feather
# Total number of source in fiducial 4188088
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_nU20/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_fiducial.feather

# ###### nU30
# Total number of source in test 4160472
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_nU30/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_test.feather
# Total number of source in fiducial 4188088
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_nU30/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_fiducial.feather
# Elapsed:2:32.93,User=141.129,System=326.862,CPU=306.0%.

# ###### qU5
# Total number of source in test 4189906
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_qU5/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_test.feather
# Total number of source in fiducial 4188088
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_qU5/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_fiducial.feather
# Elapsed:2:21.03,User=148.715,System=310.525,CPU=325.6%.

# ###### qU10
# Total number of source in test 4190250
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_qU10/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_test.feather
# Total number of source in fiducial 4188088
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_qU10/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_fiducial.feather
# Elapsed:2:17.40,User=147.366,System=334.671,CPU=350.8%.

# # ###### qU15
# Total number of source in test 4189703
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_qU15/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_test.feather
# Total number of source in fiducial 4188088
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_qU15/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_fiducial.feather
# Elapsed:2:24.10,User=141.324,System=286.299,CPU=296.7%.

# ###### sizeU5
# Total number of source in test 4177685
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_sizeU5/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_test.feather
# Total number of source in fiducial 4188088
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_sizeU5/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_fiducial.feather
# Elapsed:2:26.08,User=134.587,System=269.723,CPU=276.7%.

# # ###### sizeU10
# Total number of source in test 4159632
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_sizeU10/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_test.feather
# Total number of source in fiducial 4188088
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_sizeU10/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_fiducial.feather
# Elapsed:2:05.96,User=131.348,System=213.049,CPU=273.4%.

# # ###### sizeU15
# Total number of source in test 4134490
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_sizeU15/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_test.feather
# Total number of source in fiducial 4188088
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_sizeU15/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_fiducial.feather
# Elapsed:2:33.97,User=153.400,System=358.360,CPU=332.3%.
