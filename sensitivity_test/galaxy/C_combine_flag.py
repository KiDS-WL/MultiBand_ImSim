# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-04-03 09:33:15
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-18 18:23:54

### a combined catalogue with everything detected
### add flag columns for KiDS-like photo and lensfit selections
######## flag_gaap: 0 for good GAaP photometry
######## flag_asteroid: 1 for asteroids
######## flag_binary: 1 for binary
######## flag_LF_noWeiCut: 0 for all LF selections except for the weight selection

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ general info

# the running label
# test_label = 'skills_v07D7'
# save_label = 'skills_v07D7p1'

# test_label = 'skills_v07D7_nD'
# test_label = 'skills_v07D7_nU'
# test_label = 'skills_v07D7_qD'
test_label = 'skills_v07D7_qU'
save_label = test_label

# input 
main_dir = f'/disks/shear16/ssli/ImSim/output/{test_label}'
# output
outpath = os.path.join(main_dir, f'{save_label}_LF_321_part0_everything_col_flag.feather')

# shear and rotation info
unique_shear_tags = ['m283m283', 'm283p283', 'p283m283', 'p283p283']
unique_rots = [0., 90.]

# parts
parts = ['part0']

# +++++++++++++++++++++++++++++ workhorse

## loop over shear tags
cata_final = [] 
for part in parts:
    for shear_tag in unique_shear_tags:
        print('Running for', shear_tag)

        subdir = os.path.join(main_dir, part, shear_tag)

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

            # load the catalogue
            cata = pd.read_feather(file)

            # calculate resolution
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

            # columns for saving flags
            cata_flags = pd.DataFrame(data=np.ones((len(cata), 4)),
                                    columns=('flag_gaap', 'flag_asteroid', 'flag_binary',
                                        'flag_LF_noWeiCut'), dtype=np.int32)
            cata = pd.concat([cata, cata_flags], axis=1)
            del cata_flags

            # >>>>>>>>>> 1. 9-band photometry cut
            flag_9 = np.zeros(len(cata)) 
            for band in ['u', 'g', 'r', 'i', 'Z', 'Y', 'J', 'H', 'Ks']:
                flag_9 += cata[f'FLAG_GAAP_{band}'].values
            cata.loc[(flag_9==0), 'flag_gaap'] = 0
            del flag_9

            # >>>>>>>>>> 2. remove asteroids
            gmr = np.array(cata['MAG_GAAP_g']) - np.array(cata['MAG_GAAP_r'])
            imr = np.array(cata['MAG_GAAP_i']) - np.array(cata['MAG_GAAP_r'])
            cata.loc[((gmr <= 1.5) | (imr <= 1.5)), 'flag_asteroid'] = 0
            del gmr, imr

            # >>>>>>>>>> 3. remove binaries
            mask_binary = (np.hypot(cata['e1_LF_r'].values, cata['e2_LF_r'].values) <= 0.8) \
                        | (cata['scalelength_LF_r'] >= \
                            (0.5 * np.exp(0.65788*(24.2 - cata['MAG_GAAP_r']))))
            cata.loc[mask_binary, 'flag_binary'] = 0
            del mask_binary

            # >>>>>>>>>> 4. LF-related
            #### a) remove unmeasured
            mask_psf = (cata['psf_Q11_LF_r'] != 0.0) & (cata['psf_Q22_LF_r'] != 0.0)
            #### b) fitclass cut
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
            #### f) avoid negative variance or snr
            mask_snr = (cata['LS_variance_LF_r'].values>0)&(cata['SNR_LF_r'].values>0)
            #### g) size cut
            mask_size = (cata['scalelength_LF_r'].values 
                            - cata['scalelength_corr_LF_r'].values)>=0.5
            #### h) resolution cut
            mask_R = cata['R']<0.9
            #### combine
            cata.loc[mask_psf & mask_class & mask_mag & mask_blending & mask_snr & mask_size & mask_R, 'flag_LF_noWeiCut'] = 0
            del mask_psf, mask_class, mask_mag, mask_blending, mask_snr, mask_size, mask_R

            cata_final.append(cata)
            del cata

# combine
cata_final = pd.concat(cata_final, ignore_index=True)
Nori = len(cata_final)
print(f'Total number of source {Nori}')

# summary
Ntmp = np.sum(cata_final['flag_gaap']==0)
print('>>> number after GAaP selection', Ntmp, Ntmp/Nori)
Ntmp = np.sum(cata_final['flag_asteroid']==0)
print('>>> number after asteroid selection', Ntmp, Ntmp/Nori)
Ntmp = np.sum(cata_final['flag_binary']==0)
print('>>> number after binary selection', Ntmp, Ntmp/Nori)
Ntmp = np.sum(cata_final['flag_LF_noWeiCut']==0)
print('>>> number after LF_noWeiCut selection', Ntmp, Ntmp/Nori)
Ntmp = np.sum((cata_final['flag_gaap']==0)&(cata_final['flag_asteroid']==0)&(cata_final['flag_binary']==0)&(cata_final['flag_LF_noWeiCut']==0))
print('>>> number after all selection', Ntmp, Ntmp/Nori)

# save
cata_final.to_feather(outpath)
print('combined cata saved to', outpath)

# # ###### fiducial
# Total number of source 15061327
# >>> number after GAaP selection 15005385 0.9962857190471995
# >>> number after asteroid selection 14888088 0.9884977598587429
# >>> number after binary selection 15054772 0.9995647793849771
# >>> number after LF_noWeiCut selection 7944570 0.5274814098385886
# >>> number after all selection 7916316 0.5256054795171766
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_part0_everything_col_flag.feather
# Elapsed:5:40.41,User=315.753,System=204.902,CPU=152.9%.

# ###### nD
# Total number of source 15143303
# >>> number after GAaP selection 15087384 0.9963073445733734
# >>> number after asteroid selection 14968451 0.9884535097792073
# >>> number after binary selection 15136198 0.9995308157011716
# >>> number after LF_noWeiCut selection 7969609 0.5262794385082303
# >>> number after all selection 7941115 0.5243978146643437
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_nD/skills_v07D7_nD_LF_321_part0_everything_col_flag.feather
# Elapsed:9:18.38,User=425.660,System=258.756,CPU=122.5%.

# ###### nU
# Total number of source 14979282
# >>> number after GAaP selection 14923311 0.9962634390620325
# >>> number after asteroid selection 14806303 0.9884521167302945
# >>> number after binary selection 14972178 0.9995257449589373
# >>> number after LF_noWeiCut selection 7918292 0.5286162581090336
# >>> number after all selection 7890047 0.5267306537122407
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_nU/skills_v07D7_nU_LF_321_part0_everything_col_flag.feather
# Elapsed:7:52.91,User=431.569,System=263.835,CPU=147.0%.

# ###### qD
# Total number of source 15117568
# >>> number after GAaP selection 15061631 0.9962998678094255
# >>> number after asteroid selection 14942757 0.9884365659873334
# >>> number after binary selection 15110368 0.9995237329178873
# >>> number after LF_noWeiCut selection 7945713 0.5255946591409412
# >>> number after all selection 7917304 0.5237154547609775
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_qD/skills_v07D7_qD_LF_321_part0_everything_col_flag.feather
# Elapsed:3:08.92,User=324.802,System=584.108,CPU=481.1%.

# ###### qU
# Total number of source 15005291
# >>> number after GAaP selection 14949338 0.9962711153019291
# >>> number after asteroid selection 14832459 0.9884819294740769
# >>> number after binary selection 14998237 0.9995298991535719
# >>> number after LF_noWeiCut selection 7941955 0.5292769730357112
# >>> number after all selection 7913639 0.5273899053340585
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_qU/skills_v07D7_qU_LF_321_part0_everything_col_flag.feather
# Elapsed:3:01.28,User=349.037,System=829.198,CPU=649.9%.

