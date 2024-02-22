# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-04-17 09:59:33
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-18 22:04:45

### a combined catalogue with everything detected
### add flag columns for KiDS-like photo and lensfit selections
######## flag_gaap: 0 for good GAaP photometry
######## flag_asteroid: 1 for asteroids
######## flag_binary: 1 for binary
######## flag_LF_noWeiCut: 0 for all LF selections except for the weight selection
## has an option to remove certain tiles

#### applied to sizeU and sizeD, because one tile in skills_v07D7_sizeD/part0/m283p283 failed :(

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

# ## for fiducial run
# test_label = 'skills_v07D7'
# save_label = 'skills_v07D7p1'

## for test runs
# test_label = 'skills_v07D7_sizeU'
test_label = 'skills_v07D7_sizeD'
save_label = test_label

# which tile is discarded?
bad_tiles = ['134.0_0.5']
miss_tile_label = '_'.join(bad_tiles)

# input 
main_dir = f'/disks/shear16/ssli/ImSim/output/{test_label}'
# output
outpath = os.path.join(main_dir, f'{save_label}_LF_321_part0_no_{miss_tile_label}_everything_col_flag.feather')

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

            if tile_label in bad_tiles:
                print(f'!!! tile {tile_label} is discarded !!!')
                continue

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
# Total number of source 14288054
# >>> number after GAaP selection 14235210 0.996301525736115
# >>> number after asteroid selection 14123150 0.9884586102488134
# >>> number after binary selection 14281855 0.9995661410574176
# >>> number after LF_noWeiCut selection 7556168 0.5288451457420303
# >>> number after all selection 7528915 0.5269377481356103
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_part0_no_134.0_0.5_everything_col_flag.feather
# Elapsed:3:25.45,User=268.427,System=526.687,CPU=387.0%.

# ########## sizeU
# Total number of source 14089361
# >>> number after GAaP selection 14036493 0.9962476651709045
# >>> number after asteroid selection 13927695 0.98852566841037
# >>> number after binary selection 14082645 0.999523328275853
# >>> number after LF_noWeiCut selection 7529668 0.5344222495257237
# >>> number after all selection 7502512 0.5324948377715639
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_sizeU/skills_v07D7_sizeU_LF_321_part0_no_134.0_0.5_everything_col_flag.feather
# Elapsed:3:41.99,User=323.736,System=508.576,CPU=374.9%.

# ########### sizeD
# Total number of source 14492890
# >>> number after GAaP selection 14440061 0.9963548333010186
# >>> number after asteroid selection 14323727 0.9883278628348108
# >>> number after binary selection 14486083 0.9995303214196755
# >>> number after LF_noWeiCut selection 7575501 0.5227046503492402
# >>> number after all selection 7547850 0.520796749302589
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7_sizeD/skills_v07D7_sizeD_LF_321_part0_no_134.0_0.5_everything_col_flag.feather
# Elapsed:3:44.14,User=296.354,System=346.931,CPU=286.9%.
