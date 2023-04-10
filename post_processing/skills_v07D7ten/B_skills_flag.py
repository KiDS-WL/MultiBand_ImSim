# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-04-10 13:09:29
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-10 14:13:29

### add flag columns for KiDS-like photo and lensfit selections
#### following DR5 convention (Wright et al. 2023)
######## flag_gaap: 0 for good GAaP photometry
######## flag_asteroid: 1 for asteroids
######## flag_binary: 1 for binary
######## flag_LF_noWeiCut: 0 for all LF selections except for the weight selection

import numpy as np
import pandas as pd

# +++++++++++++++++++++++++++++ general info

# in info
inpath = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7ten_LF_321_kidsPhotometry_everything.feather'
# out info
outpath = inpath.replace('.feather', '_col_flag_dr5.feather')

bands = ['u', 'g', 'r', 'i', 'i2', 'Z', 'Y', 'J', 'H', 'Ks']

# txt file includes wanted columns
wanted_cols = './wanted_cols.txt'

# +++++++++++++++++++++++++++++ workhorse

## load and select catalogue
cata = pd.read_feather(inpath)
Nori = len(cata)
print('number original', Nori)

## wanted info
with open(wanted_cols) as f:
    cols_final = f.readlines()
cols_final = [x.strip() for x in cols_final if x[0]!='#'] 
cata = cata[cols_final]

## calculate the resolution factor
### circularised galaxy size
emod = np.hypot(cata['e1_LF_r'].values, cata['e2_LF_r'].values)
cata.loc[:, 'r_ab'] = cata['scalelength_LF_r'].values * np.sqrt((1.-emod)/(1.+emod))
del emod
### PSF size
cata.loc[:, 'PSFsize'] = (cata['psf_Q11_LF_r'].values*cata['psf_Q22_LF_r'].values
                                - cata['psf_Q12_LF_r'].values**2.)**0.5
### resolution parameter
cata.loc[:, 'R'] = cata['PSFsize'].values / (cata['r_ab'].values**2 + cata['PSFsize'].values)

## columns for saving flags
cata_flags = pd.DataFrame(data=np.ones((len(cata), 4)),
                        columns=('flag_gaap', 'flag_asteroid', 'flag_binary',
                            'flag_LF_noWeiCut'), dtype=np.int32)
cata = pd.concat([cata, cata_flags], axis=1)
del cata_flags

# >>>>>>>>>> 1. GAaP photometry cut
flag_gaap = np.zeros(len(cata)) 
for band in bands:
    flag_gaap += cata[f'FLAG_GAAP_{band}'].values
cata.loc[(flag_gaap==0), 'flag_gaap'] = 0
del flag_gaap
print('>>> removed by GAaP selection', np.sum(cata['flag_gaap']==1)/Nori)

# >>>>>>>>>> 2. remove asteroids
gmr = np.array(cata['MAG_GAAP_g']) - np.array(cata['MAG_GAAP_r'])
imr = np.array(cata['MAG_GAAP_i']) - np.array(cata['MAG_GAAP_r'])
i2mr = np.array(cata['MAG_GAAP_i2']) - np.array(cata['MAG_GAAP_r'])
cata.loc[((gmr <= 1.5) | (imr <= 1.5) | (i2mr <= 1.5)), 'flag_asteroid'] = 0
del gmr, imr, i2mr
print('>>> removed by asteroids selection', np.sum(cata['flag_asteroid']==1)/Nori)

# >>>>>>>>>> 3. remove binaries
mask_binary = (np.hypot(cata['e1_LF_r'].values, cata['e2_LF_r'].values) <= 0.8) \
            | (
                (cata['scalelength_LF_r'].values >= (0.5 * 10**((24.9 - cata['MAG_GAAP_r'].values)/3.5)))
                & (cata['MAG_GAAP_r'].values > 21.5)
                )
cata.loc[mask_binary, 'flag_binary'] = 0
del mask_binary
print('>>> removed by binary selection', np.sum(cata['flag_binary']==1)/Nori)

# >>>>>>>>>> 4. LF-related
#### a) remove unmeasured
mask_psf = (cata['psf_Q11_LF_r'] != 0.0) & (cata['psf_Q22_LF_r'] != 0.0)
print('>>> removed by unmeasured', 1 - np.sum(mask_psf)/Nori)
#### b) fitclass cut
mask_class = (cata['class_LF_r']!=-10)\
            &  (cata['class_LF_r']!=-7)\
            & (cata['class_LF_r']!=-4)\
            & (cata['class_LF_r']!=-3)\
            & (cata['class_LF_r']!=-1)\
            & (cata['class_LF_r']!=1)\
            & (cata['class_LF_r']!=2)
print('>>> removed by fitclass', 1 - np.sum(mask_class)/Nori)
#### d) magnitude cut
mask_mag = (cata['MAG_AUTO']>20.0)
print('>>> removed by mag cut', 1 - np.sum(mask_mag)/Nori)
#### e) blending cut
mask_blending = (cata['contamination_radius_LF_r']>4.25)
print('>>> removed by blending', 1 - np.sum(mask_blending)/Nori)
#### f) avoid negative variance or snr
mask_snr = (cata['LS_variance_LF_r'].values>0)&(cata['SNR_LF_r'].values>0)
#### g) size cut
mask_size = (cata['scalelength_LF_r'].values 
                - cata['scalelength_corr_LF_r'].values)>=0.5
#### h) resolution cut
mask_R = cata['R']<0.9
print('>>> removed by resolution cut', 1 - np.sum(mask_R)/Nori)
#### combine
cata.loc[mask_psf & mask_class & mask_mag & mask_blending & mask_snr & mask_size & mask_R, 'flag_LF_noWeiCut'] = 0
del mask_psf, mask_class, mask_mag, mask_blending, mask_snr, mask_size, mask_R

# >>>>>>>>>> save
cata.to_feather(outpath)
print('final cata saved to', outpath)

# ########## kids photometry
# number original 89303769
# B_skills_flag.py:51: RuntimeWarning: invalid value encountered in true_divide
#   cata.loc[:, 'R'] = cata['PSFsize'].values / (cata['r_ab'].values**2 + cata['PSFsize'].values)
# >>> removed by GAaP selection 0.00394887028788225
# >>> removed by asteroids selection 0.006760812077259584
# >>> removed by binary selection 0.0004953094420908484
# >>> removed by unmeasured 0.1704961634933907
# >>> removed by fitclass 0.27580753058697893
# >>> removed by mag cut 0.03883841677499633
# >>> removed by blending 0.22313481528422385
# >>> removed by resolution cut 0.4078605797701551
# final cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7ten_LF_321_kidsPhotometry_everything_col_flag_dr5.feather
# Elapsed:12:54.89,User=1313.204,System=6196.600,CPU=969.1%.