# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-10-10 14:20:13
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-10 17:22:01

### add flag columns for KiDS-like photo and lensfit selections
######## flag_gaap: 0 for good GAaP photometry
######## flag_asteroid: 1 for asteroids
######## flag_binary: 1 for binary
######## flag_LF_noWeiCut: 0 for all LF selections except for the weight selection

import numpy as np
import pandas as pd

# +++++++++++++++++++++++++++++ general info

# in info
# inpath = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_everything.feather'
inpath = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_everything.feather'
# out info
outpath = inpath.replace('.feather', '_col_flag.feather')

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

# >>>>>>>>>>>> summary
Ntmp = np.sum(cata['flag_gaap']==0)
print('>>> number after GAaP selection', Ntmp, Ntmp/Nori)
Ntmp = np.sum(cata['flag_asteroid']==0)
print('>>> number after asteroid selection', Ntmp, Ntmp/Nori)
Ntmp = np.sum(cata['flag_binary']==0)
print('>>> number after binary selection', Ntmp, Ntmp/Nori)
Ntmp = np.sum(cata['flag_LF_noWeiCut']==0)
print('>>> number after LF_noWeiCut selection', Ntmp, Ntmp/Nori)
Ntmp = np.sum((cata['flag_gaap']==0)&(cata['flag_asteroid']==0)&(cata['flag_binary']==0)&(cata['flag_LF_noWeiCut']==0))
print('>>> number after all selection', Ntmp, Ntmp/Nori)

# >>>>>>>>>> save
cata.to_feather(outpath)
print('final cata saved to', outpath)

# # ########## shark photometry
# number original 89303975
# B_skills_flag.py:49: RuntimeWarning: invalid value encountered in true_divide
#   cata.loc[:, 'R'] = cata['PSFsize'].values / (cata['r_ab'].values**2 + cata['PSFsize'].values)
# >>> number after GAaP selection 88951359 0.9960515083455131
# >>> number after asteroid selection 87865958 0.9838975028827104
# >>> number after binary selection 89274164 0.9996661850718291
# >>> number after LF_noWeiCut selection 48098591 0.5385940659416336
# >>> number after all selection 47929683 0.5367026831672387
# final cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_everything_col_flag.feather
# Elapsed:12:40.85,User=1371.797,System=5888.267,CPU=954.2%.

# ########## kids photometry
# number original 89303975
# B_skills_flag.py:49: RuntimeWarning: invalid value encountered in true_divide
#   cata.loc[:, 'R'] = cata['PSFsize'].values / (cata['r_ab'].values**2 + cata['PSFsize'].values)
# >>> number after GAaP selection 88951359 0.9960515083455131
# >>> number after asteroid selection 87793473 0.9830858368846404
# >>> number after binary selection 89273635 0.9996602614833214
# >>> number after LF_noWeiCut selection 48057297 0.5381316677113197
# >>> number after all selection 47870631 0.5360414360055081
# final cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_everything_col_flag.feather
# Elapsed:10:48.80,User=1478.096,System=5519.657,CPU=1078.5%.