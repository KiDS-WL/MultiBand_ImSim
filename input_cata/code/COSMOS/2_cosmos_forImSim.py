# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-01-27 18:07:58
# @Last Modified by:   lshuns
# @Last Modified time: 2021-02-01 15:10:30

### create two mock catalogues for MultiBand_ImSim
###     1. good shape measurement & have all photometry & have redshifts
###         g band from B, V bands (http://classic.sdss.org/dr5/algorithms/sdssUBVRITransform.html)
###     2. good shape measurement & r band photometry

import numpy as np
import pandas as pd
import scipy.spatial as sst
import matplotlib.pyplot as plt

from scipy import stats
from astropy.io import fits

# ++++++++++++ Files 
cata_dir = '/disks/shear10/ssli/ImSim/input/COSMOS_cata/'

## inputs
cosmos_file = cata_dir + 'cosmos_shape_z_uBVriZYJHKs.feather'

## output
outpath_9_photo = cata_dir + 'cosmos_shape_z_ugriZYJHKs_selected.feather'
outpath_r_photo = cata_dir + 'cosmos_shape_r_selected.feather'

# ++++++++++++ cosmos
cosmos_cata = pd.read_feather(cosmos_file)
print('Original number', len(cosmos_cata))

### shape from pixel to arcsec
cosmos_cata.loc[:, 'shape/Re'] = cosmos_cata['RE_GALFIT_HI'].values * 0.05

### make sersic index discrete to speed the ImSim
sersic_n_discrete = np.array([np.exp(np.linspace(np.log(0.5), np.log(6.0), 100))]).T
# build KDTree
kdt = sst.cKDTree(sersic_n_discrete, leafsize=100)
# query
bn_ori = np.array([cosmos_cata['N_GALFIT_HI']]).T
dist, ind = kdt.query(bn_ori, k=1, distance_upper_bound=np.inf)
bn_new = sersic_n_discrete[ind].flatten()
cosmos_cata.loc[:, 'shape/sersic_n'] = bn_new

### select
###### 0. discard too small or too big galaxies
mask_re = (cosmos_cata['shape/Re']>=1e-2) & (cosmos_cata['shape/Re']<=10.)
print('     good Re number', np.sum(mask_re))
###### 1. good shape
mask_galfit = (cosmos_cata['FLAG_GALFIT_HI']==0)
print('     good shape number', np.sum(mask_galfit))
###### 2. with photometry
mask_mag = (cosmos_cata['u_mag_auto']>0) & (cosmos_cata['b_mag_auto']>0) & (cosmos_cata['v_mag_auto']>0) &\
            (cosmos_cata['r_mag_auto']>0) & (cosmos_cata['ip_mag_auto']>0) & (cosmos_cata['zpp_mag_auto']>0) &\
            (cosmos_cata['y_mag_auto']>0) & (cosmos_cata['j_mag_auto']>0) & (cosmos_cata['h_mag_auto']>0) & (cosmos_cata['ks_mag_auto']>0)
print('     with photometry', np.sum(mask_mag))
###### 3. with r photometry
mask_mag_r = (cosmos_cata['r_mag_auto']>0)
print('     with r photometry', np.sum(mask_mag_r))
### apply
cosmos_cata_9 = cosmos_cata[mask_re & mask_galfit & mask_mag].copy()
cosmos_cata_9.reset_index(drop=True, inplace=True)
print('     after shape & 9 photo selection', len(cosmos_cata_9))
cosmos_cata_r = cosmos_cata[mask_re & mask_galfit & mask_mag_r].copy()
cosmos_cata_r.reset_index(drop=True, inplace=True)
print('     after shape & r photo selection', len(cosmos_cata_r))

# ++++++++++++ calculate g band
B = np.array(cosmos_cata_9['b_mag_auto'])
V = np.array(cosmos_cata_9['v_mag_auto'])

# 1. Jester 2005
g_Jester_1 = V + 0.64*(B-V) - 0.13 # Rc-Ic < 1.15 and U-B < 0
g_Jester_2 = V + 0.60*(B-V) - 0.12 # Rc-Ic < 1.15

# 2. Lupton 2005
C1 = 0.3130
C2 = 0.2271
C3 = 0.5784
C4 = 0.0038
g_Lupton = (C3/C1*B + V - C2*C3/C1 + C4) / (C3/C1+1)

# save back
cosmos_cata_9.loc[:, 'g_Jester_1'] = g_Jester_1
cosmos_cata_9.loc[:, 'g_Jester_2'] = g_Jester_2
cosmos_cata_9.loc[:, 'g_Lupton'] = g_Lupton

# ++++++++++++++ useful columns for ImSim
cols = ['OBJNO', 'RA', 'DEC', 'Z_optimal', 
            'shape/Re', 'shape/sersic_n', 'N_GALFIT_HI', 'BA_GALFIT_HI', 
            'u_mag_auto', 'g_Jester_1', 'g_Jester_2', 'g_Lupton', 'r_mag_auto', 'ip_mag_auto', 'zpp_mag_auto', 'y_mag_auto', 'j_mag_auto', 'h_mag_auto', 'ks_mag_auto']
cosmos_cata_9 = cosmos_cata_9[cols].copy()
print('     final (9 photo)', len(cosmos_cata_9))
cosmos_cata_9.to_feather(outpath_9_photo)
print('saved to', outpath_9_photo)

cols = ['OBJNO', 'RA', 'DEC', 'Z_optimal', 
            'shape/Re', 'shape/sersic_n', 'N_GALFIT_HI', 'BA_GALFIT_HI', 
            'r_mag_auto']
cosmos_cata_r = cosmos_cata_r[cols].copy()
print('     final (r photo)', len(cosmos_cata_r))
cosmos_cata_r.to_feather(outpath_r_photo)
print('saved to', outpath_r_photo)


######### info
# Original number 304688
#      good Re number 291092
#      good shape number 134836
#      with photometry 195507
#      with r photometry 204964
#      after shape & 9 photo selection 101593
#      after shape & r photo selection 102699
#      final (9 photo) 101593
#      final (r photo) 102699
