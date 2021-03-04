# @Author: ssli
# @Date:   2021-02-25, 13:10:25
# @Last modified by:   lshuns
# @Last modified time: 2021-02-26, 8:23:11

### mapping COSMOS shape info to the MICE2 catalogue

import time
import random
random.seed(94120)

import numpy as np
import pandas as pd
from scipy import stats
import pyvinecopulib as pv
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import Planck15 as astropy_Planck15
from astropy.cosmology import FlatLambdaCDM

import os
import sys
# Self-defined package
sys.path.insert(0,os.path.realpath('../SURFS/'))
from shape_mapping_func import TwoBinningMappingFunc, SimpleMappingFunc, ZBinningMappingFunc

start_time = time.time()

# ++++++++++++++ general setting
# overall cut
mag_r_max_cut = 26.
z_max_cut = 1.42
Re_cut = [1e-2, 10.] # arcsec
# region for careful mapping
mag_r_min = 20.
mag_r_max1 = 24.
mag_r_max2 = 25.

# output
outfile = '/disks/shear10/ssli/ImSim/input/MICE2_cata/MICE2_ra35-37_dec5-7_COSMOS_morphology.feather'

# ++++++++++++++ simulation catalogue
infile = '/disks/shear10/ssli/ImSim/input/MICE2_cata/MICE2_ra35-37_dec5-7_corrected.fits'
with fits.open(infile) as hdul:
    simu_cata = hdul[1].data
print('simulation original', len(simu_cata))
## general selection
mask_mag = (simu_cata['cosmos_subaru_r_true_lensed']<mag_r_max_cut)
mask_z = (simu_cata['z_cgal_v']<z_max_cut)
mask_tmp = mask_mag & mask_z
simu_cata = simu_cata[mask_tmp]
print('simulation selected', len(simu_cata))
## useful columns
# cols_in = ['z_cgal_v',
#         'cosmos_subaru_r_true_lensed', 'cosmos_subaru_i_true_lensed',
#         'cosmos_subaru_r_abs_mag_evo', 'cosmos_subaru_i_abs_mag_evo']
# cols_out = ['redshift', 'mag_r', 'mag_i', 'mag_r_abs', 'mag_i_abs']
cols_in = ['z_cgal_v',
        'cosmos_subaru_r_true_lensed', 'cosmos_subaru_i_true_lensed',
        'cosmos_subaru_r_true_lensed', 'cosmos_subaru_i_true_lensed']
cols_out = ['redshift', 'mag_r', 'mag_i', 'mag_r_abs', 'mag_i_abs']
simu_ori_cata = pd.DataFrame({'index': np.array(simu_cata['unique_gal_id']).astype(int)})
for i_col, col_in in enumerate(cols_in):
    simu_ori_cata.loc[:, cols_out[i_col]] = np.array(simu_cata[col_in]).astype(float)
## size info
simu_ori_cata.loc[:, 'bf'] = np.array(simu_cata['bulge_fraction']).astype(float)
bl = np.array(simu_cata['bulge_length']).astype(float)
dl = np.array(simu_cata['disk_length']).astype(float)
size_arcsec = np.where(bl>dl, bl, dl)
## distance for size transfer
# cosmo_mice2 = FlatLambdaCDM(H0=70, Om0=0.25, Ob0=0.044)
# simu_ori_cata.loc[:, 'Dist_A'] = cosmo_mice2.angular_diameter_distance(simu_ori_cata['redshift'].values).value # Mpc
simu_ori_cata.loc[:, 'Dist_A'] = astropy_Planck15.angular_diameter_distance(simu_ori_cata['redshift'].values).value # Mpc
simu_ori_cata.loc[:, 'size_physical'] = simu_ori_cata['Dist_A'].values * (size_arcsec/3600./180.*np.pi) * 1e3 # kpc
## three groups: key, bright, faint
# +++ key1
mask_key1 = (simu_ori_cata['mag_r']>=mag_r_min) & (simu_ori_cata['mag_r']<mag_r_max1)
simu_ori_cata_key1 = simu_ori_cata[mask_key1].copy()
simu_ori_cata_key1.reset_index(drop=True, inplace=True)
print('simulation key1', len(simu_ori_cata_key1))
# +++ key2
mask_key2 = (simu_ori_cata['mag_r']>=mag_r_max1) & (simu_ori_cata['mag_r']<mag_r_max2)
simu_ori_cata_key2 = simu_ori_cata[mask_key2].copy()
simu_ori_cata_key2.reset_index(drop=True, inplace=True)
print('simulation key2', len(simu_ori_cata_key2))
# +++ bright
mask_bright = (simu_ori_cata['mag_r']<mag_r_min)
simu_ori_cata_bright = simu_ori_cata[mask_bright].copy()
simu_ori_cata_bright.reset_index(drop=True, inplace=True)
print('simulation bright', len(simu_ori_cata_bright))
# +++ faint
mask_faint = (simu_ori_cata['mag_r']>=mag_r_max2)
simu_ori_cata_faint = simu_ori_cata[mask_faint].copy()
simu_ori_cata_faint.reset_index(drop=True, inplace=True)
print('simulation faint', len(simu_ori_cata_faint))

# ++++++++++++++ data
### cosmos catalogue
infile = '/disks/shear10/ssli/ImSim/input/COSMOS_cata/cosmos_shape_z_ugriZYJHKs_selected.feather'
cosmos_cata = pd.read_feather(infile)
print('COSMOS15 shape_z_photo_cata', len(cosmos_cata))
# ## absolute magnitude
# infile = '/disks/shear10/ssli/ImSim/input/COSMOS_cata/outside/cosmos_2015_dr_2_1_COSMOS_flagged_AbsMag.fits'
# with fits.open(infile) as hdul:
#     cosmos_cata_abs = hdul[1].data
# cosmos_cata_abs = pd.DataFrame({'id_laigle': cosmos_cata_abs['number'].astype(int),
#                                 'm_r': cosmos_cata_abs['m_r'].astype(float),
#                                 'm_i': cosmos_cata_abs['m_i'].astype(float)})
# print('COSMOS15 abs_mag_cata', len(cosmos_cata_abs))
# ### merge
# cosmos_cata = cosmos_cata.merge(cosmos_cata_abs, how='inner', on='id_laigle')
# print('COSMOS15 combined', len(cosmos_cata))
## use apparent magnitude
cosmos_cata.loc[:, 'm_r'] = cosmos_cata['r_mag_auto'].values
## general cut
mask_mag = (cosmos_cata['r_mag_auto']<mag_r_max_cut)
mask_z = (cosmos_cata['Z_optimal']<z_max_cut)
mask_size = (cosmos_cata['shape/Re']>Re_cut[0]) & (cosmos_cata['shape/Re']<Re_cut[1])
mask_tmp = mask_mag & mask_z & mask_size
cosmos_cata = cosmos_cata[mask_tmp]
cosmos_cata.reset_index(drop=True, inplace=True)
print('COSMOS15 selected', len(cosmos_cata))
## physical size for ranking
dist_z = astropy_Planck15.angular_diameter_distance(cosmos_cata['Z_optimal'].values).value # Mpc
cosmos_cata.loc[:, 'Re_physical'] = dist_z * (cosmos_cata['shape/Re'].values/3600./180.*np.pi) * 1e3 # kpc
## three groups: key, bright, faint
# +++ key1
mask_key1 = (cosmos_cata['r_mag_auto']>=mag_r_min) & (cosmos_cata['r_mag_auto']<mag_r_max1)
cosmos_cata_key1 = cosmos_cata[mask_key1].copy()
cosmos_cata_key1.reset_index(drop=True, inplace=True)
print('COSMOS15 key1', len(cosmos_cata_key1))
# +++ key2
mask_key2 = (cosmos_cata['r_mag_auto']>=mag_r_max1) & (cosmos_cata['r_mag_auto']<mag_r_max2)
cosmos_cata_key2 = cosmos_cata[mask_key2].copy()
cosmos_cata_key2.reset_index(drop=True, inplace=True)
print('COSMOS15 key2', len(cosmos_cata_key2))
# +++ bright
mask_bright = (cosmos_cata['r_mag_auto']<mag_r_min)
cosmos_cata_bright = cosmos_cata[mask_bright].copy()
cosmos_cata_bright.reset_index(drop=True, inplace=True)
print('COSMOS15 bright', len(cosmos_cata_bright))
# +++ faint
mask_faint = (cosmos_cata['r_mag_auto']>=mag_r_max2)
cosmos_cata_faint = cosmos_cata[mask_faint].copy()
cosmos_cata_faint.reset_index(drop=True, inplace=True)
print('COSMOS15 faint', len(cosmos_cata_faint))

# ++++++++++++++ running
# +++ key1
simu_ori_cata_key1 = TwoBinningMappingFunc(simu_ori_cata_key1, cosmos_cata_key1, N_zbins=20, dmagbins = 0.5, Niter = 10000)
# +++ key2
simu_ori_cata_key2 = TwoBinningMappingFunc(simu_ori_cata_key2, cosmos_cata_key2, N_zbins=20, dmagbins = 0.5, Niter = 10000)
# +++ bright
simu_ori_cata_bright = SimpleMappingFunc(simu_ori_cata_bright, cosmos_cata_bright, simu_ori_cata=None, mask_z_simu=None, mask_mag_simu=None)
# +++ faint
simu_ori_cata_faint = ZBinningMappingFunc(simu_ori_cata_faint, cosmos_cata_faint, N_zbins=20)

### save
simu_final_cata = pd.concat([simu_ori_cata_key1, simu_ori_cata_key2, simu_ori_cata_bright, simu_ori_cata_faint])
simu_final_cata.reset_index(drop=True, inplace=True)
print('simulation final', len(simu_final_cata))
simu_final_cata.to_feather(outfile)
# Table.from_pandas(simu_final_cata).write(outfile, format='fits')
print('Final catalogue saved as', outfile)
print('All finished in', time.time()-start_time)
# simulation original 421893
# simulation selected 421892
# simulation key 405496
# simulation bright 6746
# simulation faint 9650
# COSMOS15 shape_z_photo_cata 101593
# COSMOS15 abs_mag_cata 536077
# COSMOS15 combined 101593
# COSMOS15 selected 89343
# COSMOS15 key 81418
# COSMOS15 bright 1294
# COSMOS15 faint 6631
