# @Author: lshuns
# @Date:   2021-02-19, 10:49:19
# @Last modified by:   lshuns
# @Last modified time: 2021-02-21, 12:13:51

### mapping COSMOS shape info to the SURFS catalogue

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

from shape_mapping_func import TwoBinningMappingFunc, SimpleMappingFunc, ZBinningMappingFunc

start_time = time.time()

# ++++++++++++++ general setting
# overall cut
mag_r_max_cut = 27.
z_max_cut = 2.5
Re_cut = [1e-2, 10.] # arcsec
# region for careful mapping
mag_r_min = 19.
mag_r_max = 25.

# binning
### log binning for redshift
N_zbins = 20
### correlation-adaptive binning for magnitude
###### first guess of intervals
dmagbins = 0.5
###### maximum number of iterations
Niter = 10000

# output
outfile = '/disks/shear10/ssli/ImSim/input/SURFS_cata/foursqdeg_testcols_COSMOS_morphology.feather'

# ++++++++++++++ simulation catalogue
infile = '/disks/shear10/ssli/ImSim/input/SURFS_cata/SURFS_SHARK_LC_4sqdeg_testcols.fits'
with fits.open(infile) as hdul:
    simu_cata = hdul[1].data
print('simulation original', len(simu_cata))
## general selection
mask_mag = (simu_cata['r_SDSS_apparent']<mag_r_max_cut)
mask_z = (simu_cata['zobs']<z_max_cut)
mask_tmp = mask_mag & mask_z
simu_cata = simu_cata[mask_tmp]
print('simulation selected', len(simu_cata))
## useful columns
cols_in = ['zobs',
        'r_SDSS_apparent', 'i_SDSS_apparent',
        'r_SDSS_absolute', 'i_SDSS_absolute']
cols_out = ['redshift', 'mag_r', 'mag_i', 'mag_r_abs', 'mag_i_abs']
simu_ori_cata = pd.DataFrame({'index': np.arange(len(simu_cata)).astype(int)})
for i_col, col_in in enumerate(cols_in):
    simu_ori_cata.loc[:, cols_out[i_col]] = simu_cata[col_in].astype(float)
## size info
simu_ori_cata.loc[:, 'bf'] = simu_cata['mstars_bulge']/(simu_cata['mstars_bulge']+simu_cata['mstars_disk'])
simu_ori_cata.loc[:, 'size_physical'] = np.where(simu_cata['rstar_bulge_apparent']>simu_cata['rstar_disk_apparent'], simu_cata['rstar_bulge_apparent'], simu_cata['rstar_disk_apparent'])
## distance for size transfer
simu_ori_cata.loc[:, 'Dist_A'] = astropy_Planck15.angular_diameter_distance(simu_ori_cata['redshift'].values).value # Mpc
## three groups: key, bright, faint
# +++ key
mask_key = (simu_ori_cata['mag_r']>=mag_r_min) & (simu_ori_cata['mag_r']<mag_r_max)
simu_ori_cata_key = simu_ori_cata[mask_key].copy()
simu_ori_cata_key.reset_index(drop=True, inplace=True)
print('simulation key', len(simu_ori_cata_key))
# +++ bright
mask_bright = (simu_ori_cata['mag_r']<mag_r_min)
simu_ori_cata_bright = simu_ori_cata[mask_bright].copy()
simu_ori_cata_bright.reset_index(drop=True, inplace=True)
print('simulation bright', len(simu_ori_cata_bright))
# +++ faint
mask_faint = (simu_ori_cata['mag_r']>=mag_r_max)
simu_ori_cata_faint = simu_ori_cata[mask_faint].copy()
simu_ori_cata_faint.reset_index(drop=True, inplace=True)
print('simulation faint', len(simu_ori_cata_faint))

# ++++++++++++++ data
### cosmos catalogue
infile = '/disks/shear10/ssli/ImSim/input/COSMOS_cata/cosmos_shape_z_ugriZYJHKs_selected.feather'
cosmos_cata = pd.read_feather(infile)
print('COSMOS15 shape_z_photo_cata', len(cosmos_cata))
## absolute magnitude
infile = '/disks/shear10/ssli/ImSim/input/COSMOS_cata/outside/cosmos_2015_dr_2_1_COSMOS_flagged_AbsMag.fits'
with fits.open(infile) as hdul:
    cosmos_cata_abs = hdul[1].data
cosmos_cata_abs = pd.DataFrame({'id_laigle': cosmos_cata_abs['number'].astype(int),
                                'm_r': cosmos_cata_abs['m_r'].astype(float),
                                'm_i': cosmos_cata_abs['m_i'].astype(float)})
print('COSMOS15 abs_mag_cata', len(cosmos_cata_abs))
### merge
cosmos_cata = cosmos_cata.merge(cosmos_cata_abs, how='inner', on='id_laigle')
print('COSMOS15 combined', len(cosmos_cata))
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
# +++ key
mask_key = (cosmos_cata['r_mag_auto']>=mag_r_min) & (cosmos_cata['r_mag_auto']<mag_r_max)
cosmos_cata_key = cosmos_cata[mask_key].copy()
cosmos_cata_key.reset_index(drop=True, inplace=True)
print('COSMOS15 key', len(cosmos_cata_key))
# +++ bright
mask_bright = (cosmos_cata['r_mag_auto']<mag_r_min)
cosmos_cata_bright = cosmos_cata[mask_bright].copy()
cosmos_cata_bright.reset_index(drop=True, inplace=True)
print('COSMOS15 bright', len(cosmos_cata_bright))
# +++ faint
mask_faint = (cosmos_cata['r_mag_auto']>=mag_r_max)
cosmos_cata_faint = cosmos_cata[mask_faint].copy()
cosmos_cata_faint.reset_index(drop=True, inplace=True)
print('COSMOS15 faint', len(cosmos_cata_faint))

# ++++++++++++++ running
# +++ key
simu_ori_cata_key = TwoBinningMappingFunc(simu_ori_cata_key, cosmos_cata_key, N_zbins=20, dmagbins = 0.5, Niter = 10000)
# +++ bright
simu_ori_cata_bright = SimpleMappingFunc(simu_ori_cata_bright, cosmos_cata_bright, simu_ori_cata=None, mask_z_simu=None, mask_mag_simu=None)
# +++ faint
simu_ori_cata_faint = ZBinningMappingFunc(simu_ori_cata_faint, cosmos_cata_faint, N_zbins=20)

### save
simu_final_cata = pd.concat([simu_ori_cata_key, simu_ori_cata_bright, simu_ori_cata_faint])
simu_final_cata.reset_index(drop=True, inplace=True)
print('simulation final', len(simu_final_cata))
simu_final_cata.to_feather(outfile)
# Table.from_pandas(simu_final_cata).write(outfile, format='fits')
print('Final catalogue saved as', outfile)
print('All finished in', time.time()-start_time)
# All finished in 2990.9944610595703
# simulation original 2856332
# simulation selected 881394
# simulation key 251519
# simulation bright 1690
# simulation faint 628185
# COSMOS15 shape_z_photo_cata 101593
# COSMOS15 abs_mag_cata 536077
# COSMOS15 combined 101593
# COSMOS15 selected 99225
# COSMOS15 key 91622
# COSMOS15 bright 284
# COSMOS15 faint 7319
# simulation final 881394
