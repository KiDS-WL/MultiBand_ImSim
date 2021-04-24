# @Author: lshuns
# @Date:   2021-04-20, 12:25:08
# @Last modified by:   ssli
# @Last modified time: 2021-04-24, 14:13:31

### mapping COSMOS shape info to the SURFS catalogue

import os
import time

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

from shape_mapping_func import TwoBinningMappingFunc

# ++++++++++++++ I/O
# output
outdir = '/disks/shear10/ssli/ImSim/input/SURFS_cata/'
outfile_list = [f'CUM_r_Ks_patch{i_patch}.feather' for i_patch in range(4)]

# original simulation cata
indir = '/disks/shear10/ssli/ImSim/input/SURFS_cata/'
infile_list = [f'SURFS_SHARK_LC_patch{i_patch}.fits' for i_patch in range(4)]
simu_col_z = 'zobs'
simu_col_mag = 'r_SDSS_apparent'
simu_col_colour = ['r_SDSS_apparent', 'K_VISTA_apparent']
simu_col_index = 'index'
simu_col_Upsilon_corr = 'Upsilon_corr'

# selected cosmos catalogue
infile_cosmos_shape = '/disks/shear10/ssli/ImSim/input/COSMOS_cata/cosmos_shape_z_ugriZYJHKs_selected.feather'
cosmos_col_z = 'Z_optimal'
cosmos_col_mag = 'r_mag_auto'
cosmos_col_colour = ['r_mag_auto', 'ks_mag_auto']
cosmos_cols_shapes = ['shape/Re', 'N_GALFIT_HI', 'BA_GALFIT_HI']

# ++++++++++++++ general settings
# overall cut
mag_cut = 27.
z_cut = 2.5
Re_cut = [1e-2, 10.] # arcsec

# binning
N_zbins = 30
N_magbins = 40

# good names
goodname_z = 'z'
goodname_mag = 'binning_mag'
goodname_col_order = 'order_col'
goodnames_shapes = ['size', 'N_GALFIT_HI', 'BA_GALFIT_HI']

# ++++++++++++++ workhorse

# data
cosmos_cata = pd.read_feather(infile_cosmos_shape)
print('COSMOS15 shape_z_photo_cata', len(cosmos_cata))
### colours
cosmos_cata.loc[:, goodname_col_order] = np.array(cosmos_cata[cosmos_col_colour[0]] - cosmos_cata[cosmos_col_colour[1]])
### desired columns
cosmos_col_all = [cosmos_col_z, cosmos_col_mag, goodname_col_order] + cosmos_cols_shapes
cosmos_cata = cosmos_cata[cosmos_col_all]
### rename
cols_goodnames_all = [goodname_z, goodname_mag] + goodnames_shapes
cols_oldnames_all = [cosmos_col_z, cosmos_col_mag] + cosmos_cols_shapes
cols_rename = {}
for i_col, col_name in enumerate(cols_oldnames_all):
     cols_rename[col_name] = cols_goodnames_all[i_col]
cosmos_cata.rename(columns=cols_rename, inplace=True)
### general selection
mask_mag = (cosmos_cata[goodname_mag]<mag_cut)
mask_z = (cosmos_cata[goodname_z]<z_cut)
mask_size = (cosmos_cata['size']>Re_cut[0]) & (cosmos_cata['size']<Re_cut[1])
mask_tmp = mask_mag & mask_z & mask_size
cosmos_cata = cosmos_cata[mask_tmp]
cosmos_cata.reset_index(drop=True, inplace=True)
print('COSMOS15 selected', len(cosmos_cata))

# loop over simulation catalogue
for i_file, infile in enumerate(infile_list):

    start_time = time.time()

    infile_simu = os.path.join(indir, infile)
    with fits.open(infile_simu) as hdul:
        simu_cata = Table(hdul[1].data).to_pandas()
    print('simulation original', len(simu_cata))
    ### colours
    simu_cata.loc[:, goodname_col_order] = np.array(simu_cata[simu_col_colour[0]] - simu_cata[simu_col_colour[1]])
    ### desired columns
    cols_simu = [simu_col_z, simu_col_mag, goodname_col_order]
    if simu_col_index is not None:
        cols_simu += [simu_col_index]
    else:
        cols_simu += ['index']
        simu_cata.loc[:, 'index'] = np.arange(len(simu_cata)).astype(int)
    ### modify magnitudes
    if simu_col_Upsilon_corr is not None:
        simu_cata.loc[:, simu_col_mag] = np.array(simu_cata[simu_col_mag] + simu_cata[simu_col_Upsilon_corr])
    ### get desired columns
    simu_cata = simu_cata[cols_simu]
    ### rename
    cols_rename = {simu_col_z: goodname_z, simu_col_mag: goodname_mag}
    simu_cata.rename(columns=cols_rename, inplace=True)
    ### general selection
    mask_mag = (simu_cata[goodname_mag]<mag_cut)
    mask_z = (simu_cata[goodname_z]<z_cut)
    mask_tmp = mask_mag & mask_z
    simu_cata = simu_cata[mask_tmp]
    simu_cata.reset_index(drop=True, inplace=True)
    print('simulation selected', len(simu_cata))

    # shape mapping
    simu_cata = TwoBinningMappingFunc(simu_cata, cosmos_cata, N_zbins=N_zbins, N_magbins=N_magbins)
    ### rename size
    simu_cata.rename(columns={'size': 'Re_arcsec'}, inplace=True)

    ### save
    outpath = os.path.join(outdir, outfile_list[i_file])
    simu_cata.to_feather(outpath)
    # Table.from_pandas(simu_cata).write(outpath, format='fits')
    print('Final catalogue saved as', outpath)
    print(f'File {i_file} finished in', (time.time()-start_time)/60., 'min')

# COSMOS15 shape_z_photo_cata 101593
# COSMOS15 selected 99225
# simulation original 19365969
# simulation selected 8794876
# Final catalogue saved as /disks/shear10/ssli/ImSim/input/SURFS_cata/CM_r_Ks_patch0.fits
# File 0 finished in 92.76339008808137 min
# simulation original 19116707
# simulation selected 8715550
# Final catalogue saved as /disks/shear10/ssli/ImSim/input/SURFS_cata/CM_r_Ks_patch1.fits
# File 1 finished in 93.33005142211914 min
# simulation original 19039242
# simulation selected 8686136
# Final catalogue saved as /disks/shear10/ssli/ImSim/input/SURFS_cata/CM_r_Ks_patch2.fits
# File 2 finished in 90.5914852976799 min
# simulation original 19371513
# simulation selected 8829566
# Final catalogue saved as /disks/shear10/ssli/ImSim/input/SURFS_cata/CM_r_Ks_patch3.fits
# File 3 finished in 92.59597713549933 min
