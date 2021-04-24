# @Author: lshuns
# @Date:   2021-04-07, 18:42:59
# @Last modified by:   ssli
# @Last modified time: 2021-04-24, 14:41:02

### clean catalogue to be ready for ImSim

import os

import numpy as np
import pandas as pd
import scipy.spatial as sst

from astropy.io import fits
from astropy.table import Table

# ++++++++++++++ I/O
# output
outdir = '/disks/shear10/ssli/ImSim/input/SURFS_cata/'
outfile_list = [f'SS_CUM_r_Ks_patch{i_patch}.fits' for i_patch in range(4)]

# input
indir = '/disks/shear10/ssli/ImSim/input/SURFS_cata/'
## original catalogue
orifile_list = [f'SURFS_SHARK_LC_patch{i_patch}.fits' for i_patch in range(4)]
## shape mapped catalogue
shapefile_list = [f'CM_r_Ks_patch{i_patch}.feather' for i_patch in range(4)]

# ++++++++++++++ general setting
# overall cut
Re_cut = [1e-2, 10.] # arcsec

## useful columns
cols_used_mags = ['u_SDSS_apparent', 'g_SDSS_apparent', 'r_SDSS_apparent', 'i_SDSS_apparent', 'z_SDSS_apparent',
'Y_VISTA_apparent', 'J_VISTA_apparent', 'H_VISTA_apparent', 'K_VISTA_apparent']

cols_used_others = ['index', 'zobs', 'ra', 'dec',
        'Re_arcsec', 'BA', 'shape/sersic_n']

col_Upsilon_corr = 'Upsilon_corr'

# ++++++++++++++ workhorse

for i_file, outfile in enumerate(outfile_list):
    outpath = os.path.join(outdir, outfile)
    print(f'working on {outfile}')

    # original catalogue
    ori_file = os.path.join(indir, orifile_list[i_file])
    with fits.open(ori_file) as hdul:
        simu_cata_ori = Table(hdul[1].data).to_pandas()
    print('number original', len(simu_cata_ori))

    # shape
    shape_file = os.path.join(indir, shapefile_list[i_file])
    simu_cata_shape = pd.read_feather(shape_file)
    print('number shape', len(simu_cata_shape))

    ## mask nan values and size
    mask_nan = (~simu_cata_shape['sersic_n'].isnull())
    mask_size = (simu_cata_shape['Re_arcsec']>Re_cut[0]) & (simu_cata_shape['Re_arcsec']<Re_cut[1])
    simu_cata_shape = simu_cata_shape[mask_nan & mask_size]
    simu_cata_shape.reset_index(drop=True, inplace=True)

    ## discrete sersic index
    sersic_n_discrete = np.array([np.logspace(np.log10(0.5), np.log10(6.0), 100)]).T
    ### build KDTree
    kdt = sst.cKDTree(sersic_n_discrete, leafsize=100)
    ### query
    bn_ori = np.array([simu_cata_shape['sersic_n']]).T
    dist, ind = kdt.query(bn_ori, k=1, distance_upper_bound=np.inf)
    bn_new = sersic_n_discrete[ind].flatten()
    simu_cata_shape.loc[:, 'shape/sersic_n'] = bn_new

    ## delete duplicated columns
    for col in simu_cata_shape:
        if (col in simu_cata_ori) and (col != 'index'):
            simu_cata_shape.drop(columns=[col], inplace=True)

    # merge
    simu_cata_ori = simu_cata_ori.merge(simu_cata_shape, how='inner', on='index')
    print('number final', len(simu_cata_ori))

    # correct magnitudes
    if col_Upsilon_corr is not None:
        for col_mag in cols_used_mags:
            simu_cata_ori.loc[:, f'{col_mag}_corr'] = np.array(simu_cata_ori[col_mag]) + np.array(simu_cata_ori[col_Upsilon_corr])

    # desired columns
    if col_Upsilon_corr is not None:
        cols_used = cols_used_others + [f'{col_mag}_corr' for col_mag in cols_used_mags] + [col_Upsilon_corr]
    else:
        cols_used = cols_used_others + cols_used_mags
    simu_cata_ori = simu_cata_ori[cols_used]

    # save
    Table.from_pandas(simu_cata_ori).write(outpath, format='fits')
    # simu_cata_ori.to_feather(outfile)
    print(f'final data saved as {outfile}')

# working on SS_CUM_r_Ks_patch0.fits
# number original 19365969
# number shape 8794876
# number final 8794876
# final data saved as SS_CUM_r_Ks_patch0.fits
# working on SS_CUM_r_Ks_patch1.fits
# number original 19116707
# number shape 8715550
# number final 8715550
# final data saved as SS_CUM_r_Ks_patch1.fits
# working on SS_CUM_r_Ks_patch2.fits
# number original 19039242
# number shape 8686136
# number final 8686136
# final data saved as SS_CUM_r_Ks_patch2.fits
# working on SS_CUM_r_Ks_patch3.fits
# number original 19371513
# number shape 8829566
# number final 8829566
# final data saved as SS_CUM_r_Ks_patch3.fits
