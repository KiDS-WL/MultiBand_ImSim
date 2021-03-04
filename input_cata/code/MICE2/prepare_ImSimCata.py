# @Author: ssli
# @Date:   2021-02-25, 18:11:41
# @Last modified by:   ssli
# @Last modified time: 2021-02-25, 18:18:34

### to construct a catalogue as ImSim input

import numpy as np
import pandas as pd
import scipy.spatial as sst
import matplotlib.pyplot as plt

from scipy import stats
from astropy.io import fits

# overall cut
Re_cut = [1e-2, 10.] # arcsec

# ++++++++++++ Files

cata_dir = '/disks/shear10/ssli/ImSim/input/MICE2_cata/'

## +++ output
outpath_9_photo = cata_dir + 'cosmos_shape_z_ugriZYJHKs_selected.feather'

## +++ input
# original
ori_file = cata_dir + 'MICE2_ra35-37_dec5-7_corrected.fits'
with fits.open(ori_file) as hdul:
    simu_cata = hdul[1].data
## useful columns
simu_ori_cata = pd.DataFrame({'index': np.array(simu_cata['unique_gal_id']).astype(int),
                                'zobs': np.array(simu_cata['z_cgal_v']).astype(float),
                                'ra': np.array(simu_cata['ra_gal_mag']).astype(float),
                                'dec': np.array(simu_cata['dec_gal_mag']).astype(float),
                                'u': np.array(simu_cata['sdss_u_true_lensed']).astype(float),
                                'g': np.array(simu_cata['sdss_g_true_lensed']).astype(float),
                                'r': np.array(simu_cata['sdss_r_true_lensed']).astype(float),
                                'i': np.array(simu_cata['sdss_i_true_lensed']).astype(float),
                                'Z': np.array(simu_cata['sdss_z_true_lensed']).astype(float),
                                'Y': np.array(simu_cata['des_asahi_full_y_true_lensed']).astype(float),
                                'J': np.array(simu_cata['vhs_j_true_lensed']).astype(float),
                                'H': np.array(simu_cata['vhs_h_true_lensed']).astype(float),
                                'Ks': np.array(simu_cata['vhs_ks_true_lensed']).astype(float)})
print('number original', len(simu_ori_cata))

# shape
shape_file = cata_dir + 'MICE2_ra35-37_dec5-7_COSMOS_morphology.feather'
simu_cata_shape = pd.read_feather(shape_file)
## mask nan values and size
mask_nan = (~simu_cata_shape['sersic_n_bf'].isnull())
mask_size_bf = (simu_cata_shape['Re_arcsec_bf']>Re_cut[0]) & (simu_cata_shape['Re_arcsec_bf']<Re_cut[1])
mask_size_c = (simu_cata_shape['Re_arcsec_c']>Re_cut[0]) & (simu_cata_shape['Re_arcsec_c']<Re_cut[1])
mask_size_r = (simu_cata_shape['Re_arcsec_r']>Re_cut[0]) & (simu_cata_shape['Re_arcsec_r']<Re_cut[1])
simu_cata_shape = simu_cata_shape[mask_nan & mask_size_bf & mask_size_c & mask_size_r].copy()
simu_cata_shape.reset_index(drop=True, inplace=True)
## discrete sersic index
sersic_n_discrete = np.array([np.exp(np.linspace(np.log(0.5), np.log(6.0), 100))]).T
# build KDTree
kdt = sst.cKDTree(sersic_n_discrete, leafsize=100)
# query
## bf
bn_ori = np.array([simu_cata_shape['sersic_n_bf']]).T
dist, ind = kdt.query(bn_ori, k=1, distance_upper_bound=np.inf)
bn_new = sersic_n_discrete[ind].flatten()
simu_cata_shape.loc[:, 'shape/sersic_n_bf'] = bn_new
## c
bn_ori = np.array([simu_cata_shape['sersic_n_c']]).T
dist, ind = kdt.query(bn_ori, k=1, distance_upper_bound=np.inf)
bn_new = sersic_n_discrete[ind].flatten()
simu_cata_shape.loc[:, 'shape/sersic_n_c'] = bn_new
## r
bn_ori = np.array([simu_cata_shape['sersic_n_r']]).T
dist, ind = kdt.query(bn_ori, k=1, distance_upper_bound=np.inf)
bn_new = sersic_n_discrete[ind].flatten()
simu_cata_shape.loc[:, 'shape/sersic_n_r'] = bn_new
## useful columns
cols = ['index',
            'Re_arcsec_bf', 'Re_arcsec_c', 'Re_arcsec_r',
            'BA_bf', 'BA_c', 'BA_r',
            'shape/sersic_n_bf', 'shape/sersic_n_c', 'shape/sersic_n_r']
simu_cata_shape = simu_cata_shape[cols].copy()

# ++++++++++++ combine & save
simu_cata_final = simu_ori_cata.merge(simu_cata_shape, how='right', on='index')
print('number final', len(simu_cata_final))
simu_cata_final.to_feather(outpath_9_photo)
