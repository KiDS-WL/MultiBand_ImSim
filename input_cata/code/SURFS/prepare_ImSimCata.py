# @Author: lshuns
# @Date:   2021-02-19, 16:19:29
# @Last modified by:   lshuns
# @Last modified time: 2021-02-21, 15:02:00

### to construct a catalogue as ImSim input

import numpy as np
import pandas as pd
import scipy.spatial as sst
import matplotlib.pyplot as plt

from scipy import stats
from astropy.io import fits

# overall cut
mag_r_max_cut = 27.
z_max_cut = 2.5
Re_cut = [1e-2, 10.] # arcsec

# ++++++++++++ Files

cata_dir = '/disks/shear10/ssli/ImSim/input/SURFS_cata/'

## +++ output
outpath_9_photo = cata_dir + 'cosmos_shape_z_ugriZYJHKs_selected.feather'

## +++ input
# original
ori_file = cata_dir + 'SURFS_SHARK_LC_4sqdeg_testcols.fits'
with fits.open(ori_file) as hdul:
    simu_cata = hdul[1].data
## general selection
mask_mag = (simu_cata['r_SDSS_apparent']<mag_r_max_cut)
mask_z = (simu_cata['zobs']<z_max_cut)
mask_tmp = mask_mag & mask_z
simu_cata = simu_cata[mask_tmp]
## useful columns
simu_ori_cata = pd.DataFrame({'index': np.arange(len(simu_cata)).astype(int),
                                'zobs': np.array(simu_cata['zobs']).astype(float),
                                'ra': np.array(simu_cata['ra']).astype(float),
                                'dec': np.array(simu_cata['dec']).astype(float),
                                'u_SDSS_apparent': np.array(simu_cata['u_SDSS_apparent']).astype(float),
                                'g_SDSS_apparent': np.array(simu_cata['g_SDSS_apparent']).astype(float),
                                'r_SDSS_apparent': np.array(simu_cata['r_SDSS_apparent']).astype(float),
                                'i_SDSS_apparent': np.array(simu_cata['i_SDSS_apparent']).astype(float),
                                'z_SDSS_apparent': np.array(simu_cata['z_SDSS_apparent']).astype(float),
                                'Y_VISTA_apparent': np.array(simu_cata['Y_VISTA_apparent']).astype(float),
                                'J_VISTA_apparent': np.array(simu_cata['J_VISTA_apparent']).astype(float),
                                'H_VISTA_apparent': np.array(simu_cata['H_VISTA_apparent']).astype(float),
                                'K_VISTA_apparent': np.array(simu_cata['K_VISTA_apparent']).astype(float)})
print('number original', len(simu_ori_cata))

# shape
shape_file = cata_dir + 'foursqdeg_testcols_COSMOS_morphology.feather'
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
#
# number original 881394
# number final 880587
