# @Author: lshuns
# @Date:   2021-02-22, 15:13:40
# @Last modified by:   lshuns
# @Last modified time: 2021-02-23, 12:37:43

### to construct a catalogue for ImSim

import numpy as np
import pandas as pd
import scipy.spatial as sst
import matplotlib.pyplot as plt

from scipy import stats
from astropy.io import fits
from astropy.table import Table

# overall cut
mag_r_max_cut = 27.
z_max_cut = 2.5

cata_dir = '/disks/shear10/ssli/ImSim/input/Flagship_cata/'

## +++ output
outpath_9_photo = cata_dir + 'v1_10_11_ra220_222_dec20_22_selected.feather'

## +++ input
# original
ori_file = cata_dir + 'v1_10_11_ra220_222_dec20_22.fits'
with fits.open(ori_file) as hdul:
    data = hdul[1].data
## general selection
mask_mag = (data['kids_r']<mag_r_max_cut)
mask_z = (data['observed_redshift_gal']<z_max_cut)
mask_tmp = mask_mag & mask_z
data = data[mask_tmp]

# ++ output
## unique index
index = data['halo_id']*10000 + data['galaxy_id']
data_out = pd.DataFrame({'index': (data['halo_id']*10000 + data['galaxy_id']).astype(int)})
# position
incol = ['ra_mag_gal', 'dec_mag_gal', 'observed_redshift_gal']
outcol = ['position/ra/obs', 'position/dec/obs', 'position/z/obs']
for i_p in range(len(incol)):
    data_out[outcol[i_p]] = data[incol[i_p]]
# magnitudes
incol = ['kids_u', 'kids_g', 'kids_r', 'kids_i', 'lsst_z', 'lsst_y', '2mass_j', '2mass_h', '2mass_ks']
outcol = ['mags/lensed/u', 'mags/lensed/g', 'mags/lensed/r', 'mags/lensed/i', 'mags/lensed/Z', 'mags/lensed/Y', 'mags/lensed/J', 'mags/lensed/H', 'mags/lensed/Ks']
d_mu = 2.0 * data['kappa']
for i_band in range(len(incol)):
    data_out[outcol[i_band]] = ((-2.5 * np.log10(data[incol[i_band]]) - 48.6) - 2.5 * np.log10(1 + d_mu)).astype(float)
# shapes
incol = ['disk_angle', 'bulge_fraction', 'bulge_r50', 'bulge_axis_ratio', 'disk_r50', 'disk_axis_ratio']
outcol = ['shape/position_angle', 'shape/bulge/fraction', 'shape/bulge/size', 'shape/bulge/axis_ratio', 'shape/disk/size', 'shape/disk/axis_ratio']
for i_shape in range(len(incol)):
    data_out[outcol[i_shape]] = data[incol[i_shape]]

#++++++++ make sersic index discrete to speed the ImSim
sersic_n_discrete = np.array([np.exp(np.linspace(np.log(0.5), np.log(6.0), 100))]).T
# build KDTree
kdt = sst.cKDTree(sersic_n_discrete, leafsize=100)
# query
bn_ori = np.array([data['bulge_nsersic']]).T
dist, ind = kdt.query(bn_ori, k=1, distance_upper_bound=np.inf)
bn_new = sersic_n_discrete[ind].flatten()
data_out['shape/bulge/sersic_n'] = bn_new

# ++++++++ only mag_r < 25 really matter
data_out.loc[data_out['mags/lensed/r']>=25., 'shape/bulge/sersic_n'] = 4.

# ++++++++ save
data_out.to_feather(outpath_9_photo)
