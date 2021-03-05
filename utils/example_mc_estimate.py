# @Author: lshuns
# @Date:   2021-03-03, 18:19:42
# @Last modified by:   lshuns
# @Last modified time: 2021-03-04, 23:19:02

### an example to calculate the mc shear bias from simulated catalogues

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

import os
import sys
# Self-defined package
sys.path.insert(0,os.path.realpath('../'))
from modules import ShearBias

# +++++++++++++++++++++++++++++ the only variables you want to modify

# ============  general
# Bootstrap for errors
## set nboot=0, if that is not needed
nboot = 50
rng_seed_boot = 201294

# outpath
## set to 'show', if you only want to see the results on the screen
## otherwise, some file names for the txt output
# out_path = 'show'
out_path = './test/surfs_1psf.txt'


# ============ simulation catalogue info
# the input simulation catalogue
### supported formats: fits, csv, feather
###             guess from the file suffix
in_file_sim = '/disks/shear15/ssli/ImSim/output/test_surfs_onePFSmean_combined.fits'

# the name of the mask column
### 0 for good
### comment it out, if mask is not required
col_mask = 'MASK_gold'

# column name for redshift bins
### comment it out, if tomographic binnning is not required
### the output zbin_id specify bins with 0 corresponding to all samples
col_z_sim = 'Z_B'

# tomographic bin edges
### following KiDS fashion: (low, high]
### not used, if col_z is not set
z_bin_edges = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]

# column names for the desired parameters
### order: input shear (g1_in, g2_in),
###         measured e (e1_out, e2_out), measured shape (size_out),
###         SNR for shape measurement (shape_snr), weights for shape measurement (shape_weight)
###         PSF quadrupole (psf_Q11, psf_Q22, psf_Q12) or PSF size (psf_size) [depends on psf_sizeORmoments_sim below]
cols_sim = ['g1_in', 'g2_in',
            'e1_LF_r', 'e2_LF_r', 'scalelength_LF_r',
            'SNR_LF_r', 'weight_LF_r',
            'psf_Q11_LF_r', 'psf_Q22_LF_r', 'psf_Q12_LF_r']

# use PSF moments or PSF size ?
# psf_sizeORmoments_sim = 'size'
psf_sizeORmoments_sim = 'moments'

# ============ observation catalogue info
# the input observation catalogue
### comment it out, if reweighting is not required
### if provided, make sure galaxies have been selected
in_file_obs = '/disks/shear15/ssli/KV450/combined/KV450_reweight_3x4x4_v2_good_combined_selected.feather'

# column names for the desired parameters
### order: measured e (e1, e2), measured shape (size),
###         SNR for shape measurement (shape_snr), weights for shape measurement (shape_weight)
###         PSF quadrupole (psf_Q11, psf_Q22, psf_Q12) or PSF size (psf_size) [depends on psf_sizeORmoments_obs below]
cols_obs = ['bias_corrected_e1', 'bias_corrected_e2', 'bias_corrected_scalelength_pixels',
            'model_SNratio', 'recal_weight',
            'PSF_Q11', 'PSF_Q22', 'PSF_Q12']

# column name for redshift bins
### comment it out, if tomographic binnning is not required
### the output zbin_id specify bins with 0 corresponding to all samples
col_z_obs = 'Z_B'

# use PSF moments or PSF size ?
# psf_sizeORmoments_obs = 'size'
psf_sizeORmoments_obs = 'moments'

# +++++++++++++++++++++++++++++ workhorse

# load simulation catalogue
file_type = in_file_sim[-3:]
if file_type == 'csv':
    cata_sim = pd.read_csv(in_file_sim)
elif file_type == 'her':
    cata_sim = pd.read_feather(in_file_sim)
elif file_type == 'its':
    with fits.open(in_file_sim) as hdul:
        cata_sim = Table(hdul[1].data).to_pandas()
else:
    raise Exception(f'Not supported input file type! {in_file_sim}')
print('Number of sources in simulated catalogue', len(cata_sim))

# selection
if 'col_mask' in locals():
    cata_sim = cata_sim[cata_sim[col_mask]==0]
    cata_sim.reset_index(drop=True, inplace=True)
    print('Number of selected sources in simulated catalogue', len(cata_sim))

# rename column names
cols_rename = {}
cols_goodnames = ['g1_in', 'g2_in',
            'e1_out', 'e2_out', 'size_out',
            'shape_snr', 'shape_weight']
if psf_sizeORmoments_sim == 'size':
    cols_goodnames.append('psf_size')
elif psf_sizeORmoments_sim == 'moments':
    cols_goodnames.extend(['psf_Q11', 'psf_Q22', 'psf_Q12'])
for i_col, col_name in enumerate(cols_sim):
     cols_rename[col_name] = cols_goodnames[i_col]
cata_sim.rename(columns=cols_rename, inplace=True)

# tomographic binning
cata_sim_list = [cata_sim]
if 'col_z_sim' in locals():
    N_zbins = len(z_bin_edges)-1
    print(f'Use {N_zbins} tomographic bins')
    for i_zbin in range(N_zbins):
        z_min_bin = z_bin_edges[i_zbin]
        z_max_bin = z_bin_edges[i_zbin+1]
        mask_zbin = (cata_sim[col_z_sim]>z_min_bin) & (cata_sim[col_z_sim]<=z_max_bin)
        cata_sim_selec = cata_sim[mask_zbin]
        cata_sim_selec.reset_index(drop=True, inplace=True)
        cata_sim_list.append(cata_sim_selec)

# output
if out_path == 'show':
    f = sys.stdout
else:
    f = open(out_path, 'w')

# simple without data re-weighting
if not 'in_file_obs' in locals():
    print('Shear bias estimation without data re-weighting')
    for i_cata, cata_sim_tmp in enumerate(cata_sim_list):
        res = ShearBias.mcCalFunc_simple(cata_sim_tmp, nboot=nboot, rng_seed_boot=rng_seed_boot)

        # print out some comments and column names
        if i_cata == 0:
            print('# zbin_id = 0 correspond to all samples, number increase towards high z bins', file=f)
            # collect columns names
            cols = '# ' + '    '.join(list(res.keys()) + ['zbin_id'])
            print(cols, file=f)

        # print out values
        vals = '    '.join(["{0:0.4f}".format(val) for val in res.values()]) + f'    {i_cata}'
        print(vals, file=f)

    # close what is opened
    if out_path != 'show':
        f.close()
    sys.exit()

# load observation catalogue
file_type = in_file_obs[-3:]
if file_type == 'csv':
    cata_obs = pd.read_csv(in_file_obs)
elif file_type == 'her':
    cata_obs = pd.read_feather(in_file_obs)
elif file_type == 'its':
    with fits.open(in_file_obs) as hdul:
        cata_obs = Table(hdul[1].data).to_pandas()
else:
    raise Exception(f'Not supported input file type! {in_file_obs}')
print('Number of sources in observation catalogue', len(cata_obs))

# rename column names
cols_rename = {}
cols_goodnames = ['e1', 'e2', 'size',
            'shape_snr', 'shape_weight']
if psf_sizeORmoments_obs == 'size':
    cols_goodnames.append('psf_size')
elif psf_sizeORmoments_obs == 'moments':
    cols_goodnames.extend(['psf_Q11', 'psf_Q22', 'psf_Q12'])
for i_col, col_name in enumerate(cols_obs):
     cols_rename[col_name] = cols_goodnames[i_col]
cata_obs.rename(columns=cols_rename, inplace=True)

# tomographic binning
cata_obs_list = [cata_obs]
if 'col_z_obs' in locals():
    N_zbins = len(z_bin_edges)-1
    for i_zbin in range(N_zbins):
        z_min_bin = z_bin_edges[i_zbin]
        z_max_bin = z_bin_edges[i_zbin+1]
        mask_zbin = (cata_obs[col_z_obs]>z_min_bin) & (cata_obs[col_z_obs]<=z_max_bin)
        cata_obs_selec = cata_obs[mask_zbin]
        cata_obs_selec.reset_index(drop=True, inplace=True)
        cata_obs_list.append(cata_obs_selec)

# calculation
print('Shear bias estimation with data re-weighting')
for i_cata, cata_sim_tmp in enumerate(cata_sim_list):
    cata_obs_tmp = cata_obs_list[i_cata]
    res = ShearBias.mcCalFunc_reweight_R_SNR(cata_sim_tmp, cata_obs_tmp, Nbin_SNR=20, Nbin_R=20, nboot=nboot, rng_seed_boot=rng_seed_boot)

    # print out some comments and column names
    if i_cata == 0:
        print('# zbin_id = 0 correspond to all samples, number increase towards high z bins', file=f)
        # collect columns names
        cols = '# ' + '    '.join(list(res.keys()) + ['zbin_id'])
        print(cols, file=f)

    # print out values
    vals = '    '.join(["{0:0.4f}".format(val) for val in res.values()]) + f'    {i_cata}'
    print(vals, file=f)

# close what is opened
if out_path != 'show':
    print(f'results saved to {out_path}')
    f.close()
