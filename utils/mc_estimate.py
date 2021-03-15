# @Author: lshuns
# @Date:   2021-03-03, 18:19:42
# @Last modified by:   lshuns
# @Last modified time: 2021-03-15, 22:30:20

### a script to calculate the mc shear bias from combined simulated catalogues

import argparse

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

import os
import sys
# Self-defined package
sys.path.insert(0,os.path.realpath('../'))
from modules import ShearBias

# +++++++++++++++++++++++++++++ parser for command-line interfaces
parser = argparse.ArgumentParser(
    description=f"mc_estimate.py: estimate mc shear bias from simulated catalogue.",
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--in_file_sim", type=str,
    help="the input simulation catalogue.\n\
    supported formats: fits, csv, feather\n\
        guess from the file suffix")
parser.add_argument(
    "--in_file_obs", type=str, default=None,
    help="the input observation catalogue.\n\
        If not set, data-reweighting will not be performed.\n\
        If provided, make sure galaxies have already been selected.")
parser.add_argument(
    "--out_path", type=str, default='show',
    help="where to save the final results. \n\
    If not provided, results will be shown on the screen.")
parser.add_argument(
    "--use_tomographic", action="store_true",
    help="calculate shear bias tomographically")
parser.add_argument(
    "--col_z_sim", type=str, default='Z_B',
    help="column name of photoz in simulation")

## arg parser
args = parser.parse_args()
in_file_sim = args.in_file_sim
in_file_obs = args.in_file_obs
out_path = args.out_path
use_tomographic = args.use_tomographic
col_z_sim = args.col_z_sim

# +++++++++++++++++++++++++++++ the only variables may you want to modify

# Bootstrap for errors
## set nboot=0, if that is not needed
nboot = 50
rng_seed_boot = 201294

# tomographic bin edges
### following KiDS fashion: (low, high]
### not used, if use_tomographic is not set
z_bin_edges = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]

# ============ simulation catalogue info
# column names for the desired parameters
### order: input shear (g1_in, g2_in),
###         measured e (e1_out, e2_out), measured shape (size_out),
###         SNR for shape measurement (shape_snr), weights for shape measurement (shape_weight)
###         PSF quadrupole (psf_Q11, psf_Q22, psf_Q12) or PSF size (psf_size) [depends on psf_sizeORmoments_sim below]
cols_sim = ['g1_in', 'g2_in',
            'e1_LF_r', 'e2_LF_r', 'scalelength_LF_r',
            'SNR_LF_r', 'weight_global_LF_r',
            'psf_Q11_LF_r', 'psf_Q22_LF_r', 'psf_Q12_LF_r']

# use PSF moments or PSF size ?
# psf_sizeORmoments_sim = 'size'
psf_sizeORmoments_sim = 'moments'

# ============ observation catalogue info
# column names for the desired parameters
### order: measured e (e1, e2), measured shape (size),
###         SNR for shape measurement (shape_snr), weights for shape measurement (shape_weight)
###         PSF quadrupole (psf_Q11, psf_Q22, psf_Q12) or PSF size (psf_size) [depends on psf_sizeORmoments_obs below]
cols_obs = ['bias_corrected_e1', 'bias_corrected_e2', 'bias_corrected_scalelength_pixels',
            'model_SNratio', 'recal_weight',
            'PSF_Q11', 'PSF_Q22', 'PSF_Q12']

# column name for redshift bins
### the output zbin_id specify bins with 0 corresponding to all samples
### not used, if use_tomographic is not set
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
# lensfit cuts with 0 means no issue and -9 means large galaxies
fitcuts = ((cata_sim['class_LF_r']==0) | (cata_sim['class_LF_r']==-9))
# # remove invalid measurements
weight_cuts = (cata_sim['weight_global_LF_r']>0)
# remove too good sources
# snr_cuts = (cata_sim['SNR_LF_r']<210)
# remove potentially blended sources
blend_cuts = (cata_sim['contamination_radius_LF_r']>4.25)
# remove unresolved binary stars
binary_star_cuts = ((np.hypot(cata_sim['e1_LF_r'], cata_sim['e2_LF_r'])<=0.8) | (cata_sim['scalelength_LF_r']>=0.5*np.exp(0.65788*(24.2-cata_sim['MAG_AUTO']))))
# remove stars
star_cuts = (cata_sim['perfect_flag_star']==0)
# cata_sim = cata_sim[fitcuts & weight_cuts & snr_cuts & blend_cuts & binary_star_cuts & star_cuts]
cata_sim = cata_sim[fitcuts & weight_cuts & blend_cuts & binary_star_cuts & star_cuts]
cata_sim.reset_index(drop=True, inplace=True)
print('Number of selected sources in simulated catalogue', len(cata_sim))

# extract useable parameters
if use_tomographic:
    cata_sim = cata_sim[cols_sim+[col_z_sim]]
else:
    cata_sim = cata_sim[cols_sim]

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

# output
if out_path == 'show':
    f = sys.stdout
else:
    f = open(out_path, 'w')

# simple without data re-weighting
if in_file_obs is None:
    print('Shear bias estimation without data re-weighting')

    # whole sample
    res = ShearBias.mcCalFunc_simple(cata_sim, nboot=nboot, rng_seed_boot=rng_seed_boot)

    # print out some comments and column names
    print('# zbin_id = 0 correspond to all samples, number increase towards high z bins', file=f)
    # collect columns names
    cols = '# ' + '    '.join(list(res.keys()) + ['zbin_id'])
    print(cols, file=f)

    # print out values
    vals = '    '.join(["{0:0.4f}".format(val) for val in res.values()]) + '    0'
    print(vals, file=f)

    # print out mean values
    m = (res['m1']+res['m2'])/2.
    merr = (res['m1_err'] + res['m2_err'])/2.
    merror_low = (res['m1_err_BS_16'] + res['m2_err_BS_16'])/2.
    merror_high = (res['m1_err_BS_84'] + res['m2_err_BS_84'])/2.
    print(f'### whole: m = {m:.4f}, merr = {merr:.4f}, merr_BS_16 = {merror_low:.4f}, merr_BS_84 = {merror_high:.4f}')

    # tomographic binning
    if use_tomographic:
        N_zbins = len(z_bin_edges)-1
        print(f'Use {N_zbins} tomographic bins')
        for i_zbin in range(N_zbins):
            z_min_bin = z_bin_edges[i_zbin]
            z_max_bin = z_bin_edges[i_zbin+1]
            mask_zbin = (cata_sim[col_z_sim]>z_min_bin) & (cata_sim[col_z_sim]<=z_max_bin)
            cata_sim_selec = cata_sim[mask_zbin]
            cata_sim_selec.reset_index(drop=True, inplace=True)

            res = ShearBias.mcCalFunc_simple(cata_sim_selec, nboot=nboot, rng_seed_boot=rng_seed_boot)

            # print out values
            vals = '    '.join(["{0:0.4f}".format(val) for val in res.values()]) + f'    {i_zbin+1}'
            print(vals, file=f)

            # print out source number info
            print(f'### bin {i_zbin+1}: N = {len(cata_sim_selec)}, frac = {len(cata_sim_selec)/len(cata_sim)}')

            # print out mean values
            m = (res['m1']+res['m2'])/2.
            merr = (res['m1_err'] + res['m2_err'])/2.
            merror_low = (res['m1_err_BS_16'] + res['m2_err_BS_16'])/2.
            merror_high = (res['m1_err_BS_84'] + res['m2_err_BS_84'])/2.
            print(f'### m = {m:.4f}, merr = {merr:.4f}, merr_BS_16 = {merror_low:.4f}, merr_BS_84 = {merror_high:.4f}')

    # close what is opened
    if out_path != 'show':
        print(f'results saved to {out_path}')
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
# extract useable parameters
if use_tomographic:
    cata_obs = cata_obs[cols_obs+[col_z_obs]]
else:
    cata_obs = cata_obs[cols_obs]

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

# calculation
print('Shear bias estimation with data re-weighting')

# whole sample
res = ShearBias.mcCalFunc_reweight_R_SNR(cata_sim, cata_obs, Nbin_SNR=20, Nbin_R=20, nboot=nboot, rng_seed_boot=rng_seed_boot)

# print out some comments and column names
print('# zbin_id = 0 correspond to all samples, number increase towards high z bins', file=f)
# collect columns names
cols = '# ' + '    '.join(list(res.keys()) + ['zbin_id'])
print(cols, file=f)

# print out values
vals = '    '.join(["{0:0.4f}".format(val) for val in res.values()]) + '    0'
print(vals, file=f)

# print out mean values
m = (res['m1']+res['m2'])/2.
merr = (res['m1_err'] + res['m2_err'])/2.
merror_low = (res['m1_err_BS_16'] + res['m2_err_BS_16'])/2.
merror_high = (res['m1_err_BS_84'] + res['m2_err_BS_84'])/2.
print(f'### whole: m = {m:.4f}, merr = {merr:.4f}, merr_BS_16 = {merror_low:.4f}, merr_BS_84 = {merror_high:.4f}')

# tomographic binning
if use_tomographic:
    N_zbins = len(z_bin_edges)-1
    print(f'Use {N_zbins} tomographic bins')
    for i_zbin in range(N_zbins):
        z_min_bin = z_bin_edges[i_zbin]
        z_max_bin = z_bin_edges[i_zbin+1]
        mask_zbin = (cata_sim[col_z_sim]>z_min_bin) & (cata_sim[col_z_sim]<=z_max_bin)
        cata_sim_selec = cata_sim[mask_zbin]
        cata_sim_selec.reset_index(drop=True, inplace=True)

        mask_zbin = (cata_obs[col_z_obs]>z_min_bin) & (cata_obs[col_z_obs]<=z_max_bin)
        cata_obs_selec = cata_obs[mask_zbin]
        cata_obs_selec.reset_index(drop=True, inplace=True)

        res = ShearBias.mcCalFunc_reweight_R_SNR(cata_sim_selec, cata_obs_selec, Nbin_SNR=20, Nbin_R=20, nboot=nboot, rng_seed_boot=rng_seed_boot)

        # print out values
        vals = '    '.join(["{0:0.4f}".format(val) for val in res.values()]) + f'    {i_zbin+1}'
        print(vals, file=f)

        # print out source number info
        print(f'### bin {i_zbin+1}: N_sim = {len(cata_sim_selec)}, frac_sim = {len(cata_sim_selec)/len(cata_sim)}')
        print(f'                    N_obs = {len(cata_obs_selec)}, frac_obs = {len(cata_obs_selec)/len(cata_obs)}')

        # print out mean values
        m = (res['m1']+res['m2'])/2.
        merr = (res['m1_err'] + res['m2_err'])/2.
        merror_low = (res['m1_err_BS_16'] + res['m2_err_BS_16'])/2.
        merror_high = (res['m1_err_BS_84'] + res['m2_err_BS_84'])/2.
        print(f'### m = {m:.4f}, merr = {merr:.4f}, merr_BS_16 = {merror_low:.4f}, merr_BS_84 = {merror_high:.4f}')

# close what is opened
if out_path != 'show':
    print(f'results saved to {out_path}')
    f.close()
