# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-10-19 13:15:31
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-20 22:28:48

### Diagnostic plot after alphaRecal step1:
###### alpha in measurement variance as a function of R and SNR

import os
import argparse

import numpy as np
import pandas as pd 

import statsmodels.api as sm

import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rc('font', size=5)
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import AutoMinorLocator, LogLocator, NullFormatter, NullLocator

# +++++++++++++++++++++++++++++ parser for command-line interfaces
parser = argparse.ArgumentParser(
    description=f"QC_step1.py: make diagnostic plot for alphaRecal step1",
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--inpath", type=str,
    help="the catalogue after alphaRecal step1\n\
    supported formats: feather")
parser.add_argument(
    "--outpath", type=str, default='./QC_alphaRecal_step1.png',
    help="the output plot name")
parser.add_argument(
    "--col_R", type=str,
    help="column name to the R para in the catalogue")
parser.add_argument(
    "--col_snr", type=str,
    help="column name to the SNR para in the catalogue")
parser.add_argument(
    "--cols_e12_raw", type=str, nargs=2, default=None,
    help="column names for e1, e2, raw measurement before self-recal")
parser.add_argument(
    "--cols_psf_e12", type=str, nargs=2, 
    help="column names for e1_psf, e2_psf")
parser.add_argument(
    "--col_var_before", type=str,
    help="column name to the variance before alphaRecal")
parser.add_argument(
    "--col_var_after", type=str,
    help="column name to the variance after alphaRecal")
parser.add_argument(
    "--Nbins_R_snr", type=int, nargs=2, default=[30, 30],
    help="Number of bins for the check plot, in order of R and SNR")

## arg parser
args = parser.parse_args()
inpath = args.inpath
outpath = args.outpath

col_R = args.col_R
col_snr = args.col_snr
col_e1_raw, col_e2_raw = args.cols_e12_raw
col_psf_e1, col_psf_e2 = args.cols_psf_e12
col_var_before = args.col_var_before
col_var_after = args.col_var_after
N_R, N_SNR = args.Nbins_R_snr
del parser, args

# >>>>>>>>>>>>>>>>>>>>>> general info

# plot range
RANGE_alpha = [-0.35, 0.35]
RANGE_R = [0.1, 1]
RANGE_SNR = [10**0.5, 10**2.5]
## colour bar range in the 2D map
cmin = -0.4
cmax = 0.4

# >>>>>>>>>>>>>>>>>>>>>> load catalogue
cata = pd.read_feather(inpath)
## used columns
cata = cata[[col_R, col_snr, 
                col_e1_raw, col_e2_raw, 
                col_psf_e1, col_psf_e2,
                col_var_before, col_var_after]]

# >>>>>>>>>>>>>>>>>>>>>> calculate alpha as a function of R and SNR

# start with R binn
cata.loc[:, 'bin_R'] = pd.qcut(cata[col_R], N_R, 
                                    labels=False, retbins=False)

# initialize
cata.loc[:, 'bin_snr'] = -999

# in each R bin, do SNR binning
for ibin_R in range(N_R):

    # select catalogue
    mask_binR = (cata['bin_R'].values == ibin_R)

    # bin in R
    cata.loc[mask_binR, 'bin_snr'] = pd.qcut(cata.loc[mask_binR, col_snr].values, N_SNR, 
                                                labels=False, retbins=False)

# print(">>>>>>>>>>>>>>> bin R", np.sort(cata['bin_R']))
# print(">>>>>>>>>>>>>>> bin SNR", np.sort(cata['bin_snr']))

# group based on binning
cata_grouped = cata.groupby(by=['bin_R', 'bin_snr'])
Ngroups = cata_grouped.ngroups
# print('>>> number of groups', Ngroups)
del cata

# loop over bins and calcuate alpha
## wanted info
R_median = np.zeros(Ngroups)
SNR_median = np.zeros(Ngroups)
alpha_before = np.zeros(Ngroups)
alpha_before_err = np.zeros(Ngroups)
alpha_after = np.zeros(Ngroups)
alpha_after_err = np.zeros(Ngroups)

i_group = 0
for name, group in cata_grouped:

    # save median center
    R_median[i_group] = np.median(np.array(group[col_R]))
    SNR_median[i_group] = np.median(np.array(group[col_snr]))

    # out shear
    e1_out = np.array(group[col_e1_raw])
    e2_out = np.array(group[col_e2_raw])
    emod_out = np.hypot(e1_out, e2_out)

    # out PSF 
    e1_psf = np.array(group[col_psf_e1])
    e2_psf = np.array(group[col_psf_e2])

    # rotate to PSF frame
    mask_tmp = emod_out > 0
    e_psf_rotated = e1_out * e1_psf + e2_out * e2_psf
    e_psf_rotated[mask_tmp] /= emod_out[mask_tmp]
    del mask_tmp, e1_psf, e2_psf, e1_out, e2_out, emod_out

    # calculate alpha using least square
    ### before correction
    Var = np.array(group[col_var_before])
    mod_ols = sm.OLS(Var, sm.add_constant(e_psf_rotated))
    res_ols = mod_ols.fit()
    alpha_before[i_group] = res_ols.params[1]
    alpha_before_err[i_group] = (res_ols.cov_params()[1, 1])**0.5
    del Var
    ### after correction
    Var = np.array(group[col_var_after])
    mod_ols = sm.OLS(Var, sm.add_constant(e_psf_rotated))
    res_ols = mod_ols.fit()
    alpha_after[i_group] = res_ols.params[1]
    alpha_after_err[i_group] = (res_ols.cov_params()[1, 1])**0.5
    del group, Var, e_psf_rotated

    i_group += 1

del cata_grouped

# >>>>>>>>>>>>>>>>>>>>>> plot
fig, axs = plt.subplots(2, 3)
fig.suptitle(os.path.basename(inpath))

####### before correction

TITLE = 'before alphaRecal'

# 1. in 2D
norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
axs[0, 0].scatter(R_median, SNR_median, c=alpha_before, 
                s=2, cmap='coolwarm', norm=norm)

SM = ScalarMappable(norm=norm, cmap='coolwarm')
SM.set_array([])
plt.colorbar(SM, ax=axs[0, 0])
axs[0, 0].set_xlim(RANGE_R[0], RANGE_R[1])
axs[0, 0].set_ylim(RANGE_SNR[0], RANGE_SNR[1])
axs[0, 0].set_yscale('log')
axs[0, 0].xaxis.set_minor_locator(AutoMinorLocator())
axs[0, 0].yaxis.set_minor_locator(LogLocator(base=10.0, subs=None, numticks=10))

axs[0, 0].set_xlabel('median R')
axs[0, 0].set_ylabel('median SNR')
axs[0, 0].set_title(TITLE)

# 2. as a function of SNR
axs[0, 1].errorbar(SNR_median, alpha_before, 
                    yerr = np.vstack([alpha_before_err, alpha_before_err]), 
                    linestyle='', marker='o', markersize=1, elinewidth=0.6, capsize=1)

axs[0, 1].axhline(y=0, ls='--', color='k', linewidth=0.5)
axs[0, 1].set_xlim(RANGE_SNR[0], RANGE_SNR[1])
axs[0, 1].set_ylim(RANGE_alpha[0], RANGE_alpha[1])
axs[0, 1].set_xscale('log')
axs[0, 1].set_xlabel('median SNR')
axs[0, 1].set_ylabel('alphaS in variance')
axs[0, 1].set_title(TITLE)

# 3. as a function of R
axs[0, 2].errorbar(R_median, alpha_before, 
                    yerr = np.vstack([alpha_before_err, alpha_before_err]), 
                    linestyle='', marker='o', markersize=1, elinewidth=0.6, capsize=1)

axs[0, 2].axhline(y=0, ls='--', color='k', linewidth=0.5)
axs[0, 2].set_xlim(RANGE_R[0], RANGE_R[1])
axs[0, 2].set_ylim(RANGE_alpha[0], RANGE_alpha[1])
axs[0, 2].set_xlabel('median R')
axs[0, 2].set_ylabel('alphaS in variance')
axs[0, 2].set_title(TITLE)

####### after correction

TITLE = 'after alphaRecal'

# 1. in 2D
norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
axs[1, 0].scatter(R_median, SNR_median, c=alpha_after, 
                s=2, cmap='coolwarm', norm=norm)

SM = ScalarMappable(norm=norm, cmap='coolwarm')
SM.set_array([])
plt.colorbar(SM, ax=axs[1, 0])
axs[1, 0].set_xlim(RANGE_R[0], RANGE_R[1])
axs[1, 0].set_ylim(RANGE_SNR[0], RANGE_SNR[1])
axs[1, 0].set_yscale('log')
axs[1, 0].xaxis.set_minor_locator(AutoMinorLocator())
axs[1, 0].yaxis.set_minor_locator(LogLocator(base=10.0, subs=None, numticks=10))

axs[1, 0].set_xlabel('median R')
axs[1, 0].set_ylabel('median SNR')
axs[1, 0].set_title(TITLE)

# 2. as a function of SNR
axs[1, 1].errorbar(SNR_median, alpha_after, 
                    yerr = np.vstack([alpha_after_err, alpha_after_err]), 
                    linestyle='', marker='o', markersize=1, elinewidth=0.6, capsize=1)

axs[1, 1].axhline(y=0, ls='--', color='k', linewidth=0.5)
axs[1, 1].set_xlim(RANGE_SNR[0], RANGE_SNR[1])
axs[1, 1].set_ylim(RANGE_alpha[0], RANGE_alpha[1])
axs[1, 1].set_xscale('log')
axs[1, 1].set_xlabel('median SNR')
axs[1, 1].set_ylabel('alphaS in variance')
axs[1, 1].set_title(TITLE)

# 3. as a function of R
axs[1, 2].errorbar(R_median, alpha_after, 
                    yerr = np.vstack([alpha_after_err, alpha_after_err]), 
                    linestyle='', marker='o', markersize=1, elinewidth=0.6, capsize=1)

axs[1, 2].axhline(y=0, ls='--', color='k', linewidth=0.5)
axs[1, 2].set_xlim(RANGE_R[0], RANGE_R[1])
axs[1, 2].set_ylim(RANGE_alpha[0], RANGE_alpha[1])
axs[1, 2].set_xlabel('median R')
axs[1, 2].set_ylabel('alphaS in variance')
axs[1, 2].set_title(TITLE)

##### save plot
plt.tight_layout()
plt.savefig(outpath, dpi=300)
plt.close()
print("QC plot for alphaRecal step 1 saved as", outpath)

