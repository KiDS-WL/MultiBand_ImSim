# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-10-20 18:41:07
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-20 22:28:24

### Diagnostic plot after alphaRecal step2:
###### alpha in ellipticity as a function of R and SNR

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
    description=f"QC_step2.py: make diagnostic plot for alphaRecal step2",
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--inpath", type=str,
    help="the catalogue after alphaRecal step1\n\
    supported formats: feather")
parser.add_argument(
    "--outprefix", type=str, default='./QC_alphaRecal_step2',
    help="the output plot name")
parser.add_argument(
    "--col_ZB", type=str,
    help="column name to the photo-z para in the catalogue")
parser.add_argument(
    "--col_R", type=str,
    help="column name to the R para in the catalogue")
parser.add_argument(
    "--col_snr", type=str,
    help="column name to the SNR para in the catalogue")
parser.add_argument(
    "--cols_psf_e12", type=str, nargs=2, 
    help="column names for e1_psf, e2_psf")
parser.add_argument(
    "--cols_e12_before", type=str, nargs=2, default=None,
    help="column names for e1, e2, measurement before alphaRecal")
parser.add_argument(
    "--cols_e12_after", type=str, nargs=2, default=None,
    help="column names for e1, e2, measurement after alphaRecal")
parser.add_argument(
    "--col_weight_after", type=str,
    help="column name to the LF weight after alphaRecal")
parser.add_argument(
    "--Nbins_R_snr", type=int, nargs=2, default=[30, 30],
    help="Number of bins for the check plot, in order of R and SNR")
parser.add_argument(
    "--Z_B_edges", type=float, nargs='*', default=None,
    help="edges for tomographic binning")

## arg parser
args = parser.parse_args()
inpath = args.inpath
outprefix = args.outprefix

col_ZB = args.col_ZB
col_R = args.col_R
col_snr = args.col_snr
col_psf_e1, col_psf_e2 = args.cols_psf_e12
col_e1_before, col_e2_before = args.cols_e12_before
col_e1_after, col_e2_after = args.cols_e12_after
col_weight_after = args.col_weight_after

N_R, N_SNR = args.Nbins_R_snr
Z_B_edges = args.Z_B_edges
del parser, args

# >>>>>>>>>>>>>>>>>>>>>> general info

# plot range
RANGE_alpha = [-1.5, 1.0]
RANGE_R = [0.1, 1]
RANGE_SNR = [10**0.5, 10**2.5]
## colour bar range in the 2D map
cmin = -0.4
cmax = 0.4

# in ZB bins
RANGE_alpha_ZB = [0.2, -0.2]

############## as a function of SNR and R
# >>>>>>>>>>>>>>>>>>>>>> load catalogue
cata = pd.read_feather(inpath)
## used columns
cata = cata[[col_R, col_snr, 
            col_e1_before, col_e2_before, 
            col_e1_after, col_e2_after,
            col_psf_e1, col_psf_e2,
            col_weight_after
            ]]

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

# group based on binning
cata_grouped = cata.groupby(by=['bin_R', 'bin_snr'])
Ngroups = cata_grouped.ngroups
del cata

# loop over bins and calcuate alpha
## wanted info
alphaMap_df = pd.DataFrame(0, 
                index=np.arange(Ngroups), 
                columns=['R', 'SNR', 
                'alpha1_before', 'alpha1_before_err', 
                'alpha1_after', 'alpha1_after_err',
                'alpha2_before', 'alpha2_before_err',
                'alpha2_after', 'alpha2_after_err'])

i_group = 0
for name, group in cata_grouped:

    # save median center
    alphaMap_df.loc[i_group, 'R'] = np.median(np.array(group[col_R]))
    alphaMap_df.loc[i_group, 'SNR'] = np.median(np.array(group[col_snr]))

    # out PSF 
    e1_psf = np.array(group[col_psf_e1])
    e2_psf = np.array(group[col_psf_e2])

    # LF weights
    weight_out = np.array(group[col_weight_after])
    # out shear
    e1_out_list = [np.array(group[col_e1_before]), np.array(group[col_e1_after])] 
    e2_out_list = [np.array(group[col_e2_before]), np.array(group[col_e2_after])]
    col_save_list = ['before', 'after']
    del group

    # >>>>>>>>>>>>>>>> calculate alpha
    for i_e, e1_out in enumerate(e1_out_list):
        e2_out = e2_out_list[i_e]
        col_save = col_save_list[i_e]
    
        # calculate alpha using least square
        ## e1
        mod_wls = sm.WLS(e1_out, sm.add_constant(e1_psf), weights=weight_out)
        res_wls = mod_wls.fit()
        alphaMap_df.loc[i_group, f'alpha1_{col_save}'] = res_wls.params[1]
        alphaMap_df.loc[i_group, f'alpha1_{col_save}_err'] = (res_wls.cov_params()[1, 1])**0.5
        del res_wls, mod_wls, e1_out
        ## e2
        mod_wls = sm.WLS(e2_out, sm.add_constant(e2_psf), weights=weight_out)
        res_wls = mod_wls.fit()
        alphaMap_df.loc[i_group, f'alpha2_{col_save}'] = res_wls.params[1]
        alphaMap_df.loc[i_group, f'alpha2_{col_save}_err'] = (res_wls.cov_params()[1, 1])**0.5
        del res_wls, mod_wls, e2_out

    del weight_out, e1_psf, e2_psf, e1_out_list, e2_out_list

    i_group += 1

del cata_grouped

# >>>>>>>>>>>>>>>>>>>>>> plot

for i_e in range(1, 3):

    fig, axs = plt.subplots(2, 3)
    fig.suptitle(os.path.basename(inpath))
    outpath = outprefix + f'_SNR_R_alpha{i_e}.png'

    ####### before correction
    TITLE = 'before alphaRecal'

    # 1. in 2D
    norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
    axs[0, 0].scatter(alphaMap_df['R'].values, alphaMap_df['SNR'].values, c=alphaMap_df[f'alpha{i_e}_before'].values, 
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
    axs[0, 1].errorbar(alphaMap_df['SNR'].values, alphaMap_df[f'alpha{i_e}_before'].values,  
                        yerr = np.vstack([alphaMap_df[f'alpha{i_e}_before_err'].values, alphaMap_df[f'alpha{i_e}_before_err'].values]), 
                        linestyle='', marker='o', markersize=1, elinewidth=0.6, capsize=1)

    axs[0, 1].axhline(y=0, ls='--', color='k', linewidth=0.5)
    axs[0, 1].set_xlim(RANGE_SNR[0], RANGE_SNR[1])
    axs[0, 1].set_ylim(RANGE_alpha[0], RANGE_alpha[1])
    axs[0, 1].set_xscale('log')
    axs[0, 1].set_xlabel('median SNR')
    axs[0, 1].set_ylabel(f'alpha in e{i_e}')
    axs[0, 1].set_title(TITLE)

    # 3. as a function of R
    axs[0, 2].errorbar(alphaMap_df['R'].values, alphaMap_df[f'alpha{i_e}_before'].values, 
                        yerr = np.vstack([alphaMap_df[f'alpha{i_e}_before_err'].values, alphaMap_df[f'alpha{i_e}_before_err'].values]), 
                        linestyle='', marker='o', markersize=1, elinewidth=0.6, capsize=1)

    axs[0, 2].axhline(y=0, ls='--', color='k', linewidth=0.5)
    axs[0, 2].set_xlim(RANGE_R[0], RANGE_R[1])
    axs[0, 2].set_ylim(RANGE_alpha[0], RANGE_alpha[1])
    axs[0, 2].set_xlabel('median R')
    axs[0, 2].set_ylabel(f'alpha in e{i_e}')
    axs[0, 2].set_title(TITLE)

    ####### after correction
    TITLE = 'after alphaRecal'

    # 1. in 2D
    norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
    axs[1, 0].scatter(alphaMap_df['R'].values, alphaMap_df['SNR'].values, c=alphaMap_df[f'alpha{i_e}_after'].values,
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
    axs[1, 1].errorbar(alphaMap_df['SNR'].values, alphaMap_df[f'alpha{i_e}_after'].values, 
                        yerr = np.vstack([alphaMap_df[f'alpha{i_e}_after_err'].values, alphaMap_df[f'alpha{i_e}_after_err'].values]), 
                        linestyle='', marker='o', markersize=1, elinewidth=0.6, capsize=1)

    axs[1, 1].axhline(y=0, ls='--', color='k', linewidth=0.5)
    axs[1, 1].set_xlim(RANGE_SNR[0], RANGE_SNR[1])
    axs[1, 1].set_ylim(RANGE_alpha[0], RANGE_alpha[1])
    axs[1, 1].set_xscale('log')
    axs[1, 1].set_xlabel('median SNR')
    axs[1, 1].set_ylabel(f'alpha in e{i_e}')
    axs[1, 1].set_title(TITLE)

    # 3. as a function of R
    axs[1, 2].errorbar(alphaMap_df['R'].values, alphaMap_df[f'alpha{i_e}_after'].values, 
                        yerr = np.vstack([alphaMap_df[f'alpha{i_e}_after_err'].values, alphaMap_df[f'alpha{i_e}_after_err'].values]), 
                        linestyle='', marker='o', markersize=1, elinewidth=0.6, capsize=1)

    axs[1, 2].axhline(y=0, ls='--', color='k', linewidth=0.5)
    axs[1, 2].set_xlim(RANGE_R[0], RANGE_R[1])
    axs[1, 2].set_ylim(RANGE_alpha[0], RANGE_alpha[1])
    axs[1, 2].set_xlabel('median R')
    axs[1, 2].set_ylabel(f'alpha in e{i_e}')
    axs[1, 2].set_title(TITLE)

    ##### save plot
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()
    print(f"QC plot for alpha{i_e} in SNR-R plane saved as", outpath)


############## as a function of ZB
# >>>>>>>>>>>>>>>>>>>>>> load catalogue
cata = pd.read_feather(inpath)
## used columns
cata = cata[[col_ZB, 
            col_e1_before, col_e2_before, 
            col_e1_after, col_e2_after,
            col_psf_e1, col_psf_e2,
            col_weight_after
            ]]

# >>>>>>>>>>>>>>>>>>>>>> calculate alpha as a function of ZB

# in ZB bins
cata.loc[:, 'bin_ZB'] = pd.cut(cata[col_ZB], Z_B_edges, 
                                    right=True, labels=False)
# print('>>>>> bin_ZB', np.sort(cata['bin_ZB']))

# group based on binning
cata_grouped = cata.groupby(by=['bin_ZB'])
Ngroups = cata_grouped.ngroups
del cata

# loop over bins and calcuate alpha
## wanted info
alphaMap_df = pd.DataFrame(0, 
                index=np.arange(Ngroups), 
                columns=[ 
                'alpha1_before', 'alpha1_before_err', 
                'alpha1_after', 'alpha1_after_err',
                'alpha2_before', 'alpha2_before_err',
                'alpha2_after', 'alpha2_after_err'])

i_group = 0
for name, group in cata_grouped:

    # out PSF 
    e1_psf = np.array(group[col_psf_e1])
    e2_psf = np.array(group[col_psf_e2])

    # LF weights
    weight_out = np.array(group[col_weight_after])
    # out shear
    e1_out_list = [np.array(group[col_e1_before]), np.array(group[col_e1_after])] 
    e2_out_list = [np.array(group[col_e2_before]), np.array(group[col_e2_after])]
    col_save_list = ['before', 'after']
    del group

    # >>>>>>>>>>>>>>>> calculate alpha
    for i_e, e1_out in enumerate(e1_out_list):
        e2_out = e2_out_list[i_e]
        col_save = col_save_list[i_e]
    
        # calculate alpha using least square
        ## e1
        mod_wls = sm.WLS(e1_out, sm.add_constant(e1_psf), weights=weight_out)
        res_wls = mod_wls.fit()
        alphaMap_df.loc[i_group, f'alpha1_{col_save}'] = res_wls.params[1]
        alphaMap_df.loc[i_group, f'alpha1_{col_save}_err'] = (res_wls.cov_params()[1, 1])**0.5
        del res_wls, mod_wls, e1_out
        ## e2
        mod_wls = sm.WLS(e2_out, sm.add_constant(e2_psf), weights=weight_out)
        res_wls = mod_wls.fit()
        alphaMap_df.loc[i_group, f'alpha2_{col_save}'] = res_wls.params[1]
        alphaMap_df.loc[i_group, f'alpha2_{col_save}_err'] = (res_wls.cov_params()[1, 1])**0.5
        del res_wls, mod_wls, e2_out

    del weight_out, e1_psf, e2_psf, e1_out_list, e2_out_list

    i_group += 1

del cata_grouped

# >>>>>>>>>>>>>>>>>>>>>> plot

fig, axs = plt.subplots(1, 2)
fig.suptitle(os.path.basename(inpath))
outpath = outprefix + f'_ZBbins.png'

# 1. for e1
axs[0].errorbar(np.arange(Ngroups)+1, alphaMap_df[f'alpha1_before'].values,  
                    yerr = np.vstack([alphaMap_df[f'alpha1_before_err'].values, alphaMap_df[f'alpha1_before_err'].values]), 
                    color='k', label='before alphaRecal_step2',
                    linestyle='', marker='x', markersize=2, elinewidth=1, capsize=2)
axs[0].errorbar(np.arange(Ngroups)+1, alphaMap_df[f'alpha1_after'].values,  
                    yerr = np.vstack([alphaMap_df[f'alpha1_after_err'].values, alphaMap_df[f'alpha1_after_err'].values]), 
                    color='r', label='after alphaRecal_step2',
                    linestyle='', marker='o', markersize=2, elinewidth=1, capsize=2)

plt.legend(frameon=True, loc='best')
axs[0].axhline(y=0, ls='--', color='gray', linewidth=0.5)
axs[0].set_ylim(RANGE_alpha_ZB[0], RANGE_alpha_ZB[1])
axs[0].set_xlabel('Z_B bins')
axs[0].set_ylabel('alpha1')

# 2. for e2
axs[1].errorbar(np.arange(Ngroups)+1, alphaMap_df[f'alpha2_before'].values,  
                    yerr = np.vstack([alphaMap_df[f'alpha2_before_err'].values, alphaMap_df[f'alpha2_before_err'].values]), 
                    color='k', label='before alphaRecal_step2',
                    linestyle='', marker='x', markersize=2, elinewidth=1, capsize=2)
axs[1].errorbar(np.arange(Ngroups)+1, alphaMap_df[f'alpha2_after'].values,  
                    yerr = np.vstack([alphaMap_df[f'alpha2_after_err'].values, alphaMap_df[f'alpha2_after_err'].values]), 
                    color='r', label='after alphaRecal_step2',
                    linestyle='', marker='o', markersize=2, elinewidth=1, capsize=2)

plt.legend(frameon=True, loc='best')
axs[1].axhline(y=0, ls='--', color='gray', linewidth=0.5)
axs[1].set_ylim(RANGE_alpha_ZB[0], RANGE_alpha_ZB[1])
axs[1].set_xlabel('Z_B bins')
axs[1].set_ylabel('alpha2')

##### save plot
plt.tight_layout()
plt.savefig(outpath, dpi=300)
plt.close()
print(f"QC plot for alpha in ZB bins saved as", outpath)

