# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-07-09 11:42:38
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-11 15:44:29

### calculate the alpha map 
###### in e12, report the mean

import pandas as pd 
import numpy as np

import statsmodels.api as sm 

import plotting

# >>>>>>>>>>>>>>>>>>> I / O

# the catalogue for the data
outpath = './plots/alpha_2Dmap_data.pdf'
TITLE = 'In measured ellipticity'
inpath = '/disks/shear16/ssli/KiDS/K1000_forSKiLLS/kids_photo_LF_321_shear_noSG_noWeiCut_newCut_817tiles_A1.feather'
## column names
col_snr = 'model_SNratio'
col_wei = 'AlphaRecalC_weight'
col_e1 = 'autocal_e1'
col_e2 = 'autocal_e2'
col_psf_e1 = 'PSF_e1'
col_psf_e2 = 'PSF_e2'
## load and select
cata = pd.read_feather(inpath)
cata = cata.loc[(cata[col_wei]>0), 
                ['R', col_snr, col_wei,
                col_e1, col_e2, col_psf_e1, col_psf_e2]]
cata.reset_index(drop=True, inplace=True)

# >>>>>>>>>>>>>>>>>>> general setups
## number of bins
N_R = 20
N_SNR = 20

# >>>>>>>>>>>>>>>>>>>> binning
# start with R bin
cata.loc[:, 'bin_R'] = pd.qcut(cata['R'], N_R, 
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
del cata

# loop over bins and calcuate alpha
## wanted info
R_median = np.zeros(N_R*N_SNR)
SNR_median = np.zeros(N_R*N_SNR)
alpha = np.zeros(N_R*N_SNR)

i_group = 0
for name, group in cata_grouped:

    # save median center
    R_median[i_group] = np.median(np.array(group['R']))
    SNR_median[i_group] = np.median(np.array(group[col_snr]))

    # >>>>>>>>>>>>>>>> calculate alpha
    # out shear
    e1_out = np.array(group[col_e1])
    e2_out = np.array(group[col_e2])
    weight_out = np.array(group[col_wei])
    # out PSF 
    e1_psf = np.array(group[col_psf_e1])
    e2_psf = np.array(group[col_psf_e2])
    del group
    # calculate alpha using least square
    ## e1
    mod_wls = sm.WLS(e1_out, sm.add_constant(e1_psf), weights=weight_out)
    res_wls = mod_wls.fit()
    alpha1 = res_wls.params[1]
    del res_wls, mod_wls
    ## e2
    mod_wls = sm.WLS(e2_out, sm.add_constant(e2_psf), weights=weight_out)
    res_wls = mod_wls.fit()
    alpha2 = res_wls.params[1]
    del weight_out, res_wls, mod_wls

    alpha[i_group] = (alpha1 + alpha2) / 2.

    i_group += 1

# >>>>>>>>>>>>>>>>>>>>>> plot alpha1

xval = SNR_median
yval = R_median
cval = alpha

XRANGE = [10**0.7, 10**2.5]
YRANGE = [0.15, 0.95]

cmap = 'coolwarm'
cmin = -0.5
cmax = 0.5

POINTS = 15

XLABEL = r'$\nu_{\rm SN}$'
YLABEL = r'$\mathcal{R}$'
bar_label = r'$|\boldsymbol{\alpha}|$'

texPacks = [r"\usepackage{amsmath}"]

xlog = True
ylog = False
invertY = True
ytick_min_label = False

plotting.ScatterPlotFunc(outpath,
                xval, yval, POINT=None, POINTS=POINTS, alpha=None,
                cval=cval, cmap=cmap, cmin=cmin, cmax=cmax,
                bar_loc=None, bar_ori=None, bar_label=bar_label, bar_tick=None,
                XRANGE=XRANGE, YRANGE=YRANGE,
                XLABEL=XLABEL, YLABEL=YLABEL, TITLE=TITLE,
                xtick_min_label=True, xtick_spe=None, ytick_min_label=ytick_min_label, ytick_spe=None,
                vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                xlog=xlog, invertX=False, ylog=ylog, invertY=invertY, loc_legend='best',
                font_size=14, usetex=True,
                texPacks=texPacks)