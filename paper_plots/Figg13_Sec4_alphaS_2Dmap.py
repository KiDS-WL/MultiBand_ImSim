# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-05-31 12:32:06
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-11 15:45:03

### calculate the alpha map 
###### in variance

import pandas as pd 
import numpy as np

import statsmodels.api as sm 

import plotting

# >>>>>>>>>>>>>>>>>>> I / O

# the catalogue for the data
outpath = './plots/alphaS_2Dmap_data.pdf'
TITLE = 'In measurement variance'
inpath = '/disks/shear16/ssli/KiDS/K1000_forSKiLLS/kids_photo_LF_321_shear_noSG_noWeiCut_newCut_817tiles.feather'
## column names
col_snr = 'model_SNratio'
col_var = '2D_measurement_variance'
col_wei = 'weight'
col_e1 = 'raw_e1'
col_e2 = 'raw_e2'
col_psf_e1 = 'PSF_e1'
col_psf_e2 = 'PSF_e2'
## load and select
cata = pd.read_feather(inpath)
cata = cata.loc[(cata[col_wei]>0), 
                ['R', col_snr, col_var, col_wei,
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

    # variance
    Var = np.array(group[col_var])

    # out shear
    e1_out = np.array(group[col_e1])
    e2_out = np.array(group[col_e2])
    emod_out = np.hypot(e1_out, e2_out)

    # out PSF 
    e1_psf = np.array(group[col_psf_e1])
    e2_psf = np.array(group[col_psf_e2])
    del group

    # rotate to PSF frame
    mask_tmp = emod_out > 0
    e_psf_rotated = e1_out * e1_psf + e2_out * e2_psf
    e_psf_rotated[mask_tmp] /= emod_out[mask_tmp]
    del mask_tmp, e1_psf, e2_psf, e1_out, e2_out, emod_out

    # calculate alpha using least square
    mod_ols = sm.OLS(Var, sm.add_constant(e_psf_rotated))
    res_ols = mod_ols.fit()
    alpha[i_group] = res_ols.params[1]

    i_group += 1

# >>>>>>>>>>>>>>>>>>>>>> plot

xval = SNR_median
yval = R_median
cval = alpha

XRANGE = [10**0.7, 10**2.5]
YRANGE = [0.15, 0.95]

cmap = 'coolwarm'
cmin = -0.2
cmax = 0.2

POINTS = 15

XLABEL = r'$\nu_{\rm SN}$'
YLABEL = r'$\mathcal{R}$'
bar_label = r'$\alpha_S$'

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
                font_size=14, usetex=True)

