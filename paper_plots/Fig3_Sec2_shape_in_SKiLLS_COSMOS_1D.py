# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-12-17 10:56:08
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-12 14:57:10

### compare learned morphology between SKiLLS and COSMOS
###### magnitude-morphology relations in several redshift bins

import os
import sys

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy import stats
from matplotlib.ticker import AutoMinorLocator, LogLocator

############################### I/O

# SKiLLS
infile = '/disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0.feather'
cata_s = pd.read_feather(infile)
## ell
cata_s.loc[:, 'e'] = (1-cata_s['BA'].values)/(1+cata_s['BA'].values)

# COSMOS
infile = '/disks/shear10/ssli/ImSim/input/COSMOS_cata/cosmos_shape_z_ugriZYJHKs_selected_deepAW_only.feather'
cata_c = pd.read_feather(infile)
## ell
cata_c.loc[:, 'e'] = (1-cata_c['BA_GALFIT_HI'].values)/(1+cata_c['BA_GALFIT_HI'].values)

############################### magnitude vs. e (redshift bins)

# out info
outpath = 'show'
outpath = './plots/e_mag_zbins.pdf'

# binning column
## NOTE: 2x2
bin_label = r' $<z<$ '
bin_col_s = 'zobs'
bin_col_c = 'Z_optimal'
bin_mins = [0.1, 0.5, 0.9, 1.4]
bin_maxs = [0.3, 0.7, 1.2, 2.0]
N_plots = len(bin_mins)

# plot columns
x_col_s = 'r_SDSS_apparent_corr'
x_col_c = 'MAG_AUTO_deepAW'
XRANGE = (20, 25)
XTICKS = [21, 22, 23, 24] 
x_edges = np.linspace(XRANGE[0], XRANGE[1], 10)
XLABEL = r'Input magnitude ($r$-band)'

y_col_s = 'e'
y_col_c = 'e'
YRANGE = (0.1, 0.5)
YTICKS = [0.2, 0.3, 0.4]
YLABEL = 'Mean ellipticity'

# plot-related
plt.rc('font', size=16)
plt.rcParams['text.usetex'] = True

# error type
# ## bootstrap
nboot = 500
## std
# numpy_std = True

# start plot
fig, axs = plt.subplots(2, 2, sharex=True , sharey=True)
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
i_bin = 0
for i_row in range(2):
    for i_col in range(2):
        ax = axs[i_row, i_col]

        # data
        bin_min = bin_mins[i_bin]
        bin_max = bin_maxs[i_bin]
        i_bin += 1
        print('for bin', bin_min, bin_max)
        LABEL = str(bin_min) + bin_label + str(bin_max)

        # magnitude selection
        mag_mask_c = (cata_c[bin_col_c] > bin_min) & (cata_c[bin_col_c] < bin_max)
        N_c = np.sum(mag_mask_c)
        print('number in COSMOS', N_c)
        mag_mask_s = (cata_s[bin_col_s] > bin_min) & (cata_s[bin_col_s] < bin_max)
        N_s = np.sum(mag_mask_s)
        print('number in SKiLLS', N_s)

        # used parameters
        x_c_para = cata_c.loc[mag_mask_c, x_col_c].values
        y_c_para = cata_c.loc[mag_mask_c, y_col_c].values
        x_s_para = cata_s.loc[mag_mask_s, x_col_s].values
        y_s_para = cata_s.loc[mag_mask_s, y_col_s].values

        # get the average and errors
        xval = 0.5*(x_edges[:-1]+x_edges[1:])
        yval_s = np.zeros_like(xval)
        yerr_s = [np.zeros_like(xval), np.zeros_like(xval)]
        yval_c = np.zeros_like(xval)
        yerr_c = [np.zeros_like(xval), np.zeros_like(xval)]
        for i_xv in range(len(x_edges)-1):

            xmin = x_edges[i_xv]
            xmax = x_edges[i_xv+1]

            # SKiLLS
            mask_tmp = (x_s_para >= xmin) & (x_s_para < xmax)
            y_para_tmp = y_s_para[mask_tmp]
            ## bins with very few galaxies meaningless for mean or error
            if len(y_para_tmp) < 10:
                yval_s[i_xv] = None
                yerr_s[0][i_xv] = None
                yerr_s[1][i_xv] = None
            else:
                y_mean = np.mean(y_para_tmp)
                yval_s[i_xv] = y_mean
                ## errors
                if 'nboot' in locals():
                    ## get error from boots
                    Relications_index = np.array([np.random.randint(len(y_para_tmp), size=len(y_para_tmp)) for _ in range(nboot)])
                    mean_BS = np.mean(y_para_tmp[Relications_index], axis = 1)
                    sigma_BS = mean_BS - y_mean
                    err_low = np.abs(np.percentile(sigma_BS, 16))
                    err_high = np.abs(np.percentile(sigma_BS, 84))
                    ## assign
                    yerr_s[0][i_xv] = err_low
                    yerr_s[1][i_xv] = err_high
                else:
                    std_tmp = np.std(y_para_tmp)
                    yerr_s[0][i_xv] = std_tmp
                    yerr_s[1][i_xv] = std_tmp

            # COSMOS
            mask_tmp = (x_c_para >= xmin) & (x_c_para < xmax)
            y_para_tmp = y_c_para[mask_tmp]
            ## bins with very few galaxies meaningless for mean or error
            if len(y_para_tmp) < 10:
                yval_c[i_xv] = None
                yerr_c[0][i_xv] = None
                yerr_c[1][i_xv] = None
            else:
                y_mean = np.mean(y_para_tmp)
                yval_c[i_xv] = y_mean
                ## errors
                if 'nboot' in locals():
                    ## get error from boots
                    Relications_index = np.array([np.random.randint(len(y_para_tmp), size=len(y_para_tmp)) for _ in range(nboot)])
                    mean_BS = np.mean(y_para_tmp[Relications_index], axis = 1)
                    sigma_BS = mean_BS - y_mean
                    err_low = np.abs(np.percentile(sigma_BS, 16))
                    err_high = np.abs(np.percentile(sigma_BS, 84))
                    ## assign
                    yerr_c[0][i_xv] = err_low
                    yerr_c[1][i_xv] = err_high
                else:
                    std_tmp = np.std(y_para_tmp)
                    yerr_c[0][i_xv] = std_tmp
                    yerr_c[1][i_xv] = std_tmp

        # draw lines    
        line_c = ax.errorbar(xval, yval_c, yerr=yerr_c, color='r', label='COSMOS', linestyle='-', linewidth=1.5, marker='o', markersize=3, capsize=3)
        line_s = ax.errorbar(xval, yval_s, yerr=yerr_s, color='b', label='SKiLLS', linestyle='--', linewidth=1.5, marker='o', markersize=3, capsize=3)
        print('>>>>>>>>>>>>xval', xval)
        print('>>>>>>>>>>>>yval_c', yval_c)
        print('>>>>>>>>>>>>yval_s', yval_s)

        # labels set later
        ax.set(xlabel=None)
        ax.set(ylabel=None)

        # ticks are specified
        ax.set_xticks(XTICKS)
        ax.set_yticks(YTICKS)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

        # binning label
        ax.text(0.35, 0.8, LABEL, transform=ax.transAxes)

# legend on top
fig.legend([line_c, line_s], ['COSMOS', 'SKiLLS'], loc='upper center', bbox_to_anchor=(0.5, 1.0),
             fancybox=True, shadow=True, ncol=2)

# label for the whole image
fig.text(0.5, 0.02, XLABEL, ha='center')
fig.text(0.02, 0.5, YLABEL, va='center', rotation='vertical')

# plot range
plt.xlim((XRANGE[0], XRANGE[1]))
plt.ylim((YRANGE[0], YRANGE[1]))

# save or show
if outpath == 'show':
    plt.show()
else:
    plt.savefig(outpath, dpi=300)
    print('plot saved as', outpath)
plt.close()

############################### magnitude vs. size (redshift bins)

# out info
outpath = 'show'
outpath = './plots/size_mag_zbins.pdf'

# binning column
## NOTE: 2x2
bin_label = r' $<z<$ '
bin_col_s = 'zobs'
bin_col_c = 'Z_optimal'
bin_mins = [0.1, 0.5, 0.9, 1.4]
bin_maxs = [0.3, 0.7, 1.2, 2.0]
N_plots = len(bin_mins)

# plot columns
x_col_s = 'r_SDSS_apparent_corr'
x_col_c = 'MAG_AUTO_deepAW'
XRANGE = (20, 25)
XTICKS = [21, 22, 23, 24] 
x_edges = np.linspace(XRANGE[0], XRANGE[1], 10)
XLABEL = r'Input magnitude ($r$-band)'

y_col_s = 'Re_arcsec'
y_col_c = 'shape/Re'
YRANGE = (0, 2.0)
YTICKS = [0.5, 1.0, 1.5]
YLABEL = 'Mean half-light radius [arcsec]'

# plot-related
plt.rc('font', size=16)
plt.rcParams['text.usetex'] = True

# error type
# ## bootstrap
nboot = 500
## std
# numpy_std = True

# start plot
fig, axs = plt.subplots(2, 2, sharex=True , sharey=True)
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
i_bin = 0
for i_row in range(2):
    for i_col in range(2):
        ax = axs[i_row, i_col]

        # data
        bin_min = bin_mins[i_bin]
        bin_max = bin_maxs[i_bin]
        i_bin += 1
        print('for bin', bin_min, bin_max)
        LABEL = str(bin_min) + bin_label + str(bin_max)

        # magnitude selection
        mag_mask_c = (cata_c[bin_col_c] > bin_min) & (cata_c[bin_col_c] < bin_max)
        N_c = np.sum(mag_mask_c)
        print('number in COSMOS', N_c)
        mag_mask_s = (cata_s[bin_col_s] > bin_min) & (cata_s[bin_col_s] < bin_max)
        N_s = np.sum(mag_mask_s)
        print('number in SKiLLS', N_s)

        # used parameters
        x_c_para = cata_c.loc[mag_mask_c, x_col_c].values
        y_c_para = cata_c.loc[mag_mask_c, y_col_c].values
        x_s_para = cata_s.loc[mag_mask_s, x_col_s].values
        y_s_para = cata_s.loc[mag_mask_s, y_col_s].values

        # get the average and errors
        xval = 0.5*(x_edges[:-1]+x_edges[1:])
        yval_s = np.zeros_like(xval)
        yerr_s = [np.zeros_like(xval), np.zeros_like(xval)]
        yval_c = np.zeros_like(xval)
        yerr_c = [np.zeros_like(xval), np.zeros_like(xval)]
        for i_xv in range(len(x_edges)-1):

            xmin = x_edges[i_xv]
            xmax = x_edges[i_xv+1]

            # SKiLLS
            mask_tmp = (x_s_para >= xmin) & (x_s_para < xmax)
            y_para_tmp = y_s_para[mask_tmp]
            ## bins with very few galaxies meaningless for mean or error
            if len(y_para_tmp) < 10:
                yval_s[i_xv] = None
                yerr_s[0][i_xv] = None
                yerr_s[1][i_xv] = None
            else:
                y_mean = np.mean(y_para_tmp)
                yval_s[i_xv] = y_mean
                ## errors
                if 'nboot' in locals():
                    ## get error from boots
                    Relications_index = np.array([np.random.randint(len(y_para_tmp), size=len(y_para_tmp)) for _ in range(nboot)])
                    mean_BS = np.mean(y_para_tmp[Relications_index], axis = 1)
                    sigma_BS = mean_BS - y_mean
                    err_low = np.abs(np.percentile(sigma_BS, 16))
                    err_high = np.abs(np.percentile(sigma_BS, 84))
                    ## assign
                    yerr_s[0][i_xv] = err_low
                    yerr_s[1][i_xv] = err_high
                else:
                    std_tmp = np.std(y_para_tmp)
                    yerr_s[0][i_xv] = std_tmp
                    yerr_s[1][i_xv] = std_tmp

            # COSMOS
            mask_tmp = (x_c_para >= xmin) & (x_c_para < xmax)
            y_para_tmp = y_c_para[mask_tmp]
            ## bins with very few galaxies meaningless for mean or error
            if len(y_para_tmp) < 10:
                yval_c[i_xv] = None
                yerr_c[0][i_xv] = None
                yerr_c[1][i_xv] = None
            else:
                y_mean = np.mean(y_para_tmp)
                yval_c[i_xv] = y_mean
                ## errors
                if 'nboot' in locals():
                    ## get error from boots
                    Relications_index = np.array([np.random.randint(len(y_para_tmp), size=len(y_para_tmp)) for _ in range(nboot)])
                    mean_BS = np.mean(y_para_tmp[Relications_index], axis = 1)
                    sigma_BS = mean_BS - y_mean
                    err_low = np.abs(np.percentile(sigma_BS, 16))
                    err_high = np.abs(np.percentile(sigma_BS, 84))
                    ## assign
                    yerr_c[0][i_xv] = err_low
                    yerr_c[1][i_xv] = err_high
                else:
                    std_tmp = np.std(y_para_tmp)
                    yerr_c[0][i_xv] = std_tmp
                    yerr_c[1][i_xv] = std_tmp

        # number counts
        ax.hist(x=x_c_para, bins=30, range=XRANGE,
                    density=True, 
                    color='r', alpha=0.3)
        ax.hist(x=x_s_para, bins=30,
                    density=True, range=XRANGE,
                    color='b', alpha=0.3)

        # draw lines    
        line_c = ax.errorbar(xval, yval_c, yerr=yerr_c, color='r', label='COSMOS', linestyle='-', linewidth=1.5, marker='o', markersize=3, capsize=3)
        line_s = ax.errorbar(xval, yval_s, yerr=yerr_s, color='b', label='SKiLLS', linestyle='--', linewidth=1.5, marker='o', markersize=3, capsize=3)
        print('>>>>>>>>>>>>xval', xval)
        print('>>>>>>>>>>>>yval_c', yval_c)
        print('>>>>>>>>>>>>yval_s', yval_s)

        # labels set later
        ax.set(xlabel=None)
        ax.set(ylabel=None)

        # ticks are specified
        ax.set_xticks(XTICKS)
        ax.set_yticks(YTICKS)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

        # binning label
        ax.text(0.35, 0.8, LABEL, transform=ax.transAxes)

# legend on top
fig.legend([line_c, line_s], ['COSMOS', 'SKiLLS'], loc='upper center', bbox_to_anchor=(0.5, 1.0),
             fancybox=True, shadow=True, ncol=2)

# label for the whole image
fig.text(0.5, 0.02, XLABEL, ha='center')
fig.text(0.02, 0.5, YLABEL, va='center', rotation='vertical')

# plot range
plt.xlim((XRANGE[0], XRANGE[1]))
plt.ylim((YRANGE[0], YRANGE[1]))

# save or show
if outpath == 'show':
    plt.show()
else:
    plt.savefig(outpath, dpi=300)
    print('plot saved as', outpath)

plt.close()
