# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-08-29 13:45:51
# @Last Modified by:   lshuns
# @Last Modified time: 2022-08-30 15:36:48

### compare Z_B to the true redshift 
###### in magnitude bins

import numpy as np 
import pandas as pd 
from scipy import stats
import math

import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator, LogLocator, NullFormatter, NullLocator


# ++++++++++++++++++ general info

# the catalogue
inpath = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_photo_noSG_shear_m283m283_rot_0.feather'
cata = pd.read_feather(inpath)
## only galaxies
cata = cata[cata['id_input']>0]
cata.reset_index(drop=True, inplace=True)
## used columns
cata = cata[['MAG_AUTO', 'Z_B', 'redshift_input']]

# where to save
outpath = 'show'
outpath = './plots/BPZ_performance.pdf'

# the magnitude bin
mag_bins = [20, 22, 23, 24, 25]

# the redshift range
XRANGE = [0, 2.5] # the true redshift 
YRANGE = [0, 2.5] # the Z_B

# plot properties
N_rows = 1
N_cols = 4
FIGSIZE=[16, 5]
font_size = 20
usetex = True
nbins = 30
COLOR_MAP='Blues'
# XLABEL = r'$z_{\rm true}$'
# YLABEL = r'$z_{\rm B}$'
XLABEL = 'True redshift'
YLABEL = 'Photometric redshift'
xtick_spe = [[0, 1, 2], ['0', '1', '2']]
ytick_spe = [[0, 1, 2], ['0', '1', '2']]

# ++++++++++++++++++ plot 

# font size
plt.rc('font', size=font_size)
# tex
plt.rcParams["text.usetex"] = usetex

if outpath != 'show':
    backend_orig = plt.get_backend()
    plt.switch_backend("agg")

# how many sub-plots
N_plots = len(mag_bins) - 1
fig = plt.figure()
fig, axs = plt.subplots(N_rows, N_cols, sharex=True, sharey=True, figsize=FIGSIZE)
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)

# for x,y label
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
plt.xlabel(XLABEL)
plt.ylabel(YLABEL)

# some general relation
SRANGE = [[XRANGE[0], XRANGE[1]], [YRANGE[0], YRANGE[1]]]

# loop over sub-plots
i_plot = 0
for i_row in range(N_rows):
    for i_col in range(N_cols):
        if i_plot >= N_plots:
            if N_rows == 1:
                axs[i_col].axis('off')
            elif N_cols == 1:
                axs[i_row].axis('off')
            else:
                axs[i_row, i_col].axis('off')
        else:
            if (N_rows==1) and (N_cols == 1):
                ax = axs
            elif N_rows == 1:
                ax = axs[i_col]
            elif N_cols == 1:
                ax = axs[i_row]
            else:
                ax = axs[i_row, i_col]

            # get the values
            bin_min = mag_bins[i_plot]
            bin_max = mag_bins[i_plot + 1]
            i_plot += 1
            cata_selec = cata[(cata['MAG_AUTO']>bin_min) & (cata['MAG_AUTO']<=bin_max)].copy()

            # plot the distribution
            h = ax.hist2d(cata_selec['redshift_input'].values, cata_selec['Z_B'].values, 
                    bins=nbins, range=SRANGE, cmap=COLOR_MAP, 
                    density=True)

            # the statistics
            val = (cata_selec['Z_B'].values - cata_selec['redshift_input'].values) / (1+cata_selec['redshift_input'].values)
            del cata_selec
            sigma_m = stats.median_abs_deviation(val)
            eta_3 = np.sum(np.abs(val)>3*sigma_m) / len(val)
            zeta_15 = np.sum(np.abs(val)>0.15) / len(val)
            ## show
            LABEL = f'{bin_min}' + r' $<$ MAG\_AUTO $\leq$ ' + f'{bin_max}' + '\n'
            LABEL += '    ' + r'$\sigma_{\rm m}=$' + f'  {sigma_m:.3f}' + '\n'
            LABEL += '    ' + r'$\eta_3=$' + f'  {eta_3:.3f}' + '\n'
            LABEL += '    ' + r'$\zeta_{0.15}=$' + f' {zeta_15:.3f}' + '\n'
            ax.text(0.05, 0.6, LABEL, transform=ax.transAxes, fontsize=18)

            # plot the lines
            ## central 
            xval = np.arange(99)
            yval = np.arange(99)
            tmp = ax.plot(xval, yval, 
                            color='k', linestyle='--', linewidth=1.5)
            # ## outlier lines
            # yval = xval + 0.15 * (1+xval)
            # tmp = ax.plot(xval, yval, 
            #                 color='k', linestyle=':', linewidth=1.5)
            # yval = xval - 0.15 * (1+xval)
            # tmp = ax.plot(xval, yval, 
            #                 color='k', linestyle=':', linewidth=1.5)
            ## sigma lines
            yval = xval + sigma_m * (1+xval)
            tmp = ax.plot(xval, yval, 
                            color='k', linestyle=':', linewidth=1.5)
            yval = xval - sigma_m * (1+xval)
            tmp = ax.plot(xval, yval, 
                            color='k', linestyle=':', linewidth=1.5)

            # some general setting
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.set_xlim(XRANGE[0], XRANGE[1])
            ax.set_ylim(YRANGE[0], YRANGE[1])

            if xtick_spe is not None:
                ax.set_xticks(xtick_spe[0])
                ax.set_xticklabels(xtick_spe[1])
            if ytick_spe is not None:
                ax.set_yticks(ytick_spe[0])
                ax.set_yticklabels(ytick_spe[1])

fig.tight_layout()

if outpath == 'show':
    plt.show()
    plt.close()
else:
    plt.savefig(outpath, dpi=300)
    plt.close()
    plt.switch_backend(backend_orig)
    print("2D histogram plot saved as", outpath)    




