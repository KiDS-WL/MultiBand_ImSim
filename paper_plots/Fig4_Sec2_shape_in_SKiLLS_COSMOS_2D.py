# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-09-24 17:05:52
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-24 17:07:55

### compare learned morphology between SKiLLS and COSMOS
###### two-dimensional kernal density plots

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

############################################ size vs. e (mag bins)

# out info
outpath = 'show'
outpath = './plots/size_e_magbins.pdf'

# binning column
## NOTE: 2x2
bin_label = r' $<m_{\rm r}<$ '
bin_col_s = 'r_SDSS_apparent_corr'
bin_col_c = 'MAG_AUTO_deepAW'
bin_mins = [21, 22, 23, 24]
bin_maxs = [21.5, 22.5, 23.5, 24.5]
N_plots = len(bin_mins)

# plot columns
x_col_s = 'Re_arcsec'
x_col_c = 'shape/Re'
XRANGE = (0, 2.0)
XTICKS = [0.5, 1.0, 1.5]
x_edges = np.linspace(XRANGE[0], XRANGE[1], 10)
XLABEL = r'Half-light radius [arcsec]'

y_col_s = 'e'
y_col_c = 'e'
YRANGE = (0, 1)
YTICKS = [0.2, 0.4, 0.6, 0.8]
YLABEL = 'Ellipticity'

# plot-related
plt.rc('font', size=16)
plt.rcParams['text.usetex'] = True

# start plotting
fig, axs = plt.subplots(2, 2, sharex=True , sharey=True)
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
sns.set_style("white")
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

        # build the plot dataframe
        cata_sns = []
        ## COSMOS
        cata_sns.append(pd.DataFrame({
            'x': cata_c.loc[mag_mask_c, x_col_c].values, 
            'y': cata_c.loc[mag_mask_c, y_col_c].values, 
            'kind': 'COSMOS'}))
        ## SKiLLS
        # cata_sns.append(pd.DataFrame({
        #     'x': cata_s.loc[mag_mask_s, x_col_s].values, 
        #     'y': cata_s.loc[mag_mask_s, y_col_s].values, 
        #     'kind': 'SKiLLS'}).sample(n=N_c, replace=False))
        cata_sns.append(pd.DataFrame({
            'x': cata_s.loc[mag_mask_s, x_col_s].values, 
            'y': cata_s.loc[mag_mask_s, y_col_s].values, 
            'kind': 'SKiLLS'}))
        ## combine
        cata_sns = pd.concat(cata_sns)
        cata_sns.reset_index(drop=True, inplace=True)
        print('total number in sns plot', len(cata_sns))

        # draw contours
        sns.kdeplot(x=cata_sns.loc[cata_sns['kind']=='COSMOS', 'x'], 
            y=cata_sns.loc[cata_sns['kind']=='COSMOS', 'y'], 
            colors='r', levels=[0.2, 0.4, 0.6, 0.8], ax=ax, label='COSMOS', linestyles='-')
        sns.kdeplot(x=cata_sns.loc[cata_sns['kind']=='SKiLLS', 'x'], 
            y=cata_sns.loc[cata_sns['kind']=='SKiLLS', 'y'], 
            colors='b', levels=[0.2, 0.4, 0.6, 0.8], ax=ax, label='SKiLLS', linestyles='--')

        ax.set(xlabel=None)
        ax.set(ylabel=None)

        ax.set_xticks(XTICKS)
        ax.set_yticks(YTICKS)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

        ax.text(0.35, 0.8, LABEL, transform=ax.transAxes)

line_c, = ax.plot([], color='r', linestyle='-')
line_s, = ax.plot([], color='b', linestyle='--')
fig.legend([line_c, line_s], ['COSMOS', 'SKiLLS'], loc='upper center', bbox_to_anchor=(0.5, 1.0),
             fancybox=True, shadow=True, ncol=2)

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

############################################ size vs. sersic index (mag bins)

# out info
outpath = 'show'
outpath = './plots/size_n_magbins.pdf'

# binning column
## NOTE: 2x2
bin_label = r' $<m_{\rm r}<$ '
bin_col_s = 'r_SDSS_apparent_corr'
bin_col_c = 'MAG_AUTO_deepAW'
bin_mins = [21, 22, 23, 24]
bin_maxs = [21.5, 22.5, 23.5, 24.5]
N_plots = len(bin_mins)

# plot columns
x_col_s = 'Re_arcsec'
x_col_c = 'shape/Re'
XRANGE = (0, 2.0)
XTICKS = [0.5, 1.0, 1.5]
x_edges = np.linspace(XRANGE[0], XRANGE[1], 10)
XLABEL = r'Half-light radius [arcsec]'

y_col_s = 'shape/sersic_n'
y_col_c = 'shape/sersic_n'
YRANGE = (0.3, 30.)
YTICKS = [1., 10.]
YLABEL = r'S\'ersic index'

# plot-related
plt.rc('font', size=16)
plt.rcParams['text.usetex'] = True

# start plotting
fig, axs = plt.subplots(2, 2, sharex=True , sharey=True)
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
sns.set_style("white")
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

        # build the plot dataframe
        cata_sns = []
        ## COSMOS
        cata_sns.append(pd.DataFrame({
            'x': cata_c.loc[mag_mask_c, x_col_c].values, 
            'y': cata_c.loc[mag_mask_c, y_col_c].values, 
            'kind': 'COSMOS'}))
        ## SKiLLS
        # cata_sns.append(pd.DataFrame({
        #     'x': cata_s.loc[mag_mask_s, x_col_s].values, 
        #     'y': cata_s.loc[mag_mask_s, y_col_s].values, 
        #     'kind': 'SKiLLS'}).sample(n=N_c, replace=False))
        cata_sns.append(pd.DataFrame({
            'x': cata_s.loc[mag_mask_s, x_col_s].values, 
            'y': cata_s.loc[mag_mask_s, y_col_s].values, 
            'kind': 'SKiLLS'}))
        ## combine
        cata_sns = pd.concat(cata_sns)
        cata_sns.reset_index(drop=True, inplace=True)
        print('total number in sns plot', len(cata_sns))

        # draw contours
        sns.kdeplot(x=cata_sns.loc[cata_sns['kind']=='COSMOS', 'x'], 
            y=cata_sns.loc[cata_sns['kind']=='COSMOS', 'y'], 
            colors='r', linestyles='-',
            levels=[0.2, 0.4, 0.6, 0.8], ax=ax, label='COSMOS', log_scale=(False, True))
        sns.kdeplot(x=cata_sns.loc[cata_sns['kind']=='SKiLLS', 'x'], 
            y=cata_sns.loc[cata_sns['kind']=='SKiLLS', 'y'], 
            colors='b', linestyles='--',
            levels=[0.2, 0.4, 0.6, 0.8], ax=ax, label='SKiLLS', log_scale=(False, True))

        ax.set(xlabel=None)
        ax.set(ylabel=None)

        ax.set_xticks(XTICKS)
        ax.set_yticks(YTICKS)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        # ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(LogLocator(base=10., subs=None, numticks=10))

        ax.text(0.35, 0.8, LABEL, transform=ax.transAxes)

line_c, = ax.plot([], color='r', linestyle='-')
line_s, = ax.plot([], color='b', linestyle='--')
fig.legend([line_c, line_s], ['COSMOS', 'SKiLLS'], loc='upper center', bbox_to_anchor=(0.5, 1.0),
             fancybox=True, shadow=True, ncol=2)

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