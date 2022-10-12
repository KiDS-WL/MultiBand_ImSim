# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-09-29 14:04:31
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-29 16:55:05

### m difference by changing galaxy morphology

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

######### size
# get y val
inpath_list = ['../sensitivity_test/galaxy/outputs/dm_sizeU5_magauto.csv',
                '../sensitivity_test/galaxy/outputs/dm_sizeU10_magauto.csv',
                '../sensitivity_test/galaxy/outputs/dm_sizeU15_magauto.csv']
yval_s = np.zeros(4)
yerr_s = np.zeros(4)
for i_file, inpath in enumerate(inpath_list):
    cata = pd.read_csv(inpath)
    # dm as y
    yval_s[i_file+1] = (cata.loc[0, 'm1'] + cata.loc[0, 'm2'])/2.
    yerr_s[i_file+1] = (cata.loc[0, 'm1_err'] + cata.loc[0, 'm2_err'])/2.
    del cata
# get number counts
inpath_list = ['/disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0.feather',
                '/disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_sizeU5.feather',
                '/disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_sizeU10.feather',
                '/disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_sizeU15.feather']
numbers_s = []
for inpath in inpath_list:
    cata = pd.read_feather(inpath)
    numbers_s.append(cata['Re_arcsec'].values)
    del cata

######### n
# get y val
inpath_list = ['../sensitivity_test/galaxy/outputs/dm_nU10_magauto.csv',
                '../sensitivity_test/galaxy/outputs/dm_nU20_magauto.csv',
                '../sensitivity_test/galaxy/outputs/dm_nU30_magauto.csv']
yval_n = np.zeros(4)
yerr_n = np.zeros(4)
for i_file, inpath in enumerate(inpath_list):
    cata = pd.read_csv(inpath)
    # dm as y
    yval_n[i_file+1] = (cata.loc[0, 'm1'] + cata.loc[0, 'm2'])/2.
    yerr_n[i_file+1] = (cata.loc[0, 'm1_err'] + cata.loc[0, 'm2_err'])/2.
    del cata
# get number counts
inpath_list = ['/disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0.feather',
                '/disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_nU10.feather',
                '/disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_nU20.feather',
                '/disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_nU30.feather']
numbers_n = []
for inpath in inpath_list:
    cata = pd.read_feather(inpath)
    numbers_n.append(cata['shape/sersic_n'].values)
    del cata

######### q
# get y val
inpath_list = ['../sensitivity_test/galaxy/outputs/dm_qU5_magauto.csv',
                '../sensitivity_test/galaxy/outputs/dm_qU10_magauto.csv',
                '../sensitivity_test/galaxy/outputs/dm_qU15_magauto.csv']
yval_q = np.zeros(4)
yerr_q = np.zeros(4)
for i_file, inpath in enumerate(inpath_list):
    cata = pd.read_csv(inpath)
    # dm as y
    yval_q[i_file+1] = (cata.loc[0, 'm1'] + cata.loc[0, 'm2'])/2.
    yerr_q[i_file+1] = (cata.loc[0, 'm1_err'] + cata.loc[0, 'm2_err'])/2.
    del cata
# get number counts
inpath_list = ['/disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0.feather',
                '/disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_qU5.feather',
                '/disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_qU10.feather',
                '/disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_qU15.feather']
numbers_q = []
for inpath in inpath_list:
    cata = pd.read_feather(inpath)
    numbers_q.append(cata['BA'].values)
    del cata

######### collect values
xvals = [[0, 1, 2, 3]] * 3
xtick_spe = [xvals[0], ['0', r'$1\sigma$', r'$2\sigma$', r'$3\sigma$']]

yvals = [yval_s, yval_q, yval_n]
yerrs = [yerr_s, yerr_q, yerr_n]
numbers = [numbers_s, numbers_q, numbers_n]
del numbers_s, numbers_q, numbers_n

# 0, 1, 2, 3
COLORs = ['k', 'darkred', 'darkorange', 'darkgreen']
LINEs = ['solid', 'solid', 'solid', 'solid']
LWs = [1, 1, 1, 1]
LABELs = ['fiducial', r'$1\sigma$ shift', r'$2\sigma$ shift', r'$3\sigma$ shift']
POINTs = ['o', 'v', '^', '>']

# where to save
outpath = 'show'
outpath = './plots/dm_galaxy.pdf'

# the error budget
m_error_budget = 1

XRANGE_p = [-1, 4]
YRANGE_p = [-1.6, 1.6]

nbins = 40
loc_legend=None
XRANGE_h_list = [[0.05, 2], [0.05, 0.9], [0.6, 5.9]]
xlog_list = [True, False, False]

# plot properties
FIGSIZE = [9, 4]
font_size = 12
usetex = True
MS = 4

# XLABEL_p = 'Change fraction'
YLABEL_p = r'$\Delta m$ [$\times 10^{-2}$]'

YLABEL_h = 'Probability density'
XLABEL_h_list = ['half-light radius [arcsec]', 'axis ratio', r'S\'ersic index']

# ++++++++++++++++++ plot 

# font size
plt.rc('font', size=font_size)
# tex
plt.rcParams["text.usetex"] = usetex

if outpath != 'show':
    backend_orig = plt.get_backend()
    plt.switch_backend("agg")

# how many sub-plots
fig, axs = plt.subplots(2, 3, sharex=False, sharey=False, figsize=FIGSIZE)
# fig.subplots_adjust(hspace=0)
# fig.subplots_adjust(wspace=0)

# loop over sub-plots
## start with m
i_row = 0
for i_col in range(3):
    ax = axs[i_row, i_col]

    # get the values
    xval = xvals[i_col]
    yval = yvals[i_col] * 1e2
    yerr = yerrs[i_col] * 1e2

    # error budget
    ax.fill_between(np.array(XRANGE_p), 
                -m_error_budget, m_error_budget,
               edgecolor='gray', facecolor='white', alpha=0.3, hatch='/')

    # plot 
    for i_point, x_val in enumerate(xval):
        ax.errorbar(x_val, yval[i_point], yerr=np.vstack([yerr[i_point], yerr[i_point]]),
                    color=COLORs[i_point], linestyle='', marker=POINTs[i_point], markersize=MS)

    # the labels
    # if i_col == 1:
    #     ax.set_xlabel(XLABEL_p)
    if i_col == 0:
        ax.set_ylabel(YLABEL_p)

    # some general setting
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(XRANGE_p[0], XRANGE_p[1])
    ax.set_ylim(YRANGE_p[0], YRANGE_p[1])

    if 'xtick_spe' in locals():
        ax.set_xticks(xtick_spe[0])
        ax.set_xticklabels(xtick_spe[1])
    if 'ytick_spe' is locals():
        ax.set_yticks(ytick_spe[0])
        ax.set_yticklabels(ytick_spe[1])

## then number counts
i_row = 1
for i_col in range(3):
    ax = axs[i_row, i_col]

    number = numbers[i_col]
    xlog = xlog_list[i_col]
    XRANGE = XRANGE_h_list[i_col]

    for i_para, para in enumerate(number):

        if xlog:
            logbins = np.logspace(np.log10(XRANGE[0]), np.log10(XRANGE[1]), nbins)
            if i_para == 0 :
                ax.hist(x=para, 
                    bins=logbins, 
                    range=XRANGE,
                    density=True, 
                    color=COLORs[i_para], 
                    label=LABELs[i_para], 
                    alpha=0.3, 
                    linestyle=LINEs[i_para],
                    linewidth=LWs[i_para])
            else:
                ax.hist(x=para, 
                    bins=logbins, 
                    range=XRANGE,
                    density=True, 
                    color=COLORs[i_para], 
                    label=LABELs[i_para], 
                    histtype='step', 
                    linestyle=LINEs[i_para],
                    linewidth=LWs[i_para])
        else:
            if i_para == 0:
                ax.hist(x=para, 
                    bins=nbins, 
                    range=XRANGE,
                    density=True, 
                    color=COLORs[i_para], 
                    label=LABELs[i_para], 
                    alpha=0.3, 
                    linestyle=LINEs[i_para],
                    linewidth=LWs[i_para])
            else:
                ax.hist(x=para, 
                    bins=nbins, 
                    range=XRANGE,
                    density=True, 
                    color=COLORs[i_para], 
                    label=LABELs[i_para], 
                    histtype='step', 
                    linestyle=LINEs[i_para],
                    linewidth=LWs[i_para])

    if xlog:
        ax.set_xscale('log')

    # the labels
    ax.set_xlabel(XLABEL_h_list[i_col])
    if i_col == 0:
        ax.set_ylabel(YLABEL_h)

    ax.set_xlim(XRANGE[0], XRANGE[1])

    if i_col==2:
        ax.legend(frameon=True, loc=loc_legend)

# fig.tight_layout()

if outpath == 'show':
    plt.show()
    plt.close()
else:
    plt.savefig(outpath, dpi=300)
    plt.close()
    plt.switch_backend(backend_orig)
    print("plot saved as", outpath)    




