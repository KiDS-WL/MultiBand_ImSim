# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-08-03 14:14:05
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-12 14:54:48

### compare input photometry between SKiLLS and COSMOS2015

import math
import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator, LogLocator

# ++++++++++++++ I/O
# SKiLLS cata
area_sim = 108.
infile_sim = '/disks/shear10/ssli/ImSim/input/SURFS_cata/SURFS_SHARK_LC_100sqdeg_testcols_calib.feather'
cols_sim = ['u_SDSS_apparent',
            'r_SDSS_apparent',
            'i_SDSS_apparent',
            'z_SDSS_apparent',
            'Y_VISTA_apparent',
            'J_VISTA_apparent',
            'H_VISTA_apparent',
            'K_VISTA_apparent']

# cosmos cata
area_data = 1.38
infile_data = '/disks/shear10/ssli/ImSim/input/COSMOS_cata/outside/cosmos_2015_dr_2_1_138deg2_galONLY.fits'
cols_data = ['u_mag_auto',
                'r_mag_auto',
                'ip_mag_auto',
                'zpp_mag_auto',
                'y_mag_auto',
                'j_mag_auto',
                'h_mag_auto',
                'ks_mag_auto']

# ++++++++++++++ load catalogue
## commom names
cols_goodnames = ['mag_u_app',
                'mag_r_app',
                'mag_i_app',
                'mag_Z_app',
                'mag_Y_app',
                'mag_J_app',
                'mag_H_app',
                'mag_Ks_app']

## simulation
sim_cata = pd.read_feather(infile_sim)
### rename
cols_rename = {}
for i_col, col_name in enumerate(cols_sim):
     cols_rename[col_name] = cols_goodnames[i_col]
sim_cata.rename(columns=cols_rename, inplace=True)

## data 
with fits.open(infile_data) as hdul:
    data_cata = Table(hdul[1].data).to_pandas()
### desired columns
data_cata = data_cata[cols_data]
### rename
cols_rename = {}
for i_col, col_name in enumerate(cols_data):
     cols_rename[col_name] = cols_goodnames[i_col]
data_cata.rename(columns=cols_rename, inplace=True)

## the analytic magnitude distribution
###  per square degree
def logN_mrFunc(mr):
    return -8.85 + 0.71*mr - 0.008*mr**2.

# ++++++++++++++ plot function
def LinePlotFunc_subplots(outpath, N_plots,
                            xvals_list, yvals_list,
                            COLORs_list, LABELs_list=None, LINEs_list=None, LINEWs_list=None, POINTs_list=None, POINTSs_list=None, fillstyles_list=None,
                            subLABEL_list=None, subLABEL_locX=0.1, subLABEL_locY=0.8,
                            XRANGE=None, YRANGE=None,
                            XLABEL=None, YLABEL=None, TITLE=None,
                            xtick_min_label=True, xtick_spe=None, ytick_min_label=True, ytick_spe=None,
                            vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                            hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                            xlog=False, invertX=False, ylog=False, invertY=False, loc_legend='best',
                            font_size=12, usetex=False):
    """
    Line plot for multiple subplots
    """

    # font size
    plt.rc('font', size=font_size)
    # tex
    plt.rcParams["text.usetex"] = usetex

    if outpath != 'show':
        backend_orig = plt.get_backend()
        plt.switch_backend("agg")

    N_rows = math.ceil(N_plots**0.5)
    N_cols = math.ceil(N_plots/N_rows)
    fig, axs = plt.subplots(N_rows, N_cols, sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)

    i_plot = 0
    handles = []
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

                xvals = xvals_list[i_plot]
                yvals = yvals_list[i_plot]

                COLORs = COLORs_list[i_plot]

                if LABELs_list is not None:
                    LABELs = LABELs_list[i_plot]
                else:
                    LABELs = None

                if LINEs_list is not None:
                    LINEs = LINEs_list[i_plot]
                else:
                    LINEs = None
                if LINEWs_list is not None:
                    LINEWs = LINEWs_list[i_plot]
                else:
                    LINEWs = None

                if POINTs_list is not None:
                    POINTs = POINTs_list[i_plot]
                else:
                    POINTs = None
                if POINTSs_list is not None:
                    POINTSs = POINTSs_list[i_plot]
                else:
                    POINTSs = None
                if fillstyles_list is not None:
                    fillstyles = fillstyles_list[i_plot]
                else:
                    fillstyles = None

                for i, xvl in enumerate(xvals):
                    yvl = yvals[i]

                    CR = COLORs[i]

                    if LABELs is not None:
                        LAB = LABELs[i]
                    else:
                        LAB = None

                    if LINEs is not None:
                        LN = LINEs[i]
                    else:
                        LN = '--'
                    if LINEWs is not None:
                        LW = LINEWs[i]
                    else:
                        LW = 1

                    if POINTs is not None:
                        PI = POINTs[i]
                    else:
                        PI = 'o'
                    if POINTSs is not None:
                        MS = POINTSs[i]
                    else:
                        MS = 2
                    if fillstyles is not None:
                        fillstyle = fillstyles[i]
                    else:
                        fillstyle = 'full'

                    tmp = ax.plot(xvl, yvl, color=CR, label=LAB, linestyle=LN, linewidth=LW, marker=PI, markersize=MS, fillstyle=fillstyle)
                    if i_plot == 0:
                        handles.append(tmp[0])

                if subLABEL_list is not None:
                    LABEL = subLABEL_list[i_plot]
                    ax.text(subLABEL_locX, subLABEL_locY, LABEL, transform=ax.transAxes)

                if XRANGE is not None:
                    ax.set_xlim(XRANGE[0], XRANGE[1])
                if YRANGE is not None:
                    ax.set_ylim(YRANGE[0], YRANGE[1])

                if xlog:
                    ax.set_xscale('log')
                if ylog:
                    ax.set_yscale('log')

                if xtick_min_label:
                    if xlog:
                        ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=None, numticks=10))
                    else:
                        ax.xaxis.set_minor_locator(AutoMinorLocator())
                if ytick_min_label:
                    if ylog:
                        ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=None, numticks=10))
                    else:
                        ax.yaxis.set_minor_locator(AutoMinorLocator())

                if xtick_spe is not None:
                    plt.xticks(xtick_spe[0], xtick_spe[1])
                if ytick_spe is not None:
                    plt.yticks(ytick_spe[0], ytick_spe[1])

                if invertY:
                    plt.gca().invert_yaxis()
                if invertX:
                    plt.gca().invert_xaxis()

            i_plot +=1

    fig.text(0.5, 0.02, XLABEL, ha='center')
    fig.text(0.02, 0.5, YLABEL, va='center', rotation='vertical')

    if (LABELs is not None):
        fig.legend(handles, LABELs, loc = 'upper right', bbox_to_anchor=(0.92, 0.35), fancybox=True, shadow=True, fontsize=11)

    if TITLE is not None:
        fig.text(0.5, 0.90, TITLE, ha='center')

    if outpath == 'show':
        plt.show()
        plt.close()
    else:
        plt.savefig(outpath, dpi=300)
        plt.close()
        plt.switch_backend(backend_orig)
        print("Line plot saved as", outpath)

# ++++++++++++++ plot

bands = ['u', 'r', 'i', 'Z', 'Y', 'J', 'H', "Ks"]
bands_label = [r'$u$', r'$r$', r'$i$', r'$Z$', r'$Y$',  r'$J$', r'$H$', r'$K_s$']

XRANGE = (19, 27)
YRANGE = (2e2, 5e4)
dmag = 0.1
nbins = int((XRANGE[1]-XRANGE[0])/dmag)

outpath = 'show'
outpath = './plots/photo_in_skills_cosmos.pdf'

XLABEL = r'Input magnitude'
YLABEL = r'Counts per deg$^2$ per $0.1$ mag'

xvals_list = []
yvals_list = []
COLORs_list = []
LABELs_list = []
LINEs_list = []
LINEWs_list = []
POINTs_list = []

subLABEL_list = []
for i_band, band in enumerate(bands):
    subLABEL_list.append(bands_label[i_band])

    # data
    yval_data, bin_edges = np.histogram(data_cata[f'mag_{band}_app'], bins=nbins, range=XRANGE, weights=np.full(len(data_cata), 1./area_data))

    # simulation original
    yval_sim_ori, _ = np.histogram(sim_cata[f'mag_{band}_app'], bins=nbins, range=XRANGE, weights=np.full(len(sim_cata), 1./area_sim))

    # simulation calibrated
    yval_sim, _ = np.histogram(sim_cata[f'mag_{band}_app']+sim_cata['Dmag_corr'], bins=nbins, range=XRANGE, weights=np.full(len(sim_cata), 1./area_sim))

    xval = (bin_edges[1:] + bin_edges[:-1])/2.
    # analytic
    if band == 'r':
        yval_ana = logN_mrFunc(xval)
    else:
        yval_ana = np.full(len(xval), -999)

    xvals_list.append([xval, xval, xval, xval])
    yvals_list.append([yval_sim_ori, yval_sim, yval_data, np.power(10., yval_ana)])
    COLORs_list.append(['green', 'blue', 'r', 'k'])
    LABELs_list.append(['Shark (original)', 'Shark (matched)', 'COSMOS2015', 'KiDS fitting'])

    LINEs_list.append(['--', '-', '-', '--'])
    LINEWs_list.append([2, 2, 2, 2])
    POINTs_list.append(['', '', '', ''])

LinePlotFunc_subplots(outpath, len(bands),
                            xvals_list, yvals_list,
                            COLORs_list, LABELs_list=LABELs_list, LINEs_list=LINEs_list, LINEWs_list=LINEWs_list, POINTs_list=POINTs_list, POINTSs_list=None, fillstyles_list=None,
                            subLABEL_list=subLABEL_list, subLABEL_locX=0.1, subLABEL_locY=0.75,
                            XRANGE=XRANGE, YRANGE=YRANGE,
                            XLABEL=XLABEL, YLABEL=YLABEL, TITLE=None,
                            xtick_min_label=True, xtick_spe=None, ytick_min_label=True, ytick_spe=None,
                            vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                            hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                            xlog=False, invertX=False, ylog=True, invertY=False, loc_legend='best',
                            font_size=16, usetex=True)