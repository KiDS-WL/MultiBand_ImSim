# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-12-21 15:10:57
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-29 19:56:58

### Appendix A: stellar mass to light ratio

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
# original SURFS cata
area_surfs = 108.
infile_surfs = '/disks/shear10/ssli/ImSim/input/SURFS_cata/ori/SURFS_SHARK_LC_100sqdeg_testcols.fits'
with fits.open(infile_surfs) as hdul:
    surfs_cata = hdul[1].data
print('surfs ', len(surfs_cata))

# cosmos catalogue
area_cosmos = 1.38
infile_cosmos = '/disks/shear10/ssli/ImSim/input/COSMOS_cata/outside/cosmos_2015_dr_2_1_138deg2_galONLY.fits'
with fits.open(infile_cosmos) as hdul:
    cosmos_cata = hdul[1].data
print('cosmos ', len(cosmos_cata))

# ++++++++++++++ the function
def ErrorPlotFunc_subplots(outpath, N_plots,
                            xvals_list, yvals_list, yerrs_list,
                            COLORs_list, LABELs_list=None, LINEs_list=None, LINEWs_list=None, POINTs_list=None, POINTSs_list=None, ERRORSIZEs_list=None,
                            subLABEL_list=None, subLABEL_locX=0.1, subLABEL_locY=0.8,
                            XRANGE=None, YRANGE=None,
                            XLABEL=None, YLABEL=None, TITLE=None,
                            xtick_min_label=True, xtick_spe=None, ytick_min_label=True, ytick_spe=None,
                            vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                            hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                            xlog=False, invertX=False, ylog=False, invertY=False, loc_legend='best',
                            font_size=12, usetex=False):
    """
    Errorbar plot for multiple subplots
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
                if yerrs_list is not None:
                    yerrs = yerrs_list[i_plot]
                else:
                    yerrs = None

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
                if ERRORSIZEs_list is not None:
                    ERRORSIZEs = ERRORSIZEs_list[i_plot]
                else:
                    ERRORSIZEs = None

                for i, xvl in enumerate(xvals):
                    yvl = yvals[i]
                    if yerrs is not None:
                        yerr = yerrs[i]
                    else:
                        yerr = None
                    if yerr is not None:
                        yerr = np.array(yerr)
                        yerr = np.vstack([yerr[0], yerr[1]])

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

                    if ERRORSIZEs is not None:
                        ERRORSIZE = ERRORSIZEs[i]
                    else:
                        ERRORSIZE = 2

                    tmp = ax.errorbar(xvl, yvl, yerr=yerr, color=CR, label=LAB, linestyle=LN, linewidth=LW, marker=PI, markersize=MS, capsize=ERRORSIZE)
                    if i_plot == 0:
                        handles.append(tmp[0])

                if subLABEL_list is not None:
                    LABEL = subLABEL_list[i_plot]
                    ax.text(subLABEL_locX, subLABEL_locY, LABEL, transform=ax.transAxes, fontsize=14)

                if XRANGE is not None:
                    ax.set_xlim(XRANGE[0], XRANGE[1])
                if YRANGE is not None:
                    ax.set_ylim(YRANGE[0], YRANGE[1])

                if xlog:
                    ax.set_xscale('log')
                if ylog:
                    ax.set_yscale('log')

                if vlines is not None:
                    _vhlines('v', vlines, line_styles=vline_styles, line_colors=vline_colors, line_labels=vline_labels, line_widths=vline_widths, ax=ax)
                if hlines is not None:
                    _vhlines('h', hlines, line_styles=hline_styles, line_colors=hline_colors, line_labels=hline_labels, line_widths=hline_widths, ax=ax)

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
        fig.legend(handles, LABELs, loc = 'upper right', bbox_to_anchor=(0.908, 0.30), fancybox=True, shadow=True, fontsize=12)

    if TITLE is not None:
        fig.text(0.5, 0.90, TITLE, ha='center')

    if outpath == 'show':
        plt.show()
        plt.close()
    else:
        plt.savefig(outpath, dpi=300)
        plt.close()
        plt.switch_backend(backend_orig)
        print("Errorbar plot saved as", outpath)

# ++++++++++++++ run

############# Ks band
mag_sun = 5.08

# log stellar mass
logmass_cosmos = cosmos_cata[f'mass_med'] # log stellar mass
logmass_surfs = np.log10((surfs_cata['mstars_bulge'] + surfs_cata['mstars_disk'])/0.7)

# luminosity
logL_cosmos = -0.4*(cosmos_cata[f'm_k']-mag_sun) # solar L
logL_surfs = -0.4*(surfs_cata[f'K_VISTA_absolute']-mag_sun)

# mass to light ratio
Upsilon_cosmos = 10**(logmass_cosmos - logL_cosmos)
Upsilon_surfs = 10**(logmass_surfs - logL_surfs)

x_cosmos, x_surfs = [logmass_cosmos, logmass_surfs]
y_cosmos, y_surfs = [Upsilon_cosmos, Upsilon_surfs]

XRANGE = [8.5, 11.5]
YRANGE = [0., 0.8]
ytick_spe = [[0.2, 0.4, 0.6], [r'$0.2$', r'$0.4$', r'$0.6$']]

# nboot = 500

# outpath = 'show'
outpath = f'./plots/mass_to_light_Ks.pdf'

XLABEL = r'${\rm log}_{10}~({\rm M}_{\star}/{\rm M}_{\odot})$'
YLABEL = r'${\rm M}_{\star}/{\rm L}_{K_{\rm s}}~[{\rm M}_{\odot}/{\rm L}_{\odot}]$'

xvals_list = []
yvals_list = []
yerrs_list = None
COLORs_list = []
LABELs_list = []
LINEs_list = []
LINEWs_list = []
POINTs_list = []

subLABEL_list = []

zbin_edges = np.arange(0.2, 2.5, 0.2)
N_plots = 8
for i_zbin in range(N_plots):
    zbin_min = zbin_edges[i_zbin]
    zbin_max = zbin_edges[i_zbin+1]

    label = f'${zbin_min:.1f}\\leq z<{zbin_max:.1f}$'
    subLABEL_list.append(label)

    # cosmos
    mask_tmp = (cosmos_cata['photoz']>=zbin_min) & (cosmos_cata['photoz']<zbin_max)
    x_cosmos_selec = x_cosmos[mask_tmp]
    y_cosmos_selec = y_cosmos[mask_tmp]

    # surfs
    mask_tmp = (surfs_cata['zobs']>=zbin_min) & (surfs_cata['zobs']<zbin_max)
    x_surfs_selec = x_surfs[mask_tmp]
    y_surfs_selec = y_surfs[mask_tmp]

    xval_combine = [x_cosmos_selec, x_surfs_selec]
    yval_combine = [y_cosmos_selec, y_surfs_selec]

    # binning
    xvals = []
    yvals = []
    # errs = []
    edges_x_bins = np.arange(XRANGE[0], XRANGE[1], 0.2)
    for i_val, x_val in enumerate(xval_combine):
        y_val = yval_combine[i_val]

        x_list = []
        median_list = []
        # err_low_list = []
        # err_high_list = []
        for i_xbin in range(len(edges_x_bins)-1):
            x_min = edges_x_bins[i_xbin]
            x_max = edges_x_bins[i_xbin+1]
            mask_x = (x_val >= x_min) & (x_val < x_max)
            y_selec = y_val[mask_x]

            if len(y_selec)>10:
                y_median = np.median(y_selec)
                median_list.append(y_median)

                # ## get error from boots
                # Relications_index = np.array([np.random.randint(len(y_selec), size=len(y_selec)) for _ in range(nboot)])
                # med_BS = np.median(y_selec[Relications_index], axis = 1)
                # sigma_BS = med_BS - y_median
                # err_low = np.abs(np.percentile(sigma_BS, 16))
                # err_high = np.abs(np.percentile(sigma_BS, 84))
                # ## assign
                # err_low_list.append(err_low)
                # err_high_list.append(err_high)

                x_list.append((x_min+x_max)/2.)

            else:
                continue

        xvals.append(x_list)
        yvals.append(median_list)
        # errs.append([err_low_list, err_high_list])

    xvals_list.append(xvals)
    yvals_list.append(yvals)
    # yerrs_list.append(errs)

    COLORs_list.append(['r', 'green'])
    LABELs_list.append(['COSMOS2015', 'Shark'])

    LINEs_list.append(['-', '--'])
    LINEWs_list.append([2, 2])
    POINTs_list.append(['', ''])


ErrorPlotFunc_subplots(outpath, N_plots,
                            xvals_list, yvals_list, yerrs_list,
                            COLORs_list, LABELs_list=LABELs_list, LINEs_list=LINEs_list, 
                            LINEWs_list=LINEWs_list, POINTs_list=POINTs_list, POINTSs_list=None, ERRORSIZEs_list=None,
                            subLABEL_list=subLABEL_list, subLABEL_locX=0.08, subLABEL_locY=0.7,
                            XRANGE=XRANGE, YRANGE=YRANGE,
                            XLABEL=XLABEL, YLABEL=YLABEL, TITLE=None,
                            xtick_min_label=True, xtick_spe=None, ytick_min_label=True, ytick_spe=ytick_spe,
                            vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                            hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                            xlog=False, invertX=False, ylog=False, invertY=False, loc_legend='best',
                            font_size=16, usetex=True)

