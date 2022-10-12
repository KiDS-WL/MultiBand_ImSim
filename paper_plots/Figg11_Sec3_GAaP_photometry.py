# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-01-22 15:27:31
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-12 15:15:22

### GAaP colour and magnitudes

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
from matplotlib.patches import Rectangle

import plotting

# ++++++++++++++++++++++++++ I/O

# simulation 
infile = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_photo_noSG_shear_m283m283_rot_0.feather'
# infile = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_photo_noSG_shear_m283m283_rot_0.feather'
sim_whole = pd.read_feather(infile)
## get galaxy only
sim_whole = sim_whole[sim_whole['id_input']>0]
sim_whole.reset_index(drop=True, inplace=True)
## sky area
area_sim = 108.

# data
infile = '/disks/shear16/ssli/KiDS/K1000_forSKiLLS/kids_photo_LF_321_photo_noSG_817tiles.feather'
data_whole = pd.read_feather(infile)
## get galaxy only
data_whole = data_whole[(data_whole['SG_FLAG']==1) & (data_whole['SG2DPHOT']==0)]
data_whole.reset_index(drop=True, inplace=True)
## sky area
area_data = 636.2882737918676

# tile labels
infile = '../noise_info/skills_fiducial/noise_selec_combined.csv'
tile_selec = pd.read_csv(infile)
tile_labels = tile_selec['label'].to_list()
del tile_selec

# zero point for observation
infile = '../noise_info/kids_dr4_dmag_ugri.csv'
dmag_cata = pd.read_csv(infile)

# ++++++++++++++++++++++++++ general setting

bands = ['u', 'g', 'r', 'i', 'Z', 'Y', 'J', 'H', 'Ks']
bands_label = [r'$u$', r'$g$', r'$r$', r'$i$', r'$Z$', r'$Y$',  r'$J$', r'$H$', r'$K_s$']

# ++++++++++++++ plot function
def HistPlotFunc_subplots(outpath, N_plots,
                            paras_list, wgs_list, COLORs_list, LABELs_list,
                            nbins_list, XRANGE, YRANGE=None,
                            subLABEL_list=None, subLABEL_locX=0.1, subLABEL_locY=0.8,
                            XLABEL=None, YLABEL=None,
                            DENSITY=False, HISTTYPE='step', STACKED=False,
                            TITLE=None, xtick_min_label=True, ytick_min_label=True,
                            xtick_spe=None, ytick_spe=None,
                            vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                            hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                            xlog=False, ylog=False,
                            loc_legend='best',
                            font_size=12, usetex=False,
                            LABELs_loc='top'):
    """
    Histogram plot for multiple subplots
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

    if DENSITY and (wgs_list is not None):
        logger.warning('DENSITY and wgs are provided simultaneously!!!')

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

                paras = paras_list[i_plot]
                COLORs = COLORs_list[i_plot]
                if LABELs_list is not None:
                    LABELs = LABELs_list[i_plot]
                else:
                    LABELs = None

                nbins = nbins_list[i_plot]
                if wgs_list is not None:
                    wgs = wgs_list[i_plot]
                else:
                    wgs = None

                if xlog:
                    logbins = np.logspace(np.log10(XRANGE[0]), np.log10(XRANGE[1]), nbins)
                    ax.hist(x=paras, bins=logbins, density=DENSITY, weights=wgs, color=COLORs, label=LABELs, histtype=HISTTYPE, stacked=STACKED)
                else:
                    ax.hist(x=paras, bins=nbins, range=XRANGE, density=DENSITY, weights=wgs, color=COLORs, label=LABELs, histtype=HISTTYPE, stacked=STACKED)

                if subLABEL_list is not None:
                    LABEL = subLABEL_list[i_plot]
                    ax.text(subLABEL_locX, subLABEL_locY, LABEL, transform=ax.transAxes)

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

            i_plot +=1

    fig.text(0.5, 0.02, XLABEL, ha='center')
    fig.text(0.02, 0.5, YLABEL, va='center', rotation='vertical')

    if (LABELs is not None):

        handles = [Rectangle((0,0),1,1,color='white', ec=c) for c in COLORs]

        if LABELs_loc == 'top':
            fig.legend(handles, LABELs, 
                    loc = 'center', ncol=2,
                    bbox_to_anchor=(0.5, 0.95), fancybox=True, shadow=True, fontsize=12)

        elif LABELs_loc == 'right':
            fig.legend(handles, LABELs, 
                loc = 'upper right', ncol=1,
                bbox_to_anchor=(0.92, 0.30), fancybox=True, shadow=True, fontsize=12)

    if TITLE is not None:
        fig.text(0.5, 0.90, TITLE, ha='center')

    if outpath == 'show':
        plt.show()
        plt.close()
    else:
        plt.savefig(outpath, dpi=300)
        plt.close()
        plt.switch_backend(backend_orig)
        print("Histogram plot saved as", outpath)

######### magnitude

XRANGE = [18., 26]
dmag = 0.1
nbins = int((XRANGE[1]-XRANGE[0])/dmag)

outpath = './plots/gaap_mag.pdf'
# outpath = 'show'
N_plots = len(bands)
paras_list = []
wgs_list = []
COLORs_list = []
LABELs_list = []
nbins_list = []
subLABEL_list = []
for i_plot in range(N_plots):

    band = bands[i_plot]
    band_label = bands_label[i_plot]

    subLABEL_list.append(band_label)

    paras_list.append([data_whole[f'MAG_GAAP_{band}'], 
                        sim_whole[f'MAG_GAAP_{band}']])
    COLORs_list.append(['r', 'b'])
    LABELs_list.append(['KiDS DR4', 'SKiLLS fiducial'])
    wgs_list.append([np.ones(len(data_whole))/area_data, np.ones(len(sim_whole))/area_sim])
    nbins_list.append(nbins)

xtick_spe = [[20, 22, 24], ['20', '22', '24']]

YRANGE = [1e2, 7e3]
ytick_spe = None

XLABEL = r'\textsc{GAaP} magnitude'
YLABEL = r'Counts per deg$^2$ per $0.1$ mag'
DENSITY = False
ylog = True

HistPlotFunc_subplots(outpath, N_plots,
                            paras_list, wgs_list, COLORs_list, LABELs_list,
                            nbins_list, XRANGE, YRANGE=YRANGE,
                            subLABEL_list=subLABEL_list, subLABEL_locX=0.1, subLABEL_locY=0.75,
                            XLABEL=XLABEL, YLABEL=YLABEL,
                            DENSITY=DENSITY, HISTTYPE='step', STACKED=False,
                            TITLE=None, xtick_min_label=True, ytick_min_label=True,
                            xtick_spe=xtick_spe, ytick_spe=ytick_spe,
                            vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                            hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                            xlog=False, ylog=ylog,
                            loc_legend='best',
                            font_size=16, usetex=True,
                            LABELs_loc='top')

######### colours

XRANGE = [-2., 2.]
dmag = 0.1
nbins = int((XRANGE[1]-XRANGE[0])/dmag)

outpath = './plots/gaap_colour.pdf'
# outpath = 'show'
N_plots = len(bands) - 1 
paras_list = []
wgs_list = []
COLORs_list = []
LABELs_list = []
nbins_list = []
subLABEL_list = []
for i_plot in range(N_plots):

    band0 = bands[i_plot]
    band1 = bands[i_plot + 1]
    band0_label = bands_label[i_plot]
    band1_label = bands_label[i_plot + 1]

    subLABEL_list.append(band0_label + '-' + band1_label)

    # data
    para_data = data_whole[f'MAG_GAAP_{band0}']-data_whole[f'MAG_GAAP_{band1}']
    para_data = para_data[para_data!=0]

    # simulation
    para_sim = sim_whole[f'MAG_GAAP_{band0}']-sim_whole[f'MAG_GAAP_{band1}']
    para_sim = para_sim[para_sim!=0]

    paras_list.append([para_data, para_sim])
    COLORs_list.append(['r', 'b'])
    LABELs_list.append(['KiDS DR4', 'SKiLLS fiducial'])
    wgs_list.append([np.ones(len(para_data))/area_data, np.ones(len(para_sim))/area_sim])
    nbins_list.append(nbins)

xtick_spe = [[-1, 0, 1], ['-1', '0', '1']]

YRANGE = [0, 12000]
ytick_spe = [[3000, 6000, 9000], 
        ['3', '6', '9']]

XLABEL = r'\textsc{GAaP} colour'
YLABEL = r'Counts per deg$^2$ per $0.1$ mag [$\times 10^{3}$]'
DENSITY = False
ylog = False

HistPlotFunc_subplots(outpath, N_plots,
                            paras_list, wgs_list, COLORs_list, LABELs_list,
                            nbins_list, XRANGE, YRANGE=YRANGE,
                            subLABEL_list=subLABEL_list, subLABEL_locX=0.1, subLABEL_locY=0.75,
                            XLABEL=XLABEL, YLABEL=YLABEL,
                            DENSITY=DENSITY, HISTTYPE='step', STACKED=False,
                            TITLE=None, xtick_min_label=True, ytick_min_label=True,
                            xtick_spe=xtick_spe, ytick_spe=ytick_spe,
                            vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                            hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                            xlog=False, ylog=ylog,
                            loc_legend='best',
                            font_size=16, usetex=True,
                            LABELs_loc='right')
