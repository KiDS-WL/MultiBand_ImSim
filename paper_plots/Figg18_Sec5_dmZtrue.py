# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-09-04 21:30:42
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-30 11:44:56

#### how Dm VS zTrue look like

import os
import pandas as pd 
import numpy as np

import plotting

# the catalogues
inpath_dm_list = [f'../correction_varShear_PSFmodelling/outputs/dm_forBlending_dz0p1_ZtrueBins_ZBbin{i_bin}.csv' for i_bin in range(6)]

subLABEL_list = [r'$0.1< z_{\rm B} \leq 0.3$', 
                r'$0.3< z_{\rm B} \leq 0.5$', 
                r'$0.5< z_{\rm B} \leq 0.7$', 
                r'$0.7< z_{\rm B} \leq 0.9$', 
                r'$0.9< z_{\rm B} \leq 1.2$',
                r'$1.2< z_{\rm B} \leq 2.0$']

# plot related
outpath = './plots/dmZtrue.pdf'
# outpath = 'show'

N_rows = 2
N_cols = 3
FIGSIZE = [8, 4]

XLABEL = 'True redshift'
YLABEL = r'$\Delta m_{\rm blending}$ + 1'
XRANGE = [-0.05, 2.55]
YRANGE = [0, 1.5]
hlines = [1]
hline_styles = [':']
hline_colors = ['gray']

font_size = 12

# get values
N_plots = 6
xvals_list = []
yvals_list = []
yerrs_list = []
COLORs_list = []
LINEs_list = []
LINEWs_list = []
POINTs_list = []
POINTSs_list = []
ERRORSIZEs_list = []
fill_between_xs_list = None
fill_between_yLows_list = None
fill_between_yHighs_list = None
fill_between_COLORs_list = None
fill_between_alphas_list = None
for i_bin in range(N_plots):

    inpath = inpath_dm_list[i_bin]
    cata = pd.read_csv(inpath)

    # select with measurements
    cata = cata[cata['TotalWei']>0]

    # z true as the x
    xval = cata['zTrueMean'].values

    # normalise the weights
    y_wei = cata['TotalWei'].values / np.max(cata['TotalWei'].values)

    # 1+m as y
    yval = (cata['m1'].values + cata['m2'].values)/2. + 1
    yerr = (cata['m1_err'].values + cata['m2_err'].values)/2.

    # plot
    xvals_list.append([xval, xval])
    yvals_list.append([yval, y_wei])
    yerrs_list.append([[yerr, yerr], None])

    COLORs_list.append(['orange', 'k'])
    LINEWs_list.append([1, 1.5])
    POINTSs_list.append([3, 0])
    ERRORSIZEs_list.append([1.5, 0])
    LINEs_list.append(['', '--'])
    POINTs_list.append(['o', ''])

plotting.ErrorPlotFunc_subplots(outpath, N_plots,
                            xvals_list, yvals_list, yerrs_list,
                            COLORs_list, LABELs_list=None, LINEs_list=LINEs_list, LINEWs_list=LINEWs_list, 
                            POINTs_list=POINTs_list, POINTSs_list=POINTSs_list, ERRORSIZEs_list=ERRORSIZEs_list,
                            subLABEL_list=subLABEL_list, subLABEL_locX=0.45, subLABEL_locY=0.85,
                            XRANGE=XRANGE, YRANGE=YRANGE,
                            XLABEL=XLABEL, YLABEL=YLABEL, TITLE=None,
                            xtick_min_label=True, xtick_spe=None, ytick_min_label=True, ytick_spe=None,
                            vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                            hlines=hlines, hline_styles=hline_styles, hline_colors=hline_colors, hline_labels=None, hline_widths=None,
                            xlog=False, invertX=False, ylog=False, invertY=False, 
                            loc_legend='best', legend_frame=False,
                            font_size=font_size, usetex=True,
                            fill_between_xs_list=fill_between_xs_list, 
                            fill_between_yLows_list=fill_between_yLows_list, fill_between_yHighs_list=fill_between_yHighs_list,
                            fill_between_COLORs_list=fill_between_COLORs_list, fill_between_alphas_list=fill_between_alphas_list,
                            LABEL_position='inSub', LABEL_position_SUBid=0,
                            LABEL_cols=1,
                            FIGSIZE=FIGSIZE,
                            TIGHT=True, 
                            N_rows=N_rows, N_cols=N_cols)