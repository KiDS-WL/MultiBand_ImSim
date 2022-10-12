# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-06-28 13:38:35
# @Last Modified by:   lshuns
# @Last Modified time: 2022-08-25 12:18:48

### compare the GAaP limiting magnitude

import os
import sys
import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.io import fits

import plotting

# ++++++++++++++++++++++++++ I/O

# simulation 
infile = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_photo_noSG.feather'
sim_whole = pd.read_feather(infile)

# data
infile = '/disks/shear16/ssli/KiDS/K1000_forSKiLLS/skills_fiducial_kids_photo_noSG.feather'
data_whole = pd.read_feather(infile)

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

# plot directory
outpath = 'show'
# outpath = './plots/gaap_limiting.pdf'

# ++++++++++++++++++++++++++ workhorse

# collect points
yvals_list = []
mag_auto_list = []
for tile_label in tile_labels:

    # +++++++++++++++ observation
    data = data_whole[data_whole['THELI_NAME']== 'KIDS_'+tile_label.replace('.', 'p').replace('-', 'm')]
    ## zeropoints
    zero_ugri = dmag_cata[dmag_cata.label==tile_label]

    # ++++++++++++++++ simulation
    dataSim = sim_whole[sim_whole['tile_label']==tile_label]

    para_sims = []
    para_reals = []
    # ++++++++++++++++ calculation
    for band in bands:

        # observation
        flux_tmp = data['FLUXERR_GAAP_{:}'.format(band)]
        if band in ['u', 'g', 'r', 'i']:
            zero = zero_ugri['dmag_{:}'.format(band)].values[0]
            gaap_tmp = zero - 2.5*np.log10(flux_tmp)
        else:
            gaap_tmp = 30 - 2.5*np.log10(flux_tmp)
        para_reals.append(gaap_tmp)

        # simulation
        flux_tmp = dataSim['FLUXERR_GAAP_{:}'.format(band)]
        gaap_tmp = 30 - 2.5*np.log10(flux_tmp)
        para_sims.append(gaap_tmp)

    # median difference as y
    yvals = np.nanmedian(para_sims, axis=1) - np.nanmedian(para_reals, axis=1)
    yvals_list.append(yvals)

    del data, dataSim

# get y values
yvals_median = np.nanmedian(yvals_list, axis=0)
yvals_low = np.nanpercentile(yvals_list, 16, axis=0)
yvals_high = np.nanpercentile(yvals_list, 84, axis=0)

## for plot
y_forplot = [yvals_median, yvals_low, yvals_high]
yerr_forplot = None

LABELs = ['50', '16', '84']
LINEs = ['-', '--', ':']

# get x values
bands = [r'$u$', r'$g$', r'$r$', r'$i$', r'$Z$', r'$Y$',  r'$J$', r'$H$', r'$K_s$']
x_forplot = [np.arange(1,len(bands)+1, 1)] * 3
xtick_spe = [np.arange(1,len(bands)+1,1), bands]

# plot related
COLORs = ['k'] * 3
XLABEL = 'Band'
YLABEL = r'$\Delta~1\sigma$ \textsc{GAaP} limiting magnitude'
YRANGE = [-0.3, 0.3]
hlines = [0.0]
hline_styles = ['dotted']
hline_colors = ['gray']

xtick_min_label = False

plotting.ErrorPlotFunc(outpath,
                x_forplot, y_forplot, yerr_forplot,
                COLORs, LABELs=LABELs, LINEs=LINEs, LINEWs=None, POINTs=None, POINTSs=None, ERRORSIZEs=None,
                XRANGE=None, YRANGE=YRANGE,
                XLABEL=XLABEL, YLABEL=YLABEL, TITLE=None,
                xtick_min_label=xtick_min_label, xtick_spe=xtick_spe, ytick_min_label=True, ytick_spe=None,
                vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                hlines=hlines, hline_styles=hline_styles, hline_colors=hline_colors, hline_labels=None, hline_widths=None,
                xlog=False, invertX=False, ylog=False, invertY=False, loc_legend='best',
                font_size=16, usetex=True)