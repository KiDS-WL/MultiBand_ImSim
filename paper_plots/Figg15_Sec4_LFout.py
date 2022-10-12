# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-06-06 16:52:28
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-12 15:19:06

### plot results from LF measurements

import pandas as pd 
import numpy as np

import plotting

# ++++++++++++++++++++++++++ I/O

outprefix = 'show'
outprefix = './plots/LFcomp'

# data
infile = '/disks/shear16/ssli/KiDS/K1000_forSKiLLS/kids_photo_LF_321_shear_noSG_noWeiCut_newCut_817tiles.feather'
data_whole = pd.read_feather(infile)
## only weight > 0
data_whole = data_whole[data_whole['weight']>0]

# simulation 
infile = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_shear_noSG_noWeiCut_newCut.feather'
sim_whole = pd.read_feather(infile)
## only weight > 0
sim_whole = sim_whole[sim_whole['oldweight_LF_r']>0]

########## general 
font_size = 20

# ++++++++++++++++++++++++++ plot SNR
if outprefix == 'show':
    outpath = 'show'
else:
    outpath = outprefix + '_SNR.pdf'

paras = [data_whole['model_SNratio'].values, sim_whole['SNR_LF_r'].values]
wgs = [data_whole['weight'].values, sim_whole['oldweight_LF_r'].values]
COLORs = ['r', 'b']
LABELs = ['KiDS DR4', 'SKiLLS fiducial']
nbins = 60
XRANGE = [5, 200]

XLABEL = r'$\nu_{\rm SN}$'
YLABEL = r'\textit{lens}fit-weighted distribution'
DENSITY = True

xlog = True

plotting.HistPlotFunc(outpath,
                paras, wgs, COLORs, LABELs,
                nbins, XRANGE, YRANGE=None,
                XLABEL=XLABEL, YLABEL=YLABEL,
                DENSITY=DENSITY, HISTTYPE='step', STACKED=False,
                TITLE=None, xtick_min_label=True, ytick_min_label=True,
                xtick_spe=None, ytick_spe=None,
                vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                xlog=xlog, ylog=False,
                loc_legend='best', 
                font_size=font_size, usetex=True,
                cumulative=False)

# ++++++++++++++++++++++++++ plot LF weights
if outprefix == 'show':
    outpath = 'show'
else:
    outpath = outprefix + '_weight.pdf'

paras = [data_whole['weight'].values, sim_whole['oldweight_LF_r'].values]
wgs = None
COLORs = ['r', 'b']
LABELs = None
nbins = 60
XRANGE = [0, np.max(data_whole['weight'].values)+0.01]

XLABEL = r'\textit{lens}fit weight'
YLABEL = r'Distribution'
DENSITY = True

xlog = False

plotting.HistPlotFunc(outpath,
                paras, wgs, COLORs, LABELs,
                nbins, XRANGE, YRANGE=None,
                XLABEL=XLABEL, YLABEL=YLABEL,
                DENSITY=DENSITY, HISTTYPE='step', STACKED=False,
                TITLE=None, xtick_min_label=True, ytick_min_label=True,
                xtick_spe=None, ytick_spe=None,
                vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                xlog=xlog, ylog=False,
                loc_legend='best', 
                font_size=font_size, usetex=True,
                cumulative=False)

# ++++++++++++++++++++++++++ plot scalelength
if outprefix == 'show':
    outpath = 'show'
else:
    outpath = outprefix + '_size.pdf'

paras = [data_whole['autocal_scalelength_pixels'].values * 0.214, sim_whole['scalelength_LF_r'].values * 0.214]
wgs = [data_whole['weight'].values, sim_whole['oldweight_LF_r'].values]
COLORs = ['r', 'b']
LABELs = None
nbins = 60
XRANGE = [0.08, 1.5]

XLABEL = 'Scalelength [arcsec]'
YLABEL = r'\textit{lens}fit-weighted distribution'
DENSITY = True

xlog = True

plotting.HistPlotFunc(outpath,
                paras, wgs, COLORs, LABELs,
                nbins, XRANGE, YRANGE=None,
                XLABEL=XLABEL, YLABEL=YLABEL,
                DENSITY=DENSITY, HISTTYPE='step', STACKED=False,
                TITLE=None, xtick_min_label=True, ytick_min_label=True,
                xtick_spe=None, ytick_spe=None,
                vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                xlog=xlog, ylog=False,
                loc_legend='best', 
                font_size=font_size, usetex=True,
                cumulative=False)

# ++++++++++++++++++++++++++ plot Ellipticity
if outprefix == 'show':
    outpath = 'show'
else:
    outpath = outprefix + '_e.pdf'

paras = [np.hypot(data_whole['autocal_e1'].values, data_whole['autocal_e2'].values), 
            np.hypot(sim_whole['e1_LF_r'].values, sim_whole['e2_LF_r'].values)]
wgs = [data_whole['weight'].values, sim_whole['oldweight_LF_r'].values]
COLORs = ['r', 'b']
LABELs = None
nbins = 60
XRANGE = [0., 0.9]

XLABEL = 'Ellipticity'
YLABEL = r'\textit{lens}fit-weighted distribution'
DENSITY = True

xlog = False

plotting.HistPlotFunc(outpath,
                paras, wgs, COLORs, LABELs,
                nbins, XRANGE, YRANGE=None,
                XLABEL=XLABEL, YLABEL=YLABEL,
                DENSITY=DENSITY, HISTTYPE='step', STACKED=False,
                TITLE=None, xtick_min_label=True, ytick_min_label=True,
                xtick_spe=None, ytick_spe=None,
                vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                xlog=xlog, ylog=False,
                loc_legend='best', 
                font_size=font_size, usetex=True,
                cumulative=False)

# ++++++++++++++++++++++++++ plot neighbour magnitude
if outprefix == 'show':
    outpath = 'show'
else:
    outpath = outprefix + '_nm.pdf'

paras = [data_whole['neighbour_mag'].values, sim_whole['nm_LF_r'].values]
wgs = [data_whole['weight'].values, sim_whole['oldweight_LF_r'].values]
COLORs = ['r', 'b']
LABELs = None
nbins = 60
XRANGE = [-6, 5]

XLABEL = 'Neighbour magnitude difference'
YLABEL = r'\textit{lens}fit-weighted distribution'
DENSITY = True

xlog = False

plotting.HistPlotFunc(outpath,
                paras, wgs, COLORs, LABELs,
                nbins, XRANGE, YRANGE=None,
                XLABEL=XLABEL, YLABEL=YLABEL,
                DENSITY=DENSITY, HISTTYPE='step', STACKED=False,
                TITLE=None, xtick_min_label=True, ytick_min_label=True,
                xtick_spe=None, ytick_spe=None,
                vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                xlog=xlog, ylog=False,
                loc_legend='best', 
                font_size=font_size, usetex=True,
                cumulative=False)

# ++++++++++++++++++++++++++ plot neighbour distance
if outprefix == 'show':
    outpath = 'show'
else:
    outpath = outprefix + '_nd.pdf'

paras = [data_whole['neighbour_distance'].values, sim_whole['nd_LF_r'].values]
wgs = [data_whole['weight'].values, sim_whole['oldweight_LF_r'].values]
COLORs = ['r', 'b']
LABELs = None
nbins = 60
XRANGE = [0.8, 18]

XLABEL = 'Neighbour distance [arcsec]'
YLABEL = r'\textit{lens}fit-weighted distribution'
DENSITY = True

xlog = True

plotting.HistPlotFunc(outpath,
                paras, wgs, COLORs, LABELs,
                nbins, XRANGE, YRANGE=None,
                XLABEL=XLABEL, YLABEL=YLABEL,
                DENSITY=DENSITY, HISTTYPE='step', STACKED=False,
                TITLE=None, xtick_min_label=True, ytick_min_label=True,
                xtick_spe=None, ytick_spe=None,
                vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                xlog=xlog, ylog=False,
                loc_legend='best', 
                font_size=font_size, usetex=True,
                cumulative=False)
