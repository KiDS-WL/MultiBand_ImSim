# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-07-28 13:58:32
# @Last Modified by:   lshuns
# @Last Modified time: 2021-11-05 15:43:16

### compare stellar density

import numpy as np
import pandas as pd

from matplotlib.pyplot import cm

import plotting

# ++++++++++++++++++++++++ I/O

## COllege
old = '/disks/shear15/KiDS/ImSim/pipeline/utils/bsc_model_175_0_m18_25'
old = np.loadtxt(old).flatten()

## SKiLLS
infile = '/disks/shear10/ssli/ImSim/input/star_cata/trilegal_9band_-25.0_-31.0_10deg2_default.feather'
new25 = pd.read_feather(infile)
infile = '/disks/shear10/ssli/ImSim/input/star_cata/trilegal_9band_10.0_-31.0_10deg2_default.feather'
new10 = pd.read_feather(infile)
infile = '/disks/shear10/ssli/ImSim/input/star_cata/trilegal_9band_45.0_-31.0_10deg2_default.feather'
new45 = pd.read_feather(infile)
infile = '/disks/shear10/ssli/ImSim/input/star_cata/trilegal_9band_135.0_0.0_10deg2_default.feather'
new135 = pd.read_feather(infile)
infile = '/disks/shear10/ssli/ImSim/input/star_cata/trilegal_9band_180.0_0.0_10deg2_default.feather'
new180 = pd.read_feather(infile)
infile = '/disks/shear10/ssli/ImSim/input/star_cata/trilegal_9band_230.0_0.0_10deg2_default.feather'
new230 = pd.read_feather(infile)

# +++++++++++++++ plot

XRANGE = [14.5, 25.5]
YRANGE = [5, 600]

dmag = 0.1
nbins = int((XRANGE[1]-XRANGE[0])/dmag)

# outpath = 'show'
outpath = f'./plots/stars.pdf'
XLABEL = r'Input magnitude ($r$-band)'
YLABEL = r'Counts per deg$^2$ per 0.1 mag'

# old 
yval_old, bin_edges = np.histogram(old, bins=nbins, range=XRANGE, weights=np.full(len(old), 1))
# new 
yval_new25, _ = np.histogram(new25['r'], bins=nbins, range=XRANGE, weights=np.full(len(new25), 1/10.))
# new 
yval_new10, _ = np.histogram(new10['r'], bins=nbins, range=XRANGE, weights=np.full(len(new10), 1/10.))
# new 
yval_new45, _ = np.histogram(new45['r'], bins=nbins, range=XRANGE, weights=np.full(len(new45), 1/10.))
# new 
yval_new135, _ = np.histogram(new135['r'], bins=nbins, range=XRANGE, weights=np.full(len(new135), 1/10.))
# new 
yval_new180, _ = np.histogram(new180['r'], bins=nbins, range=XRANGE, weights=np.full(len(new180), 1/10.))
# new 
yval_new230, _ = np.histogram(new230['r'], bins=nbins, range=XRANGE, weights=np.full(len(new230), 1/10.))

xval = (bin_edges[1:] + bin_edges[:-1])/2.

LABELs = ['COllege', r'$(-25.0, -31.0)$', r'$(10.0, -31.0)$', r'$(45.0, -31.0)$', r'$(135.0, 0.0)$', r'$(180.0, 0.0)$', r'$(230.0, 0.0)$']
# COLORs = ['orange', 'purple', 'royalblue', 'darkgreen', 'gold', 'chocolate', 'brown']
COLORs = cm.rainbow(np.linspace(0,1,7))

xvals = [xval, xval, xval, xval, xval, xval, xval]
yvals = [yval_old, yval_new25, yval_new10, yval_new45, yval_new135, yval_new180, yval_new230]

LINEs = ['-', '--', '--', '--', '--', '--', '--']
LINEWs = [2, 2, 2, 2, 2, 2, 2]
POINTs = ['', '', '', '', '', '', '']

plotting.LinePlotFunc(outpath,
                xvals, yvals,
                COLORs, LABELs=LABELs, LINEs=LINEs, LINEWs=LINEWs, POINTs=POINTs, POINTSs=None, fillstyles=None,
                XRANGE=XRANGE, YRANGE=YRANGE,
                XLABEL=XLABEL, YLABEL=YLABEL, TITLE=None,
                xtick_min_label=True, xtick_spe=None, ytick_min_label=True, ytick_spe=None,
                vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                xlog=False, invertX=False, ylog=False, invertY=False, loc_legend='upper left',
                font_size=16, usetex=True)