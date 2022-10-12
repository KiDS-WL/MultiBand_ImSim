# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-07-07 09:37:00
# @Last Modified by:   lshuns
# @Last Modified time: 2022-08-25 18:01:33

### cosmic shear vs. redshift

import plotting

from matplotlib.pyplot import cm

import os
import glob

import numpy as np
import pandas as pd
import astropy.units as u

from scipy.stats import truncnorm

from astropy.io import fits
from astropy import constants
from astropy.cosmology import LambdaCDM

# >>>>>>>>>>>>>>>>>>> I/O

indir = '/disks/shear10/ssli/ImSim/input/MICE2_cata/'

# redshift bins
inpath = os.path.join(indir, 'zB_centres_measured.npy')
zB_centres_measured = np.load(inpath, allow_pickle=True)
print('zB measured', zB_centres_measured)
N_measured = len(zB_centres_measured)
inpath = os.path.join(indir, 'zB_centres_learned.npy')
zB_centres_learned = np.load(inpath, allow_pickle=True)
print('zB learned', zB_centres_learned)
N_learned = len(zB_centres_learned)
# gI
inpath = os.path.join(indir, 'gI_sigma_measured.npy')
gI_sigma_measured = np.load(inpath, allow_pickle=True)
print('gI measured', gI_sigma_measured)
inpath = os.path.join(indir, 'gI_sigma_learned.npy')
gI_sigma_learned = np.load(inpath, allow_pickle=True)
print('gI learned', gI_sigma_learned)


# redshift
inpath = os.path.join(indir, 'zbins_centres_measured.npy')
zbins_centres_measured = np.load(inpath)
print('zbins measured', zbins_centres_measured)
inpath = os.path.join(indir, 'zbins_centres_learned.npy')
zbins_centres_learned = np.load(inpath)
print('zbins learned', zbins_centres_learned)
# gF
inpath = os.path.join(indir, 'gF_mean_measured.npy')
gF_mean_measured = np.load(inpath)
print('gF measured', gF_mean_measured)
inpath = os.path.join(indir, 'gF_mean_learned.npy')
gF_mean_learned = np.load(inpath)
print('gF learned', gF_mean_learned)

# >>>>>>>>>>>>>>>>>> plot
outpath = './plots/Zshear_CS.pdf'
# outpath = 'show'

xvals = [zbins_centres_measured, zbins_centres_learned] + zB_centres_measured[[4, 8]].tolist() + zB_centres_learned[[4, 8, 14, 18]].tolist()

yvals = [gF_mean_measured, gF_mean_learned] + gI_sigma_measured[[4, 8]].tolist() + gI_sigma_learned[[4, 8, 14, 18]].tolist()

yerrs = None

N_measured = 2
N_learned = 4
COLORs_N = np.linspace(0,1,N_measured + N_learned)
COLORs = ['k', 'k'] + np.vstack([cm.rainbow(COLORs_N[:N_measured]), cm.rainbow(COLORs_N)]).tolist()

LABELs = [r'$|\gamma_{\rm F}|$', None,
            r'$\sigma_{\gamma_{\rm I}}~(z_{\rm F}=0.65)$',
            r'$\sigma_{\gamma_{\rm I}}~(z_{\rm F}=1.05)$',
            None, None,
            r'$\sigma_{\gamma_{\rm I}}~(z_{\rm F}=1.75)$',
            r'$\sigma_{\gamma_{\rm I}}~(z_{\rm F}=2.15)$']

LINEs = ['-', '-'] + ['--'] * N_measured + ['--'] * N_learned
LINEWs = [1, 2] + [1] * N_measured + [2] * N_learned

POINTs = ['o', ''] + ['o'] * N_measured + [''] * N_learned
POINTSs = [5] * len(xvals)

XLABEL = 'Redshift'
YLABEL = 'Mean amplitude or dispersion'

TITLE = None

XRANGE = None
YRANGE = None

plotting.ErrorPlotFunc(outpath,
                xvals, yvals, yerrs,
                COLORs, LABELs=LABELs, LINEs=LINEs, LINEWs=LINEWs, POINTs=POINTs, POINTSs=POINTSs, ERRORSIZEs=None,
                XRANGE=XRANGE, YRANGE=YRANGE,
                XLABEL=XLABEL, YLABEL=YLABEL, TITLE=TITLE,
                xtick_min_label=True, xtick_spe=None, ytick_min_label=True, ytick_spe=None,
                vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                xlog=False, invertX=False, ylog=False, invertY=False, loc_legend='best',
                font_size=16, usetex=True)
