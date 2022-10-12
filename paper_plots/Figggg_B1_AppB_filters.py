# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-08-30 15:01:39
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-29 19:58:09

### compare the filters

import os
import pandas as pd 
import numpy as np 

import math

import numpy as np
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator, LogLocator, NullFormatter, NullLocator


# ++++++++++++++ general info

# where to save
outpath = 'show'
outpath = './plots/filtersComp.pdf'

# where to find the filters
indir = '/disks/shear10/ssli/ImSim/code/bpz-1.99.3_expanded/FILTER/'

# bands 
bands = [r'$u$', r'$g$', r'$r$', r'$i$', r'$Z$']
filters_kids = ['KiDSVIKING_u_positive.res', 
                'KiDSVIKING_g_positive.res', 
                'KiDSVIKING_r_positive.res', 
                'KiDSVIKING_i_positive.res', 
                'KiDSVIKING_Z2_positive.res']
filters_shark = ['shark_u_SDSS.res',
                'shark_g_SDSS.res',
                'shark_r_SDSS.res',
                'shark_i_SDSS.res',
                'shark_z_SDSS.res']

# plot related
font_size = 16
usetex = True
LW = 1.5
LINES = ['-', '--']
COLORs = ['r', 'green']
LABELs = ['KiDS/VIKING', 'SDSS'] 
XLABEL = 'Wavelength [Angstrom]'
YLABEL = 'Normalised transmission'
XRANGE = [2500, 11500]
YRANGE = [0, 1.2]
xtick_spe = None
ytick_spe = [[0, 0.5, 1.0], ['0.0', '0.5', '1.0']]

# >>>>>>>> get data
N_plots = len(bands)
xvals_k = []
yvals_k = []
xvals_s = []
yvals_s = []
for i_band, band in enumerate(bands):

    # kids
    inpath = os.path.join(indir, filters_kids[i_band])
    cata_kids = np.loadtxt(inpath)
    xvals_k.append(cata_kids[:, 0])
    yvals_k.append(cata_kids[:, 1]/np.max(cata_kids[:, 1]))

    # shark
    inpath = os.path.join(indir, filters_shark[i_band])
    cata_shark = np.loadtxt(inpath)
    xvals_s.append(cata_shark[:, 0])
    yvals_s.append(cata_shark[:, 1]/np.max(cata_shark[:, 1]))


# >>>>>>>> plot
# font size
plt.rc('font', size=font_size)
# tex
plt.rcParams["text.usetex"] = usetex

if outpath != 'show':
    backend_orig = plt.get_backend()
    plt.switch_backend("agg")

fig, ax = plt.subplots()
for i, xvl_k in enumerate(xvals_k):

    # kids
    yvl_k = yvals_k[i]
    ix_kids = np.argmax(yvl_k)
    xvl_max = xvl_k[ix_kids]
    if i == 0:
        plt.plot(xvl_k, yvl_k, color=COLORs[0], label=LABELs[0], 
                    linestyle=LINES[0], linewidth=LW, marker='')
    else:
        plt.plot(xvl_k, yvl_k, color=COLORs[0],
                    linestyle=LINES[0], linewidth=LW, marker='')
    del xvl_k, yvl_k

    # shark
    xvl_s = xvals_s[i]
    yvl_s = yvals_s[i]
    if i == 0:
        plt.plot(xvl_s, yvl_s, color=COLORs[1], label=LABELs[1], 
                    linestyle=LINES[1], linewidth=LW, marker='')
    else:
        plt.plot(xvl_s, yvl_s, color=COLORs[1], 
                    linestyle=LINES[1], linewidth=LW, marker='')

    del xvl_s, yvl_s

    # set band label
    plt.text(x=xvl_max, y=1.1, s=bands[i])

plt.xlim(XRANGE[0], XRANGE[1])
plt.ylim(YRANGE[0], YRANGE[1])
plt.legend(frameon=True, loc='center', ncol=2, bbox_to_anchor=(0.5, 1.08))

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

if xtick_spe is not None:
    plt.xticks(xtick_spe[0], xtick_spe[1])
if ytick_spe is not None:
    plt.yticks(ytick_spe[0], ytick_spe[1])

plt.xlabel(XLABEL)
plt.ylabel(YLABEL)

if outpath=='show':
    plt.show()
    plt.close()
else:
    plt.savefig(outpath, dpi=300)
    plt.close()
    plt.switch_backend(backend_orig)
    print("Line plot saved as", outpath)
