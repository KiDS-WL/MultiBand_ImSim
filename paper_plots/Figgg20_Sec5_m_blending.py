# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-09-11 13:34:23
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-12 15:22:58

### m as a function of the neighbour distance

import numpy as np 
import pandas as pd 

import matplotlib.pyplot as plt 
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

from matplotlib.ticker import AutoMinorLocator, LogLocator, NullFormatter, NullLocator

# +++++++++++++ I/O

outpath = 'show'
outpath = './plots/m_blending.pdf'

inpath_list = ['../correction_varShear_PSFmodelling/outputs/m_forBlending_neighbour_dist_var.csv',
                    '../correction_varShear_PSFmodelling/outputs/m_forBlending_neighbour_dist_const.csv']
LABELs = ['varShear', 'constShear']
COLORs = ['magenta', 'darkgreen']
SYMBOLs = ['v', '^']
MSs = [6, 6]
ELWs = [1, 1]
LSs = ['--', ':']
LWs = [1, 1]

inpath_const_whole = '../correction_varShear_PSFmodelling/outputs/m_forBlending_neighbour_dist_const_whole.csv'

inpath_dm = '../correction_varShear_PSFmodelling/outputs/dm_forBlending_neighbour_dist.csv'
CRr = 'k'
MKr = 'o'
MSr = 6
ELWr = 1
LSr = 'dashdot'
LWr = 1

col_binning = 'nnb_dist'
XLABEL = 'Input distance to the nearest neighbour [arcsec]'

blending_cut = 4

# +++++++++++++ plot

plt.rc('font', size=16)
plt.rcParams["text.usetex"] =True

#### the upper plot for m
rect = [0.13, 0.37, 0.8, 0.6] ### rect = [left, bottom, width, height]
ax1 = plt.axes(rect)
for i, inpath in enumerate(inpath_list):

    data = pd.read_csv(inpath)
    print('for', inpath)

    # first point is for the whole sample
    data = data[1:]

    # x position from mean
    x_val = data[f'{col_binning}_mean'].values
    print(col_binning, x_val)

    # y position from m
    y_val = (data['m1'].values + data['m2'].values)/2.
    # error
    y_err = (data['m1_err'].values + data['m2_err'].values)/2.
    print('m', y_val)
    print('err', y_err)

    # >>>> get extra points from constant whole
    if LABELs[i] == 'constShear':
        data = pd.read_csv(inpath_const_whole)
        # first point is for the whole sample
        data = data[1:]

        # x position from mean
        x_val_tmp = data[f'{col_binning}_mean'].values[len(x_val):]
        print('extra points', col_binning, x_val_tmp)

        # y position from m
        y_val_tmp = ((data['m1'].values + data['m2'].values)/2.)[len(x_val):]
        # error
        y_err_tmp = ((data['m1_err'].values + data['m2_err'].values)/2.)[len(x_val):]
        print('extra points', 'm', y_val_tmp)
        print('extra points', 'err', y_err_tmp)

        # combine
        x_val = np.concatenate([x_val, x_val_tmp])
        y_val = np.concatenate([y_val, y_val_tmp])
        y_err = np.concatenate([y_err, y_err_tmp])

    LAB = LABELs[i]
    CR = COLORs[i]
    MK = SYMBOLs[i]
    MS = MSs[i]
    ELW = ELWs[i]
    LS = LSs[i]
    LW = LWs[i]
    ax1.errorbar(x_val, y_val, yerr=y_err, label=LAB,
                    color=CR, marker=MK, markersize=MS, elinewidth=ELW, ls=LS, linewidth=LW)
ax1.axhline(y=0, color='gray', ls = (0, (3, 5, 1, 5, 1, 5)), lw=1)

ax1.axvline(x=blending_cut, color='red', ls = '--', lw=1)

ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())

ax1.legend(loc='lower right', frameon=True, framealpha=1)

ax1.set_ylabel(r'$m$')
# ax1.set_xlim(XLIM[0], XLIM[1])
# ax1.set_ylim(YLIM[0], YLIM[1])

# ax1.set_xticklabels(" ")

#### the lower panel for residuals
rect = [0.13, 0.12, 0.8, 0.2]
ax2 = plt.axes(rect)
ax2.sharex(ax1)

data = pd.read_csv(inpath_dm)
print('for', inpath_dm)
# first point is for the whole sample
data = data[1:]
# x position from mean
x_val = data[f'{col_binning}_mean'].values
print(col_binning, x_val)
# y position from m
y_val = (data['m1'].values + data['m2'].values)/2.
# error
y_err = (data['m1_err'].values + data['m2_err'].values)/2.
print('dm', y_val)
print('dm err', y_err)
ax2.errorbar(x_val, y_val, yerr=y_err,
                color=CRr, marker=MKr, markersize=MSr, elinewidth=ELWr, ls=LSr, linewidth=LWr)
ax2.axhline(y=0, color='gray', ls = (0, (3, 5, 1, 5, 1, 5)), lw=1)

ax2.set_xlabel(XLABEL)
ax2.set_ylabel(r'$\Delta m$')

ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())

ax2.axvline(x=blending_cut, color='red', ls = '--', lw=1)

ax1.tick_params('x', labelbottom=False)

if outpath == 'show':
    plt.show()
else:
    plt.savefig(outpath, dpi=300)
    print('plot saved in', outpath)
plt.close()

# for ../correction_varShear_PSFmodelling/outputs/m_forBlending_neighbour_dist_var.csv
# nnb_dist [0.82334992 1.67265598 2.21997374 2.59934614 2.93364848 3.24863022
#  3.5527029  3.90818523]
# m [-0.3119  -0.08275 -0.02325 -0.00295  0.0106   0.0131   0.03515  0.02925]
# err [0.01835 0.0173  0.0144  0.01335 0.01275 0.0124  0.01225 0.0121 ]
# for ../correction_varShear_PSFmodelling/outputs/m_forBlending_neighbour_dist_const.csv
# nnb_dist [0.82296831 1.67248382 2.2198454  2.59925156 2.93358792 3.24862149
#  3.55271334 3.91143631]
# m [-0.26305 -0.0849  -0.0157  -0.001    0.01855  0.0227   0.0291   0.03585]
# err [0.0073  0.00665 0.0057  0.0054  0.0052  0.0051  0.005   0.005  ]
# extra points nnb_dist [4.2479106  4.74720296]
# extra points m [0.03355 0.0362 ]
# extra points err [0.004 0.004]
# for ../correction_varShear_PSFmodelling/outputs/dm_forBlending_neighbour_dist.csv
# nnb_dist [0.82318394 1.67245006 2.21986366 2.59923835 2.93360199 3.24861963
#  3.5527218  3.9082138 ]
# dm [-0.04885  0.00205 -0.0076  -0.00185 -0.00805 -0.00955  0.00605 -0.0066 ]
# dm err [0.01975 0.0185  0.0155  0.01445 0.01375 0.0134  0.01325 0.0131 ]
# plot saved in ./plots/m_blending.pdf
# Elapsed:0:07.80,User=7.060,System=22.964,CPU=384.8%.