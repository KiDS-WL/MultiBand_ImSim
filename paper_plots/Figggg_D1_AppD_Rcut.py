# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-09-26 12:45:36
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-26 16:35:14

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
outpath = './plots/Rcut.pdf'

inpath = '../biasEstimation/outputs/alpha_kids_photo_LF_321_shear_noSG_noWeiCut_noRcut_817tiles.csv'
cata = pd.read_csv(inpath)

Rcut = 0.9

# XLIM = [0, 1]

LABELs = [r'$\epsilon_1$', r'$\epsilon_2$']
COLORs = ['darkred', 'darkorange']
SYMBOLs = ['v', '^']
LSs = ['--', ':']
MS = 6
ELW = 1
LW = 1

XLABEL = r'$\mathcal{R}$'

# +++++++++++++ plot

plt.rc('font', size=16)
plt.rcParams["text.usetex"] =True

#### the upper plot for alpha
rect = [0.13, 0.5, 0.8, 0.4] ### rect = [left, bottom, width, height]
ax1 = plt.axes(rect)

# x position from mean
x_val = cata.loc[1:, 'R_mean'].values
print('x', x_val)

# e1 results
i_col = 0
y_val = cata.loc[1:, 'alpha1'].values
y_err = cata.loc[1:, 'alpha1_err'].values
ax1.errorbar(x_val, y_val, yerr=y_err, label=LABELs[i_col],
                color=COLORs[i_col], marker=SYMBOLs[i_col], 
                ls=LSs[i_col],
                markersize=MS, elinewidth=ELW, linewidth=LW)

# e2 results
i_col = 1
y_val = cata.loc[1:, 'alpha2'].values
y_err = cata.loc[1:, 'alpha2_err'].values
ax1.errorbar(x_val, y_val, yerr=y_err, label=LABELs[i_col],
                color=COLORs[i_col], marker=SYMBOLs[i_col], 
                ls=LSs[i_col],
                markersize=MS, elinewidth=ELW, linewidth=LW)

ax1.axhline(y=0, color='gray', ls = (0, (3, 5, 1, 5, 1, 5)), lw=1)

ax1.axvline(x=Rcut, color='red', ls = '--', lw=1)

ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())

ax1.legend(loc='upper left', frameon=True, framealpha=1)

ax1.set_ylabel(r'$\alpha$')
# ax1.set_xlim(XLIM[0], XLIM[1])
# ax1.set_ylim(YLIM[0], YLIM[1])

# ax1.set_xticklabels(" ")

#### the lower panel for residuals
rect = [0.13, 0.12, 0.8, 0.36]
ax2 = plt.axes(rect)
ax2.sharex(ax1)

# y from cumulative sum
y_val = np.cumsum(cata.loc[1:, 'Nwei'].values) / cata.loc[0, 'Nwei']
print('cumsum', y_val)
ax2.errorbar(x_val, y_val, yerr=None,
                color='k', marker='', ls='-', linewidth=1)
ax2.axhline(y=1, color='gray', ls = (0, (3, 5, 1, 5, 1, 5)), lw=1)

ax2.set_xlabel(XLABEL)
ax2.set_ylabel('Cumulative distribution')

ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())

ax2.axvline(x=Rcut, color='red', ls = '--', lw=1)

ax1.tick_params('x', labelbottom=False)

# reverse x
# plt.gca().invert_xaxis()

if outpath == 'show':
    plt.show()
else:
    plt.savefig(outpath, dpi=300)
    print('plot saved in', outpath)
plt.close()

# for ../correction_varShear_PSFmodelling/outputs/m_forBlending_neighbour_dist_var.csv
# nnb_dist [0.82339372 1.67269599 2.22001468 2.59939236 2.93370191 3.24868718
#  3.55275972 3.911338  ]
# m [-0.3105  -0.07885 -0.02495 -0.00395  0.0079   0.01575  0.03305  0.0294 ]
# err [0.0255  0.0242  0.0203  0.01885 0.01795 0.0175  0.01725 0.0171 ]
# for ../correction_varShear_PSFmodelling/outputs/m_forBlending_neighbour_dist_const.csv
# nnb_dist [0.82296831 1.67248382 2.2198454  2.59925156 2.93358792 3.24862149
#  3.55271334 3.91143631]
# m [-0.26305 -0.0849  -0.0157  -0.001    0.01855  0.0227   0.0291   0.03585]
# err [0.0073  0.00665 0.0057  0.0054  0.0052  0.0051  0.005   0.005  ]
# extra points nnb_dist [4.2479106  4.74720296]
# extra points m [0.03355 0.0362 ]
# extra points err [0.004 0.004]
# for ../correction_varShear_PSFmodelling/outputs/dm_forBlending_neighbour_dist.csv
# nnb_dist [0.82319264 1.67243876 2.21987034 2.59923926 2.93360201 3.24862077
#  3.55272399 3.91132941]
# dm [-0.04755  0.00615 -0.0091  -0.00305 -0.0107  -0.00695  0.0039  -0.00645]
# dm err [0.02655 0.0251  0.0211  0.01965 0.0187  0.0182  0.01795 0.0178 ]
# plot saved in ./plots/m_blending.pdf
# Elapsed:0:04.51,User=6.976,System=20.768,CPU=614.8%.
