# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-08-26 10:30:10
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-12 15:24:04

### m difference by PSF modelling

import numpy as np
import pandas as pd

from matplotlib.pyplot import cm
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams["text.usetex"] =True

# +++ general settings for plot
# mpl.use('Agg')
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
plt.rc('font', size=16)

# +++ I/O
outpath = "./plots/dm_PSFmodelling.pdf"
# outpath = "show"
inpath_list = ['../sensitivity_test/PSFmodelling/outputs/dm_ZB_test_fiducial_rewei_41.csv']

# +++ general info
m_error_budget_list = [0.02, 0.01, 0.01, 0.01, 0.01, 0.02]

# +++ plot related
TITLE = None
XLIM = [-0.025, 0.025]
COLORs = ['k']
SYMBOLs = ['o']
MSs = [7]
ELWs = [1]

YTICK = [1, 2, 3, 4, 5, 6]
YTICKLABELS = [r'$0.1< z_{\rm B} \leq 0.3$', 
                r'$0.3< z_{\rm B} \leq 0.5$', 
                r'$0.5< z_{\rm B} \leq 0.7$', 
                r'$0.7< z_{\rm B} \leq 0.9$', 
                r'$0.9< z_{\rm B} \leq 1.2$',
                r'$1.2< z_{\rm B} \leq 2.0$']
XLABEL = r"$\Delta m_{\rm PSF}$ (model PSF - input PSF)"
YLIM = [0.5, 6.5]

# +++ plot
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.1, right=0.98, top=0.88, wspace=0, hspace=0)

for i, inpath in enumerate(inpath_list):

    data = pd.read_csv(inpath)
    print('for', inpath)

    # first point is for the whole sample
    data = data[1:]

    # x position as m
    x_val = (data['m1'].values + data['m2'].values)/2.
    # error
    x_err = (data['m1_err'].values + data['m2_err'].values)/2.
    print('m', x_val)
    print('err', x_err)

    # y position as redshift bins
    y_cen = np.arange(len(data)) + 1
    y_val = y_cen 

    CR = COLORs[i]
    MK = SYMBOLs[i]
    MS = MSs[i]
    ELW = ELWs[i]

    # error budget shadow
    if (m_error_budget_list is not None) and (i == 0):
        for kk in range(len(x_val)):
            x0 = 0
            y0 = y_cen[kk]

            x = np.array([x0 - m_error_budget_list[kk], x0 + m_error_budget_list[kk]])
            y1 = y0 - 0.5
            y2 = y0 + 0.5
            plt.fill_between(x, y1, y2, edgecolor='gray', facecolor='white', alpha=0.3, hatch='/')

    plt.errorbar(x_val, y_val, xerr=x_err,
                    color=CR, marker=MK, markersize=MS, elinewidth=ELW, ls='none')

for ibin in range(6):
    plt.axhline(y=1.5+ibin, color='black', ls='-', lw=1)

plt.axvline(x=0, color='gray', ls = '--', lw=1)

plt.yticks(ticks=YTICK, labels=YTICKLABELS)
plt.tick_params(axis='y', length=0, width=0)

plt.xlim(XLIM[0], XLIM[1])
plt.ylim(YLIM[0], YLIM[1])
plt.xlabel(XLABEL)

# invert y-axis
plt.gca().invert_yaxis()

plt.title(TITLE)

plt.tight_layout()

if outpath == 'show':
    plt.show()
else:
    plt.savefig(outpath, dpi=300)
    print('plot saved in', outpath)
plt.close()


# for ../sensitivity_test/PSFmodelling/outputs/dm_ZB_test_fiducial_rewei_41.csv
# m [0.00205 0.00405 0.0039  0.00325 0.005   0.00735]
# err [0.0013  0.00115 0.00125 0.00105 0.0012  0.00175]
# plot saved in ./plots/dm_PSFmodelling.pdf
