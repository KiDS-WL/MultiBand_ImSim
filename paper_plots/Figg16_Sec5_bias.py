# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-08-17 15:53:26
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-12 15:20:57

### the final m-bias results

import numpy as np
import pandas as pd
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

# # +++ input
inpath_list = ['../biasEstimation/outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut.csv', 
                '../correction_varShear_PSFmodelling/outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_varCorr_dz0p1_2D_PSFmodelling.csv']

# m_error_budget_list = None
m_error_budget_list = [0.02, 0.01, 0.01, 0.01, 0.01, 0.02]

# output directory
outpath = "show"
outpath = './plots/m_ZB.pdf'

# plot related
LABELS = ['constant shear', 'variable shear']
COLORS = ['gray', 'red']
SYMBOLS = ['o', 'D']
# marker size
MSs = [3, 5]
# linewidth of the errorbar lines
ELWs = [1, 1]

YTICK = [1, 2, 3, 4, 5, 6]
YTICKLABELS = [r'$0.1< z_{\rm B} \leq 0.3$', 
                r'$0.3< z_{\rm B} \leq 0.5$', 
                r'$0.5< z_{\rm B} \leq 0.7$', 
                r'$0.7< z_{\rm B} \leq 0.9$', 
                r'$0.9< z_{\rm B} \leq 1.2$',
                r'$1.2< z_{\rm B} \leq 2.0$']

XLABEL = r"$m$"
XLIM = [-0.05, 0.1]
YLIM = [0.5, 6.5]

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
    if i == 0 :
        y_val = y_cen - 0.05
    else:
        y_val = y_cen + 0.05

    CR = COLORS[i]
    MK = SYMBOLS[i]
    MS = MSs[i]
    ELW = ELWs[i]
    LABEL = LABELS[i]

    # error budget shadow
    if (m_error_budget_list is not None) and (i == 1):
        for kk in range(len(x_val)):
            x0 = x_val[kk]
            y0 = y_cen[kk]

            x = np.array([x0 - m_error_budget_list[kk], x0 + m_error_budget_list[kk]])
            y1 = y0 - 0.5
            y2 = y0 + 0.5
            plt.fill_between(x, y1, y2, edgecolor='gray', facecolor='white', alpha=0.3, hatch='/')

    plt.errorbar(x_val, y_val, xerr=x_err,
                    color=CR, marker=MK, markersize=MS, elinewidth=ELW, ls='none', label=LABEL)

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

plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
             fancybox=True, shadow=True, ncol=3)

if outpath == 'show':
    plt.show()
else:
    plt.savefig(outpath, dpi=300)
    print('plot saved in', outpath)
plt.close()

# # # # # ########### out info
# for ../biasEstimation/outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut.csv
# m [-0.0124  -0.0205  -0.0064   0.02185  0.03335  0.0636 ]
# err [0.0058 0.0035 0.0037 0.0042 0.0049 0.0069]
# for ../correction_varShear_PSFmodelling/outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_varCorr_dz0p1_2D_PSFmodelling.csv
# m [-0.01262962 -0.01814822 -0.0079384   0.01890082  0.0336984   0.07181115]
# err [0.01663706 0.00707949 0.00665579 0.00554103 0.00585646 0.00809717]
# plot saved in ./plots/m_ZB.pdf
