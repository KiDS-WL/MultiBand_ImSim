# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-07-10 16:13:25
# @Last Modified by:   lshuns
# @Last Modified time: 2022-08-26 10:30:46

### m difference by splitting fiducial results

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

# +++ input

############### noise level
outpath = "./plots/dm_noise.pdf"
# outpath = "show"
label_title = "Noise level :"
inpath_list = ['../sensitivity_test/variance_within_skills/outputs/dm_noise_part_low.csv',
            '../sensitivity_test/variance_within_skills/outputs/dm_noise_part_median.csv',
            '../sensitivity_test/variance_within_skills/outputs/dm_noise_part_high.csv']

# ############### PSFe
# outpath = "./plots/dm_PSFe.pdf"
# # outpath = 'show'
# label_title = r"PSF $\epsilon$ :"
# inpath_list = ['../sensitivity_test/variance_within_skills/outputs/dm_PSFe_part_low.csv',
#             '../sensitivity_test/variance_within_skills/outputs/dm_PSFe_part_median.csv',
#             '../sensitivity_test/variance_within_skills/outputs/dm_PSFe_part_high.csv']

# ############### PSFsize
# outpath = "./plots/dm_PSFsize.pdf"
# # outpath = "show"
# label_title = "PSF size :"
# inpath_list = ['../sensitivity_test/variance_within_skills/outputs/dm_PSFsize_part_low.csv',
#             '../sensitivity_test/variance_within_skills/outputs/dm_PSFsize_part_median.csv',
#             '../sensitivity_test/variance_within_skills/outputs/dm_PSFsize_part_high.csv']

# ############### stellar density
# outpath = "./plots/dm_stellar.pdf"
# # outpath = "show"
# label_title = "Stellar density :"
# inpath_list = ['../sensitivity_test/variance_within_skills/outputs/dm_stellar_part_low.csv',
#             '../sensitivity_test/variance_within_skills/outputs/dm_stellar_part_median.csv',
#             '../sensitivity_test/variance_within_skills/outputs/dm_stellar_part_high.csv']

# +++ general info
m_error_budget_list = [0.02, 0.01, 0.01, 0.01, 0.01, 0.02]

# +++ plot related
XLIM = [-0.025, 0.025]
LABELs = ['low', 'median', 'high']
COLORs = ['orange', 'm', 'cyan']
SYMBOLs = ['<', 'o', '>']
MSs = [7, 7, 7]
ELWs = [1, 1, 1]
Nsets = len(LABELs)

YTICK = [1, 2, 3, 4, 5, 6]
YTICKLABELS = [r'$0.1< z_{\rm B} \leq 0.3$', 
                r'$0.3< z_{\rm B} \leq 0.5$', 
                r'$0.5< z_{\rm B} \leq 0.7$', 
                r'$0.7< z_{\rm B} \leq 0.9$', 
                r'$0.9< z_{\rm B} \leq 1.2$',
                r'$1.2< z_{\rm B} \leq 2.0$']
XLABEL = r"$\Delta m$"
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
    y_val = y_cen + 0.1 * (i - 0.5*Nsets)

    CR = COLORs[i]
    MK = SYMBOLs[i]
    MS = MSs[i]
    ELW = ELWs[i]
    LABEL = LABELs[i]

    # error budget shadow
    if (m_error_budget_list is not None) and (i == 1):
        for kk in range(len(x_val)):
            x0 = 0
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

# legend
h, l = ax.get_legend_handles_labels() # Extracting handles and labels
ph = [plt.plot([],marker="", ls="")[0]] # Canvas
handles = ph + h
labels = [label_title] + l  # Merging labels
leg = plt.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.15),
                handletextpad=0.2, columnspacing=0.6, 
             ncol=Nsets+1)
for vpack in leg._legend_handle_box.get_children()[:1]:
    for hpack in vpack.get_children():
        hpack.get_children()[0].set_width(0)

if outpath == 'show':
    plt.show()
else:
    plt.savefig(outpath, dpi=300)
    print('plot saved in', outpath)
plt.close()

