# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-01-05 15:05:29
# @Last Modified by:   lshuns
# @Last Modified time: 2022-01-05 16:02:26

### footprint of the selected tiles

import os
import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

# >>>>>>>>>>>>>>>>>> I/O

# all
noise_all = pd.read_csv('../noise_info/kids_dr4_noise_ugriZYJHKs.csv')

# selected
infile = '../noise_info/skills_fiducial/noise_selec_combined.csv'
noise_selec = pd.read_csv(infile)

# save
outpath = './plots/noise_footprint.pdf'

# >>>>>>>>>>>>>>>>>> setting

noise_info_list = [noise_all, noise_selec]
color_list = ['gray', 'blue']

title_list = ['KiDS-North', 'KiDS-South']
xrange_list = [[125, 245], [-35, 55]]
yrange_list = [[-8, 8], [-38, -22]]

# >>>>>>>>>>>>>>>>>> plot
# font size
plt.rc('font', size=16)
# tex
plt.rcParams["text.usetex"] = True

if outpath != 'show':
    backend_orig = plt.get_backend()
    plt.switch_backend("agg")

fig, axs = plt.subplots(2, 1)

fig.subplots_adjust(hspace=0.4)

# loop through sky 
for i_plot in range(2):

    ax = axs[i_plot]
    title = title_list[i_plot]
    XRANGE = xrange_list[i_plot]
    YRANGE = yrange_list[i_plot]

    ## loop through data
    for i_noise, noise_info in enumerate(noise_info_list):

        colour = color_list[i_noise]

        ### location info
        loc_info = noise_info['label'].str.split(r"_", expand=True)
        loc_info.columns = ['ra', 'dec']
        loc_info = loc_info.astype({'ra': 'float', 'dec': 'float'})
        loc_info.loc[loc_info.ra>300., 'ra'] = loc_info.loc[loc_info.ra>300., 'ra']-360.

        ### create patches
        x_low = loc_info['ra'].values-0.5
        y_low = loc_info['dec'].values-0.5
        xy_list = np.empty((x_low.size, 2), dtype=x_low.dtype)
        xy_list[:,0] = x_low
        xy_list[:,1] = y_low
        patches = []
        for xy in xy_list:
            rect = Rectangle(xy, 1, 1)
            patches.append(rect)
        pc_loc = PatchCollection(patches, facecolor=colour, alpha=1,
                                 edgecolor='None')

        ### plot
        ax.add_collection(pc_loc)

        ### Star points
        ra_stars = [135, 180, 230, 10, 45, -25]
        dec_stars = [0, 0, 0, -31, -31, -31]
        for i_ra, ra_star in enumerate(ra_stars):
            dec_star = dec_stars[i_ra]
            ax.plot(ra_stars, dec_stars, ls='', color='black', ms=10, marker='*')

        ### title
        ax.set_title(title, fontsize=12)

        ### range
        ax.set_xlim(XRANGE[0], XRANGE[1])
        ax.set_ylim(YRANGE[0], YRANGE[1])
        ax.invert_xaxis()

### labels
fig.text(0.5, 0.02, 'R.A. [deg]', ha='center')
fig.text(0.02, 0.5, 'Dec [deg]', va='center', rotation='vertical')

# save
if outpath == 'show':
    plt.show()
    plt.close()
else:
    plt.savefig(outpath, dpi=300)
    plt.close()
    plt.switch_backend(backend_orig)
    print("plot saved as", outpath)