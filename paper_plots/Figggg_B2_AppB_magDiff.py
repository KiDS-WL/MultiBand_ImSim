# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-08-30 09:11:01
# @Last Modified by:   lshuns
# @Last Modified time: 2022-08-30 17:28:26

### compare magnitude before and after modification

import os 
import glob

import pandas as pd 
import numpy as np 

import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator, LogLocator, NullFormatter, NullLocator

# >>>>>>>> I/O

# where to save
outpath = 'show'
outpath = './plots/magDiff.pdf'

# where to find the original photometry catalogues
maindir = '/disks/shear16/ssli/ImSim/output/skills_v07D7'
parts = ['part0', 'part1', 'part2', 'part3', 'part4', 'part5']
shear_tag = 'm283m283'
rot_tag = 'rot0'

# which bands are modified
tag_bands = ['MAG_GAAP_u', 'MAG_GAAP_g', 'MAG_GAAP_r', 'MAG_GAAP_i', 'MAG_GAAP_Z', 'MAG_AUTO']
label_bands = ['MAG\_GAAP\_ u', 'MAG\_GAAP\_ g', 'MAG\_GAAP\_ r', 'MAG\_GAAP\_ i', 'MAG\_GAAP\_ Z', 'MAG\_AUTO']

# plot related
FIGSIZE = (10, 6)
N_rows = 2
N_cols = 3
font_size = 20
usetex = True
nbins = 60
COLOR_MAP='Blues'
line_colour = 'b'
XLABEL = r'Measured magnitude (in SDSS filters)'
YLABEL = r'Mag difference (KiDS - SDSS)'
XRANGE = [19, 25]  
YRANGE = [-0.09, 0.09] 

# >>>>>>>> workhorse

# loop through all files to get values
cata_list = []
for part in parts:
    print('>>> in', part)
    subdir_ori = os.path.join(maindir, part, shear_tag, 'catalogues')
    subdir_corr = os.path.join(maindir, part, shear_tag, 'catalogues')

    # find all files
    file_list = glob.glob(os.path.join(subdir_ori, f'*_{rot_tag}_combined.feather'))
    print('number of files found', len(file_list))
    ## loop through all files
    for infile in file_list:

        # load the cata
        #### original
        cata_ori = pd.read_feather(infile)
        ## only used good measurements
        flag_9 = np.zeros(len(cata_ori)) 
        for band in ['u', 'g', 'r', 'i', 'Z', 'Y', 'J', 'H', 'Ks']:
            flag_9 += cata_ori[f'FLAG_GAAP_{band}'].values
        mask_gaap = (flag_9==0)
        del flag_9
        ## used columns
        cata_ori = cata_ori[tag_bands]
        ## rename 
        col_names = {tag_band: tag_band+'_ori' for tag_band in tag_bands}
        cata_ori.rename(columns=col_names, inplace=True)
        #### modified
        cata_new = pd.read_feather(infile.replace('_combined.feather', '_combined_kids_filters.feather'))
        ## used columns
        cata_new = cata_new[tag_bands]
        ## rename 
        col_names = {tag_band: tag_band+'_new' for tag_band in tag_bands}
        cata_new.rename(columns=col_names, inplace=True)
        del col_names

        # combine
        cata = pd.concat([cata_ori, cata_new], axis=1)
        del cata_ori, cata_new

        # get good measurements
        cata = cata[mask_gaap]
        del mask_gaap

        # save 
        cata_list.append(cata)
        del cata

# combine
cata_list = pd.concat(cata_list, ignore_index=True)
print('total number of objects', len(cata_list))

# >>> plot
# font size
plt.rc('font', size=font_size)
# tex
plt.rcParams["text.usetex"] = usetex

if outpath != 'show':
    backend_orig = plt.get_backend()
    plt.switch_backend("agg")

# how many sub-plots
fig, axs = plt.subplots(N_rows, N_cols, figsize=FIGSIZE, sharex=True, sharey=True)
fig.subplots_adjust(hspace=0.1)
fig.subplots_adjust(wspace=0.1)

# some general relation
SRANGE = [[XRANGE[0], XRANGE[1]], [YRANGE[0], YRANGE[1]]]

# loop over sub-plots
i_plot = 0
for i_row in range(N_rows):
    for i_col in range(N_cols):
        if (N_rows==1) and (N_cols == 1):
            ax = axs
        elif N_rows == 1:
            ax = axs[i_col]
        elif N_cols == 1:
            ax = axs[i_row]
        else:
            ax = axs[i_row, i_col]

        # target
        tag_band = tag_bands[i_plot]
        label_band = label_bands[i_plot]
        i_plot += 1

        # plot the distribution
        xval_ori = cata_list[tag_band+'_ori'].values
        yval_ori = cata_list[tag_band+'_new'].values - cata_list[tag_band+'_ori'].values
        h = ax.hist2d(xval_ori, yval_ori, 
                bins=nbins, range=SRANGE, cmap=COLOR_MAP, 
                density=True)

        # labels
        LABEL = label_band
        ax.text(0.1, 0.85, LABEL, transform=ax.transAxes, fontsize=font_size)

        # zero lines
        ax.axhline(y=0, ls=':', linewidth=1.0, color='k')

        # quantile lines
        xbin_edges = np.linspace(XRANGE[0], XRANGE[1], 20)
        yval_low = []
        yval_high = []
        xval = (xbin_edges[1:] + xbin_edges[:-1])/2.
        for i_magbin in range(len(xbin_edges)-1):
            x_min = xbin_edges[i_magbin]
            x_max = xbin_edges[i_magbin + 1]
            mask_tmp = (xval_ori >= x_min) & (xval_ori < x_max)

            yval_low.append(np.nanpercentile(yval_ori[mask_tmp], 16))
            yval_high.append(np.nanpercentile(yval_ori[mask_tmp], 84))

        tmp = ax.plot(xval, yval_low, 
                        color=line_colour, linestyle='--', linewidth=1.5)
        tmp = ax.plot(xval, yval_high, 
                        color=line_colour, linestyle='--', linewidth=1.5)

        # some general setting
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_xlim(XRANGE[0], XRANGE[1])
        ax.set_ylim(YRANGE[0], YRANGE[1])

fig.text(0.5, 0.01, XLABEL, ha='center', fontsize=22)
fig.text(0.01, 0.5, YLABEL, va='center', rotation='vertical', fontsize=22)

fig.tight_layout()

if outpath == 'show':
    plt.show()
    plt.close()
else:
    plt.savefig(outpath, dpi=300)
    plt.close()
    plt.switch_backend(backend_orig)
    print("2D histogram plot saved as", outpath)    



