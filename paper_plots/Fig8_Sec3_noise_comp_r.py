# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-06-27 09:23:48
# @Last Modified by:   lshuns
# @Last Modified time: 2022-07-06 16:47:20

### PSF and noise info in r band

import math
import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator, LogLocator
from matplotlib.patches import Rectangle

# ++++++++++++++++++++++++++++++++ I/O

## >>> full noise file 
### PSF
infile = '/disks/shear16/ssli/ImSim/input/KiDS_PSF_ima/r_tukey/PSFinfo_summarised.csv'
data_full_PSF = pd.read_csv(infile)
### noise
infile = '../noise_info/kids_dr4_noise_ugriZYJHKs.csv'
data_full_noise = pd.read_csv(infile)

## >>> selected noise file
infile = '../noise_info/skills_fiducial/noise_selec_combined.csv'
data_selec = pd.read_csv(infile)

# ++++++++++++++++++++++++++++++++ get PSF info

# the full cata
PSFsigma_r_all = (data_full_PSF['moment_size_pixel'].values)**0.5 * 0.214 # arcsec
PSFe_r_all = np.hypot(data_full_PSF['moment_e1'].values, data_full_PSF['moment_e2'].values)
print('total PSF models', len(PSFsigma_r_all))

# the selected cata
all_tiles = (data_full_PSF['file'].str.extract(r'tile(.*)_bandr').values).flatten()
mask_selected = np.isin(all_tiles, data_selec['label'].values)
PSFsigma_r_selected = PSFsigma_r_all[mask_selected]
PSFe_r_selected = PSFe_r_all[mask_selected]
del all_tiles, mask_selected, data_full_PSF
print('selected PSF models', len(PSFsigma_r_selected))

# ++++++++++++++++++++++++++++++++ plot

N_plots = 3

paras_list = [
            [data_full_noise['rmsExpo_r'].values, data_selec['rmsExpo_r'].values],
            [PSFsigma_r_all, PSFsigma_r_selected], 
            [PSFe_r_all, PSFe_r_selected]
            ]

COLORs_list = [
            ['r', 'b'],
            ['r', 'b'],
            ['r', 'b']
                ]
 
XLABEL_list = ['Pixel noise (RMS value)', 'PSF size [arcsec]', r'PSF $\epsilon$']

xrange_list = [[4, 8], [0.2, 0.4], [0.0, 0.1]]
xtick_spe_list = [  [4, 6, 8],
                    [0.2, 0.3, 0.4],
                    [0.0, 0.05, 0.1]
                ]

yrange_list = [[0, 0.9], [0, 16], [0, 40]]
ytick_spe_list = [  [0.3, 0.6],
                    [4, 8, 12], 
                    [10, 20, 30]
                    ]

LABELs = ['KiDS DR4', 'SKiLLS fiducial']
LABELs_c = ['r', 'b']

###### settings
outpath = 'show'
outpath = './plots/skills_bandr.pdf'

font_size = 14
usetex = True

nbins = 20
DENSITY = True
YLABEL = 'Probability density'

###### plotting
# font size
plt.rc('font', size=font_size)
# tex
plt.rcParams["text.usetex"] = usetex

if outpath != 'show':
    backend_orig = plt.get_backend()
    plt.switch_backend("agg")

N_rows = math.ceil(N_plots**0.5)
N_cols = math.ceil(N_plots/N_rows)
fig, axs = plt.subplots(N_rows, N_cols)

fig.subplots_adjust(hspace=0.35)
fig.subplots_adjust(wspace=0.2)

i_plot = 0
for i_row in range(N_rows):
    for i_col in range(N_cols):
        if i_plot >= N_plots:
            if N_rows == 1:
                axs[i_col].axis('off')
            elif N_cols == 1:
                axs[i_row].axis('off')
            else:
                axs[i_row, i_col].axis('off')
        else:
            if (N_rows==1) and (N_cols == 1):
                ax = axs
            elif N_rows == 1:
                ax = axs[i_col]
            elif N_cols == 1:
                ax = axs[i_row]
            else:
                ax = axs[i_row, i_col]

            paras = paras_list[i_plot]
            COLORs = COLORs_list[i_plot]

            XRANGE = xrange_list[i_plot]
            YRANGE = yrange_list[i_plot]

            XLABEL = XLABEL_list[i_plot]

            xtick_spe = xtick_spe_list[i_plot]
            ytick_spe = ytick_spe_list[i_plot]

            ax.hist(x=paras, bins=nbins, range=XRANGE, density=DENSITY, color=COLORs, histtype='step')

            ax.set_xlim(XRANGE[0], XRANGE[1])
            if YRANGE is not None:
                ax.set_ylim(YRANGE[0], YRANGE[1])

            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())

            if xtick_spe is not None:
                ax.set_xticks(xtick_spe)
            if ytick_spe is not None:
                ax.set_yticks(ytick_spe)

            ax.set_xlabel(XLABEL, fontsize=font_size)

        i_plot +=1

fig.text(0.02, 0.5, YLABEL, va='center', rotation='vertical')

if (LABELs is not None):

    handles = [Rectangle((0,0),1,1,color='white', ec=c) for c in LABELs_c]

    fig.legend(handles, LABELs, 
        loc = 'center', ncol=1,
        bbox_to_anchor=(0.72, 0.28), fancybox=True, shadow=True, fontsize=font_size)

if outpath == 'show':
    plt.show()
    plt.close()
else:
    plt.savefig(outpath, dpi=300)
    plt.close()
    plt.switch_backend(backend_orig)
    print("Histogram plot saved as", outpath)
