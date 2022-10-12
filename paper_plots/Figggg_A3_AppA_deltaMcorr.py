# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-12-21 17:08:51
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-25 22:12:26

### Appendix A: magnitude modification factor

import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

from matplotlib.ticker import AutoMinorLocator, LogLocator

from scipy import stats

from astropy.io import fits

# ++++++++++++++ I/O

# SURFS cata with Dmag_corr
area_surfs = 108.
infile_surfs = '/disks/shear10/ssli/ImSim/input/SURFS_cata/SURFS_SHARK_LC_100sqdeg_testcols_calib.feather'
surfs_cata = pd.read_feather(infile_surfs)
print('surfs ', len(surfs_cata))

# ++++++++++++++ plot

# font size
plt.rc('font', size=16)
# tex
plt.rcParams["text.usetex"] = True

# redshift
z_surfs = surfs_cata['zobs'].values
zbin_edges = np.arange(0.1, 2.5, 0.1)

# log stellar mass
logmass_surfs = np.log10((surfs_cata['mstars_bulge'] + surfs_cata['mstars_disk'])/0.7)
logmass_edges = np.append(np.arange(8., 12, 0.2), 12)

# Dmag_corr
dmag_list = surfs_cata['Dmag_corr'].values
dmag_matrix = np.zeros((len(logmass_edges)-1, len(zbin_edges)-1))
for i_zbin in range(len(zbin_edges)-1):
    zbin_min = zbin_edges[i_zbin]
    zbin_max = zbin_edges[i_zbin+1]

    print(f'redshift bin: {zbin_min}, {zbin_max}')

    for i_Mbin in range(len(logmass_edges)-1):

        Mbin_min = logmass_edges[i_Mbin]
        Mbin_max = logmass_edges[i_Mbin+1]

        print(f'mass bin: {Mbin_min}, {Mbin_max}')

        mask_surfs = (z_surfs>=zbin_min) & (z_surfs<zbin_max) & (logmass_surfs>=Mbin_min) & (logmass_surfs<Mbin_max)
        try:
            dmag_matrix[i_Mbin, i_zbin] = dmag_list[mask_surfs][0]
        except IndexError:
            pass

fig, ax = plt.subplots()
count_scale = [-1, 1]
X, Y = np.meshgrid(zbin_edges, logmass_edges)
qm = plt.pcolormesh(X, Y, dmag_matrix, vmin=count_scale[0], vmax=count_scale[1], cmap='RdBu')
cbar = plt.colorbar(qm, ticks=[-1, -0.5, 0, 0.5, 1])
cbar.ax.set_ylabel(r'$\Delta {\rm mag}$', rotation=270)
plt.xlabel(r'Redshift')
plt.ylabel(r'${\rm log}_{10}~({\rm M}_{\star}/{\rm M}_{\odot})$')

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

plt.xticks([0.5, 1.0, 1.5, 2.0])
plt.yticks([8, 9, 10, 11, 12])

# plt.show()
plt.savefig('./plots/Upsilon_corr.pdf', dpi=300)