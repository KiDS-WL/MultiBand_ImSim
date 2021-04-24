# @Author: lshuns
# @Date:   2021-04-19, 14:45:54
# @Last modified by:   lshuns
# @Last modified time: 2021-04-22, 11:30:21

### empirically calibrate SURFS photometry by learning Mstar/L from COSMOS15
###     reference: Laigle et al. 2016

import os
import time

import numpy as np

from astropy.io import fits
from astropy.table import Table

# ++++++++++++++ I/O

# target files
indir = '/disks/shear10/ssli/ImSim/input/SURFS_cata/'
file_list = [f'SURFS_SHARK_LC_patch{i_patch}.fits' for i_patch in range(4)]

# COSMOS15 cata
infile_cosmos = '/disks/shear10/ssli/ImSim/input/COSMOS_cata/outside/cosmos_2015_dr_2_1_138deg2_galONLY.fits'

# ++++++++++++++ general settings

# magnitude columns used to calculate the luminosity
mag_cosmos_col = 'm_k'
mag_simu_col = 'K_VISTA_absolute'

# absolute magnitude of the Sun in chosen band
mag_sun = 5.08

# col name for saving magnitude modification factor
col_name_final = 'Upsilon_corr'

# binning in z and log stellar mass
zbin_edges = np.append(np.arange(0., 2.5, 0.1), 2.5)
logmass_edges = np.append(np.arange(8., 12, 0.2), 12)

# ++++++++++++++ workhorse
start_time = time.time()
# the cosmos catalogue
with fits.open(infile_cosmos) as hdul:
    cosmos_cata = hdul[1].data
print(f'COSMOS15 cata loaded from {infile_cosmos}')
## used columns
z_cosmos = cosmos_cata['photoz']
logmass_cosmos = cosmos_cata[f'mass_med'] # log stellar mass
logL_cosmos = -0.4*(cosmos_cata[mag_cosmos_col]-mag_sun) # solar L
Upsilon_cosmos = 10**(logmass_cosmos - logL_cosmos)
del logL_cosmos
del cosmos_cata

# loop over all target catalogues
for file in file_list:

    infile_tmp = os.path.join(indir, file)
    with fits.open(infile_tmp) as hdul:
        surfs_cata = Table(hdul[1].data).to_pandas()
    print(f'working on {file}')

    # for saving
    surfs_cata.loc[:, col_name_final] = 0

    # redshift
    z_surfs = np.array(surfs_cata['zobs'])

    # log stellar mass
    logmass_surfs = np.log10((surfs_cata['mstars_bulge'] + surfs_cata['mstars_disk'])/0.7)

    # luminosity
    logL_surfs = -0.4*(np.array(surfs_cata[mag_simu_col])-mag_sun)

    # mass to light ratio
    Upsilon_surfs = 10**(logmass_surfs - logL_surfs)

    # calculate modification factor for each z and mass bin
    N_empty_cosmos = 0
    N_empty_surfs = 0
    for i_zbin in range(len(zbin_edges)-1):
        zbin_min = zbin_edges[i_zbin]
        zbin_max = zbin_edges[i_zbin+1]
        # print(f'redshift bin: {zbin_min:.1f}, {zbin_max:.1f}')

        for i_Mbin in range(len(logmass_edges)-1):
            Mbin_min = logmass_edges[i_Mbin]
            Mbin_max = logmass_edges[i_Mbin+1]
            # print(f'mass bin: {Mbin_min:.1f}, {Mbin_max:.1f}')

            # select cosmos
            mask_tmp = (z_cosmos>=zbin_min) & (z_cosmos<zbin_max) & (logmass_cosmos>=Mbin_min) & (logmass_cosmos<Mbin_max)
            Upsilon_cosmos_selec = Upsilon_cosmos[mask_tmp]

            # select surfs
            mask_surfs = (z_surfs>=zbin_min) & (z_surfs<zbin_max) & (logmass_surfs>=Mbin_min) & (logmass_surfs<Mbin_max)
            Upsilon_surfs_selec = Upsilon_surfs[mask_surfs]

            # calculate the modification factor
            if (len(Upsilon_cosmos_selec) > 0) and (len(Upsilon_surfs_selec) > 0):
                dUpsilon = np.nanmedian(Upsilon_surfs_selec) / np.nanmedian(Upsilon_cosmos_selec)
                dmag = -2.5*np.log10(dUpsilon)
                # print(f'modification factor: {dmag}')
                surfs_cata.loc[mask_surfs, col_name_final] = dmag
            elif (len(Upsilon_cosmos_selec) == 0):
                N_empty_cosmos += 1
                # print('Empty bin for cosmos')
            elif (len(Upsilon_surfs_selec) == 0):
                N_empty_surfs += 1
                # print('Empty bin for surfs')

    print('Number of empty bins for cosmos', N_empty_cosmos)
    print('Number of empty bins for surfs', N_empty_surfs)

    # save back to the original catalogue
    os.remove(infile_tmp)
    Table.from_pandas(surfs_cata).write(infile_tmp, format='fits')

print('All finished in', (time.time()-start_time)/60., 'min')

# COSMOS15 cata loaded from /disks/shear10/ssli/ImSim/input/COSMOS_cata/outside/cosmos_2015_dr_2_1_138deg2_galONLY.fits
# working on SURFS_SHARK_LC_patch0.fits
# Number of empty bins for cosmos 38
# Number of empty bins for surfs 10
# working on SURFS_SHARK_LC_patch1.fits
# Number of empty bins for cosmos 38
# Number of empty bins for surfs 10
# working on SURFS_SHARK_LC_patch2.fits
# Number of empty bins for cosmos 38
# Number of empty bins for surfs 11
# working on SURFS_SHARK_LC_patch3.fits
# Number of empty bins for cosmos 38
# Number of empty bins for surfs 10
# All finished in 40.26756705840429 min
