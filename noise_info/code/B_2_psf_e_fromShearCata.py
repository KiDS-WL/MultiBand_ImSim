# @Author: lshuns
# @Date:   2021-05-06, 13:22:00
# @Last modified by:   ssli
# @Last modified time: 2021-05-11, 15:28:05

### extract PSF ellipticity from the shear catalogue

import time
start_time = time.time()

import pandas as pd
import numpy as np

from astropy.io import fits

# +++++++++++++ I/O
## the public Shear catalogue
infile = '/disks/shear15/KiDS/KiDS-1000/K1000-SHEAR-9bandPZ-CATALOGUES/KiDS_DR4.1_ugriZYJHKs_SOM_gold_WL_cat.fits'

## where to save results
outpath = '../test_psf_e_fromShearCata.csv'

# +++++++++++++ general settings

## column name for the tile label
col_tile = 'THELI_NAME'

# +++++++++++++ workhorse
# data
with fits.open(infile) as hdul:
    data = hdul[1].data

# labels
tile_labels = np.unique(data[col_tile])
n_tiles = len(tile_labels)
print('Number of tiles', n_tiles)

# initialize the catalogue
columns = ['label', 'N_sources', 'PSF_e_mean', 'PSF_e_median', 'PSF_e_p16', 'PSF_e_p84'] + [f'PSF_e_mean_expo{i_exp}' for i_exp in range(5)] + [f'PSF_e_median_expo{i_exp}' for i_exp in range(5)] + [f'PSF_e_p16_expo{i_exp}' for i_exp in range(5)] + [f'PSF_e_p84_expo{i_exp}' for i_exp in range(5)]
cata_final = pd.DataFrame(index = np.arange(n_tiles), columns=columns)
cata_final.loc[:, 'label'] = tile_labels

# loop over all tiles
for i_tile, tile_label in enumerate(tile_labels):

    print('Working on tile:', tile_label, i_tile)
    data_selec = data[(data[col_tile] == tile_label)]
    N_sources = len(data_selec)
    print('Number of sources:', N_sources)
    cata_final.loc[i_tile, 'N_sources'] = N_sources

    # loop over exposures
    for i_exp in range(5):

        # mean value for all exposures
        e_tmp = np.hypot(data_selec[f'PSF_e1'], data_selec[f'PSF_e2'])
        cata_final.loc[i_tile, f'PSF_e_mean'] = np.nanmean(e_tmp)
        cata_final.loc[i_tile, f'PSF_e_median'] = np.nanmedian(e_tmp)
        cata_final.loc[i_tile, f'PSF_e_p16'] = np.nanpercentile(e_tmp, 16)
        cata_final.loc[i_tile, f'PSF_e_p84'] = np.nanpercentile(e_tmp, 84)

        # value for each exposure
        ### only look into those with meaningful values
        mask_tmp = (data_selec[f'PSF_e1_exp{i_exp+1}']<99.) & (data_selec[f'PSF_e2_exp{i_exp+1}']<99.)
        e_tmp = np.hypot(data_selec[mask_tmp][f'PSF_e1_exp{i_exp+1}'], data_selec[mask_tmp][f'PSF_e2_exp{i_exp+1}'])
        cata_final.loc[i_tile, f'PSF_e_mean_expo{i_exp}'] = np.nanmean(e_tmp)
        cata_final.loc[i_tile, f'PSF_e_median_expo{i_exp}'] = np.nanmedian(e_tmp)
        cata_final.loc[i_tile, f'PSF_e_p16_expo{i_exp}'] = np.nanpercentile(e_tmp, 16)
        cata_final.loc[i_tile, f'PSF_e_p84_expo{i_exp}'] = np.nanpercentile(e_tmp, 84)

cata_final.to_csv(outpath, index=False)
print('psf e info saved in', outpath)
print('All finished in', (time.time()-start_time)/3600., 'h')
