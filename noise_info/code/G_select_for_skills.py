# @Author: ssli
# @Date:   2021-04-24, 14:57:05
# @Last modified by:   ssli
# @Last modified time: 2021-05-18, 16:10:47

### select tiles for SKiLLS simulations
###     Total number of pointing: 108 
###     6 contiguous regions with variable star density (each with area: 108/6 = 18 sqdeg)
###     pointing missing chip psf info is discard (based on lensfit psf modelling: require starnum > 10)

import os

import numpy as np
import pandas as pd

# out info
outdir = '../skills_v034/'
out_prefix = 'noise_selec'
if not os.path.isdir(outdir):
    os.mkdir(outdir)

# for psf selection
min_starnum = 10

# contiguous pointing
## associated with star density
ra_stars = [135, 180, 230, 10, 45, 360-25]
dec_stars = [0, 0, 0, -31, -31, -31]

# number of tiles in each point
N_tiles = [18, 18, 18, 18, 18, 18]

# all available noise tiles
tile_info = pd.read_csv('../kids_dr4_noise_ugriZYJHKs.csv')
print('Total tiles', len(tile_info))
## associated chip psfs
chip_psf_info = pd.read_csv('../kids_dr4_psf_moffat_fromcoeffs.csv')
### discard those not meet selection requirement
chip_psf_info = chip_psf_info[chip_psf_info['starnum']>=min_starnum]
chip_psf_info.reset_index(drop=True, inplace=True)

# get the sky positions
ra_dec = tile_info['label'].str.split(r"_", expand=True)
ra_dec.columns = ['ra', 'dec']
tile_info = pd.concat([tile_info, ra_dec], axis=1)
tile_info = tile_info.astype({'ra': 'float', 'dec': 'float'})

tile_selec_all = []
for i_part, ra_center in enumerate(ra_stars):

    dec_center = dec_stars[i_part]
    N_tile = N_tiles[i_part]

    # distance to the center
    tile_info.loc[:, 'dist2center'] = np.hypot(np.array(tile_info['ra']-ra_center), np.array(tile_info['dec']-dec_center))
    ## sort tiles based on distance
    tile_info_sorted = tile_info.sort_values(by=['dist2center'])
    tile_info_sorted.reset_index(drop=True, inplace=True)

    # check each tile to make sure it contains enough chips
    tile_selec = pd.DataFrame(columns=tile_info_sorted.columns)
    i_row = 0
    for _, tile_info_row in tile_info_sorted.iterrows():

        tile_lable = tile_info_row['label']
        # get corresponding chip psf info
        mask_chips = chip_psf_info['label'] == tile_lable
        ## should contain: 5 * 32 = 160 chips 
        if np.sum(mask_chips)==160:
            tile_selec.loc[i_row] = tile_info_row
            i_row += 1
        ## break
        if i_row == N_tile:
            break

    # save
    tile_selec.drop(columns=['ra', 'dec'], inplace=True)
    print('final selected cata', len(tile_selec))
    outpath = os.path.join(outdir, f'{out_prefix}_{ra_center}_{dec_center}_N{N_tile}_part{i_part}.csv')
    tile_selec.to_csv(outpath, index=False)
    tile_selec_all.append(tile_selec)
    print(f'saved as {outpath}')
# combine to one file
tile_selec_all = pd.concat(tile_selec_all)
outpath = os.path.join(outdir, f'{out_prefix}_combined.csv')
tile_selec_all.to_csv(outpath, index=False)
print(f'saved as {outpath}')

# Total tiles 979
# final selected cata 18
# saved as ../skills_v034/noise_selec_135_0_N18_part0.csv
# final selected cata 18
# saved as ../skills_v034/noise_selec_180_0_N18_part1.csv
# final selected cata 18
# saved as ../skills_v034/noise_selec_230_0_N18_part2.csv
# final selected cata 18
# saved as ../skills_v034/noise_selec_10_-31_N18_part3.csv
# final selected cata 18
# saved as ../skills_v034/noise_selec_45_-31_N18_part4.csv
# final selected cata 18
# saved as ../skills_v034/noise_selec_335_-31_N18_part5.csv
# saved as ../skills_v034/noise_selec_combined.csv