# @Author: ssli
# @Date:   2021-04-24, 14:57:05
# @Last modified by:   ssli
# @Last modified time: 2021-05-18, 16:10:47

### select tiles for cmousse_v01 photometric catalogues
###     based on the star pointing (108 sqdeg / 6)

import os

import numpy as np
import pandas as pd

# out info
outdir = '../skills_v03/'
out_prefix = 'noise_selec'

# Star points
ra_stars = [135, 180, 230, 10, 45, 360-25]
dec_stars = [0, 0, 0, -31, -31, -31]

# number of tiles in each point
N_tiles = [18, 18, 18, 18, 18, 18]

# all available noise tiles
tile_info = pd.read_csv('../kids_dr4_noise_ugriZYJHKs.csv')
print('Total tiles', len(tile_info))

# get the sky positions
ra_dec = tile_info['label'].str.split(r"_", expand=True)
ra_dec.columns = ['ra', 'dec']
tile_info = pd.concat([tile_info, ra_dec], axis=1)
tile_info = tile_info.astype({'ra': 'float', 'dec': 'float'})
tile_info = tile_info.sort_values(by=['ra','dec'])

tile_selec_all = []
for i_part, ra_center in enumerate(ra_stars):

    dec_center = dec_stars[i_part]
    N_tile = N_tiles[i_part]

    # distance to the center
    dist_tmp = np.hypot(np.array(tile_info['ra']-ra_center), np.array(tile_info['dec']-dec_center))

    # select the nearest tiles
    dist_max_tmp = np.sort(dist_tmp)[N_tile]
    mask_tmp = dist_tmp<=dist_max_tmp
    tile_selec = tile_info[mask_tmp].copy()
    tile_selec.drop(columns=['ra', 'dec'], inplace=True)
    print('original selected cata', len(tile_selec))
    tile_selec = tile_selec.sample(n=N_tile, random_state=940120)
    print('final selected cata', len(tile_selec))

    # save
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
# original selected cata 20
# final selected cata 18
# saved as ../skills_v02/noise_selec_135_0_N18_part0.csv
# original selected cata 20
# final selected cata 18
# saved as ../skills_v02/noise_selec_180_0_N18_part1.csv
# original selected cata 20
# final selected cata 18
# saved as ../skills_v02/noise_selec_230_0_N18_part2.csv
# original selected cata 19
# final selected cata 18
# saved as ../skills_v02/noise_selec_10_-31_N18_part3.csv
# original selected cata 19
# final selected cata 18
# saved as ../skills_v02/noise_selec_45_-31_N18_part4.csv
# original selected cata 19
# final selected cata 18
# saved as ../skills_v02/noise_selec_335_-31_N18_part5.csv
# saved as ../skills_v02/noise_selec_combined.csv
