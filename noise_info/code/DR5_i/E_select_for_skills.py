# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-01-20 09:02:36
# @Last Modified by:   lshuns
# @Last Modified time: 2023-01-20 09:21:09

### select tiles for SKiLLS simulations
###     following DR4 selections

import os
import re
import glob

import numpy as np
import pandas as pd

# >>>>>>>>>>>> I/O

# where to find DR4 tiles
refDir = '../../skills_fiducial'
refPrefix = 'noise_selec'

# where to save
outDir = '../../skills_fiducial_dr5_iS/'
outPrefix = 'noise_selec_dr5_iS'
if not os.path.isdir(outDir):
    os.mkdir(outDir)

# where to find all DR5 tiles
inpath = '../../test_noise_dr5_iS.csv'

# >>>>>>>>>>> workhorse

## load all available noise tiles
tile_info = pd.read_csv(inpath)
print('Total tiles', len(tile_info))

## find all skills tiles
ref_path_list = glob.glob(os.path.join(refDir, f'{refPrefix}_*_part*.csv'))
print('number of parts found', len(ref_path_list))

## select 
tile_selec_all = []
for ref_path in ref_path_list:

    ref_point = re.search(refPrefix + r'_(.*)_N', os.path.basename(ref_path))[1]
    N_tiles = re.search(r'N(\d+)', os.path.basename(ref_path))[1]
    ref_part = re.search(r'part(\d+)', os.path.basename(ref_path))[1]

    # get the tile labels
    tile_lables = pd.read_csv(ref_path)['label'].values
    print(f'Number tiles in filename {N_tiles}')
    print(f'Number tiles found {len(tile_lables)}')

    # select the tiles
    tile_selec = pd.DataFrame({'label': tile_lables})
    tile_selec = tile_selec.merge(tile_info, on='label', how='left')

    # save
    outpath = os.path.join(outDir, f'{outPrefix}_{ref_point}_N{N_tiles}_part{ref_part}.csv')
    tile_selec.to_csv(outpath, index=False)
    tile_selec_all.append(tile_selec)
    del tile_selec
    print(f'saved as {outpath}')
# combine to one file
tile_selec_all = pd.concat(tile_selec_all)
outpath = os.path.join(outDir, f'{outPrefix}_combined.csv')
tile_selec_all.to_csv(outpath, index=False)
print(f'saved as {outpath}')

# Total tiles 1346
# number of parts found 6
# Number tiles in filename 18
# Number tiles found 18
# saved as ../../skills_fiducial_dr5_iS/noise_selec_dr5_iS_45_-31_N18_part4.csv
# Number tiles in filename 18
# Number tiles found 18
# saved as ../../skills_fiducial_dr5_iS/noise_selec_dr5_iS_135_0_N18_part0.csv
# Number tiles in filename 18
# Number tiles found 18
# saved as ../../skills_fiducial_dr5_iS/noise_selec_dr5_iS_335_-31_N18_part5.csv
# Number tiles in filename 18
# Number tiles found 18
# saved as ../../skills_fiducial_dr5_iS/noise_selec_dr5_iS_230_0_N18_part2.csv
# Number tiles in filename 18
# Number tiles found 18
# saved as ../../skills_fiducial_dr5_iS/noise_selec_dr5_iS_10_-31_N18_part3.csv
# Number tiles in filename 18
# Number tiles found 18
# saved as ../../skills_fiducial_dr5_iS/noise_selec_dr5_iS_180_0_N18_part1.csv
# saved as ../../skills_fiducial_dr5_iS/noise_selec_dr5_iS_combined.csv