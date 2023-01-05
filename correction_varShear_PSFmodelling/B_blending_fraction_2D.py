# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-09-29 20:43:26
# @Last Modified by:   lshuns
# @Last Modified time: 2023-01-05 17:46:25

### get the blending fraction in 2D surface of dome

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ general info

# the whole sample
inpath = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut.feather'
col_weight = 'oldweight_LF_r'

# the binning bounds 
inpath_prefix_bounds = '../biasEstimation/outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut'
bin1_col = 'SNR_LF_r'
bin2_col = 'R'

# the bins
col_binning = 'Z_B'
binning_edges = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2, 2.0]

# parts
parts = ['part0', 'part1', 'part2', 'part3', 'part4', 'part5']

# the blending-only inputs
inpath_blend_list = [f'/disks/shear10/ssli/ImSim/input/SURFS_cata/varShear/skills_v07Ds_input_{part}_blended4_magCut25_4ImSim.feather'
                            for part in parts]

# where to save
outpath_prefix = './outputs/blending_fraction'

# +++++++++++++++++++++++++++++ load catalogues

# get the blended input ID
id_blending = []
for inpath_blend in inpath_blend_list:
    cata = pd.read_feather(inpath_blend)
    id_blending.append(cata['index'].values)
    del cata
id_blending = np.concatenate(id_blending)
print('number in blending-only cata', len(id_blending))

# load the catalogue
cata_final = pd.read_feather(inpath)
## used columns
cata_final = cata_final[['id_input', col_weight, col_binning, bin1_col, bin2_col]]
## weight selection
cata_final = cata_final[cata_final[col_weight]>0]
cata_final.reset_index(drop=True)

# +++++++++++++++++++++++++++ get fraction
for i_bin in range(len(binning_edges)):

    # the whole results
    if i_bin == 0:
        outpath = outpath_prefix + f'_whole.csv'
        bin_min = binning_edges[0]
        bin_max = binning_edges[-1]
        # 2D bins
        bin1_bounds = np.load(inpath_prefix_bounds + '_whole_bin1.npy')
        bin2_bounds = np.load(inpath_prefix_bounds + '_whole_bin2.npy')
    # others
    else:
        outpath = outpath_prefix + f'_bin{i_bin-1}.csv'
        bin_min = binning_edges[i_bin - 1]
        bin_max = binning_edges[i_bin]
        # 2D bins
        bin1_bounds = np.load(inpath_prefix_bounds + f'_bin{i_bin-1}_bin1.npy')
        bin2_bounds = np.load(inpath_prefix_bounds + f'_bin{i_bin-1}_bin2.npy')

    bin1_Nbins = len(bin1_bounds) - 1
    bin2_Nbins = len(bin2_bounds[0]) - 1
    print('>> Number of 2D bins', bin1_Nbins, bin2_Nbins)

    ## select
    cata_selec = cata_final[(cata_final[col_binning]>bin_min)&(cata_final[col_binning]<=bin_max)].copy()
    cata_selec.reset_index(drop=True, inplace=True)

    ## bin in the first col
    cata_selec.loc[:, 'bin1_id'] = pd.cut(cata_selec[bin1_col].values, bin1_bounds, 
                                right=True, labels=False)
    del bin1_bounds
    ## bin in the second col
    for ibin1, bin2_bound in enumerate(bin2_bounds):
        # mask in bin1
        mask_bin1 = cata_selec['bin1_id'].values == ibin1
        # bin in bin2
        cata_selec.loc[mask_bin1, 'bin2_id'] = \
                            pd.cut(cata_selec.loc[mask_bin1, bin2_col].values, 
                                        bin2_bound, 
                                        right=True, labels=False)
        del mask_bin1

    # group
    N0 = len(cata_selec)
    cata_selec.dropna(inplace=True)
    print('var remaining after 2D binning', len(cata_selec)/N0)
    cata_selec = cata_selec.astype({'bin1_id': int, 'bin2_id': int})
    cata_selec = cata_selec.groupby(by=['bin1_id', 'bin2_id'])

    # loop over bins and get fraction
    i_group = 0 
    surface = pd.DataFrame(-999, 
                        index = np.arange(bin1_Nbins*bin2_Nbins), 
                        columns = ['bin1_id', 'bin2_id', 
                                    'bin1_mean', 'bin2_mean',
                                    'WeiTot', 'WeiBlend'])
    for i in range(bin1_Nbins):
        for j in range(bin2_Nbins):

            # the name for group 
            name = (i, j)
            group = cata_selec.get_group(name)

            surface.loc[i_group, 'bin1_id'] = i
            surface.loc[i_group, 'bin2_id'] = j
            surface.loc[i_group, 'bin1_mean'] = np.average(group[bin1_col].values, weights=group[col_weight].values)
            surface.loc[i_group, 'bin2_mean'] = np.average(group[bin2_col].values, weights=group[col_weight].values)

            ### weights
            Wei0 = np.sum(group[col_weight].values)
            ### get blending
            mask_blend = np.isin(group['id_input'].values, id_blending)
            Wei = np.sum(group.loc[mask_blend, col_weight].values)
            del mask_blend, group
            ### get ratio
            ratio_blend = Wei / Wei0
            print('Weights', Wei0, Wei, ratio_blend)

            surface.loc[i_group, 'WeiTot'] = Wei0
            surface.loc[i_group, 'WeiBlend'] = Wei
            i_group += 1

    # save
    surface.to_csv(outpath, index=False)
    del surface
    print('saved to', outpath)

# ...
# saved to ./outputs/blending_fraction_bin5.csv
# Elapsed:27:57.97,User=992.120,System=523.290,CPU=90.3%.
