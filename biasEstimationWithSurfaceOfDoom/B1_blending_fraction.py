# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-05-01 16:31:58
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-03 15:17:59

### get the blending fraction

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ general info

# the surface of doom
inpath_doom = './SurfaceOfDoom/m_surface_ZB0p1_SNR20_R20_skills_v07D7p1_LF_321_kidsPhotometry_gold.csv'

# the fiducial SKiLLS-gold catalogue
inpath = '/disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat.A2.feather'
col_weight = 'AlphaRecalC_weight'

# the blending-only SKiLLS inputs
parts = ['part0', 'part1', 'part2', 'part3', 'part4', 'part5']
inpath_blend_list = [f'/disks/shear10/ssli/ImSim/input/SURFS_cata/varShear/skills_v07Ds_input_{part}_blended4_magCut25_4ImSim.feather'
                            for part in parts]

# ### running info
# number in blending-only cata 4605425
# number ori 28885361
# number within doom 28885361
# results saved back to ./SurfaceOfDoom/m_surface_ZB0p1_SNR20_R20_skills_v07D7p1_LF_321_kidsPhotometry_gold.csv
# Elapsed:24:07.54,User=1026.685,System=543.857,CPU=108.4%.

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
cata_final = cata_final[['id_input', col_weight, 'Z_B', 'SNR_LF_r', 'R']]
cata_final.rename(columns={'SNR_LF_r': "SNR"}, inplace=True)
## weight selection
cata_final = cata_final[cata_final[col_weight]>0].copy()
cata_final.reset_index(drop=True, inplace=True)
print('number ori', len(cata_final))
## for saving binning info
cata_final.loc[:, 'binZB_id'] = -999
cata_final.loc[:, 'binSNR_id'] = -999
cata_final.loc[:, 'binR_id'] = -999

# load surface of doom
mc_surface = pd.read_csv(inpath_doom)
## for saving fraction
mc_surface.loc[:, 'WeiTot'] = -999
mc_surface.loc[:, 'WeiBlend'] = -999

# +++++++++++++++++++++++++++ get fraction

# bin galaxies
Z_B_edges = np.unique(mc_surface['binZB_min']).tolist() + [np.max(mc_surface['binZB_max'])]
cata_final.loc[:, 'binZB_id'] = pd.cut(cata_final['Z_B'].values, Z_B_edges, 
                                right=True, labels=False)
for i_zbin in range(len(Z_B_edges)-1):

    SNR_edges = np.unique(mc_surface.loc[(mc_surface['binZB_id']==i_zbin), 
                    'binSNR_min']).tolist() \
                    + [np.max(mc_surface.loc[(mc_surface['binZB_id']==i_zbin), 
                        'binSNR_max'])]

    mask_binZB = cata_final['binZB_id'].values == i_zbin
    cata_final.loc[mask_binZB, 'binSNR_id'] = pd.cut(cata_final.loc[mask_binZB, 'SNR'].values, SNR_edges, 
                                right=True, labels=False)

    for i_SNR in range(len(SNR_edges)-1):

        R_edges = np.unique(mc_surface.loc[(mc_surface['binZB_id']==i_zbin)&(mc_surface['binSNR_id']==i_SNR), 
                    'binR_min']).tolist() \
                    + [np.max(mc_surface.loc[(mc_surface['binZB_id']==i_zbin)&(mc_surface['binSNR_id']==i_SNR), 
                        'binR_max'])]

        mask_binSNR = cata_final['binSNR_id'].values == i_SNR
        cata_final.loc[mask_binZB&mask_binSNR, 'binR_id'] = pd.cut(
                                    cata_final.loc[mask_binZB&mask_binSNR, 'R'].values, 
                                    R_edges, 
                                    right=True, labels=False)
        del mask_binSNR
    del mask_binZB

# group
## drop -999 bins
cata_final = cata_final[(cata_final['binZB_id']>-999)&(cata_final['binSNR_id']>-999)&(cata_final['binR_id']>-999)].copy()
cata_final.reset_index(drop=True, inplace=True)
print('number within doom', len(cata_final))
## sort to speed up
cata_final = cata_final.astype({'binZB_id': int, 'binSNR_id': int, 'binR_id': int})
cata_final.sort_values(by=['binZB_id', 'binSNR_id', 'binR_id'], inplace=True)
cata_final = cata_final.groupby(by=['binZB_id', 'binSNR_id', 'binR_id'])

# loop over groups to get fraction
for name, group in cata_final:
    binZB_id, binSNR_id, binR_id = name

    # total weights in each bin
    Wei0 = np.sum(group[col_weight].values)
    ### get blending
    mask_blend = np.isin(group['id_input'].values, id_blending)
    Wei = np.sum(group.loc[mask_blend, col_weight].values)
    del mask_blend, group
    ### get ratio
    ratio_blend = Wei / Wei0
    # print('Weights', Wei0, Wei, ratio_blend)

    # save results
    mask_doom = (mc_surface['binZB_id']==binZB_id)\
            &(mc_surface['binSNR_id']==binSNR_id)\
            &(mc_surface['binR_id']==binR_id)
    mc_surface.loc[mask_doom, 'WeiTot'] = Wei0
    mc_surface.loc[mask_doom, 'WeiBlend'] = Wei

# save back
# print('>>> sorted WeiTot', np.sort(mc_surface['WeiTot'].values))
# print('>>> sorted WeiBlend', np.sort(mc_surface['WeiBlend'].values))
# print('>>> sorted fraction', np.sort(mc_surface['WeiBlend'].values/mc_surface['WeiTot'].values))
mc_surface.to_csv(inpath_doom, index=False, float_format='%.6f')
print(f'results saved back to {inpath_doom}')