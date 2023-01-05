# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-09-29 20:23:37
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-30 13:20:38

### apply the dm from varShear to the whole sample results
####### in 2D bins (R and SNR)

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ general info

# the ZB bins
ZB_mins = [-999, 0.1, 0.3, 0.5, 0.7, 0.9, 1.2]
ZB_maxs = [-999, 0.3, 0.5, 0.7, 0.9, 1.2, 2.0]

# the raw results
inpath_raw_list = ['../biasEstimation/outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_whole.csv']
inpath_raw_list += [f'../biasEstimation/outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin{i_bin}.csv'
                    for i_bin in range(6)]

# the blending-only correction
inpath_blend_list = ['./outputs/dm_forBlending_dz0p1_whole.csv']
inpath_blend_list += [f'./outputs/dm_forBlending_dz0p1_bin{i_bin}.csv' for i_bin in range(6)]

# the blending fraction
inpath_fraction_list = ['./outputs/blending_fraction_whole.csv']
inpath_fraction_list += [f'./outputs/blending_fraction_bin{i_bin}.csv' for i_bin in range(6)]

# where to save
outpath = './outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_varCorr_dz0p1_2D.csv'

# +++++++++++++++++++++++++++ correct the m 

res_df = pd.DataFrame(-999, 
                index = np.arange(7), 
                columns = ['m1', 'm2', 'c1', 'c2',
                        'm1_err', 'm2_err', 'c1_err', 'c2_err',
                        'Z_B_min', 'Z_B_max', 'mean_blendFraction'])

for i_path, inpath_raw in enumerate(inpath_raw_list):

    # the range info
    res_df.loc[i_path, 'Z_B_min'] = ZB_mins[i_path]
    res_df.loc[i_path, 'Z_B_max'] = ZB_maxs[i_path]

    inpath_blend = inpath_blend_list[i_path]
    inpath_fraction = inpath_fraction_list[i_path]

    cata_raw = pd.read_csv(inpath_raw)
    cata_blend = pd.read_csv(inpath_blend)
    cata_fraction = pd.read_csv(inpath_fraction)

    # sort for easy combine
    cata_raw.sort_values(by=['bin1_id', 'bin2_id'], inplace=True, ignore_index=True)
    cata_blend.sort_values(by=['bin1_id', 'bin2_id'], inplace=True, ignore_index=True)
    cata_fraction.sort_values(by=['bin1_id', 'bin2_id'], inplace=True, ignore_index=True)

    # results after correction
    ratio_blend = cata_fraction['WeiBlend'].values/cata_fraction['WeiTot'].values
    for col_m in ['m1', 'm2']:
        res_df.loc[i_path, col_m] = np.average(cata_raw[col_m].values \
                                        + cata_blend[col_m].values * ratio_blend,
                                        weights=cata_raw['dataWei'].values)
        res_df.loc[i_path, f'{col_m}_err'] = np.sqrt(
                                                np.sum(
                                                    np.square(cata_raw['dataWei'].values \
                                                        * np.hypot(cata_raw[f'{col_m}_err'].values, 
                                                            cata_blend[f'{col_m}_err'].values * ratio_blend)
                                                        )
                                                    )
                                                ) / np.sum(cata_raw['dataWei'].values)
    ## c is not changed
    for col_c in ['c1', 'c2']:
        res_df.loc[i_path, col_c] = np.average(cata_raw[col_c].values,
                                        weights=cata_raw['dataWei'].values)
        res_df.loc[i_path, f'{col_c}_err'] = np.sqrt(
                                                np.sum(
                                                    np.square(cata_raw['dataWei'].values \
                                                        * cata_raw[f'{col_c}_err'].values
                                                        )
                                                    )
                                                ) / np.sum(cata_raw['dataWei'].values)
    ## the average blending fraction
    res_df.loc[i_path, 'mean_blendFraction'] = np.average(ratio_blend,
                                        weights=cata_raw['simWei'].values)

# save
print(res_df)
res_df.to_csv(outpath, index=False)
print('saved to', outpath)

############# out info
#          m1        m2        c1  ...  Z_B_min  Z_B_max  mean_blendFraction
# 0 -0.002906  0.002916  0.000322  ...   -999.0   -999.0            0.353935
# 1 -0.011517 -0.017842  0.000111  ...      0.1      0.3            0.345336
# 2 -0.026941 -0.017456 -0.000026  ...      0.3      0.5            0.331540
# 3 -0.012706 -0.010971  0.000319  ...      0.5      0.7            0.364526
# 4  0.009503  0.021798  0.000622  ...      0.7      0.9            0.365566
# 5  0.024692  0.032705  0.000681  ...      0.9      1.2            0.370393
# 6  0.066028  0.062894  0.000463  ...      1.2      2.0            0.358295

# [7 rows x 11 columns]
# saved to ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_varCorr_dz0p1_2D.csv
# Elapsed:0:02.55,User=3.159,System=7.530,CPU=418.8%.
