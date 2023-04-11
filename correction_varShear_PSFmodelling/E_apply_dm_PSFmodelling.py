# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-09-07 09:08:02
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-11 17:04:11

### apply the dm from PSF modelling to the whole sample results

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ general info

# the dm caused by PSF modelling
inpath_dm = './results/dm_PSFmodelling_41_rewei.csv'

# the m results after varCorr
inpath_m = './results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_varCorr_dz0p1_2D.csv'

# where to save
outpath = os.path.join('./results', os.path.basename(inpath_m).replace('.csv', '_PSFmodellingCorr.csv'))

# +++++++++++++++++++++++++++++ load catalogues
# load m results
cata_m = pd.read_csv(inpath_m)
cata_dm = pd.read_csv(inpath_dm)

# +++++++++++++++++++++++++++ correct the m 
for col_m in ['m1', 'm2']:
    cata_m.loc[:, col_m] += cata_dm.loc[:, col_m]
    cata_m.loc[:, f'{col_m}_err'] = (cata_m.loc[:, f'{col_m}_err']**2 
                                    + cata_dm.loc[:, f'{col_m}_err']**2
                                    )**0.5  
# save
print(cata_m)
cata_m.to_csv(outpath, index=False)
print('saved to', outpath)

# # ############# out info
#          m1        m2        c1  ...  Z_B_min  Z_B_max  mean_blendFraction
# 0  0.001124  0.005683  0.000317  ...      0.1      2.0            0.353972
# 1 -0.011073 -0.015410  0.000100  ...      0.1      0.3            0.345218
# 2 -0.023123 -0.013834 -0.000007  ...      0.3      0.5            0.331529
# 3 -0.009403 -0.008299  0.000313  ...      0.5      0.7            0.364396
# 4  0.015380  0.023772  0.000608  ...      0.7      0.9            0.365923
# 5  0.030848  0.036118  0.000686  ...      0.9      1.2            0.370566
# 6  0.072979  0.068929  0.000429  ...      1.2      2.0            0.358490

# [7 rows x 11 columns]
# saved to ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_varCorr_dz0p1_2D_PSFmodellingCorr.csv
