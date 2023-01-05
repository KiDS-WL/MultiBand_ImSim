# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-09-07 09:08:02
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-30 13:21:21

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

# the m results from PSF modelling
inpath_dm = '../sensitivity_test/PSFmodelling/outputs/dm_ZB_test_fiducial_rewei_41.csv'

# the m results after varCorr
inpath_m = './outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_varCorr_dz0p1_2D.csv'

# where to save
outpath = os.path.join('./outputs', os.path.basename(inpath_m).replace('.csv', '_PSFmodelling.csv'))

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
cata_m.to_csv(outpath, index=False)
print('saved to', outpath)

# saved to ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_varCorr_dz0p1_2D_PSFmodelling.csv
