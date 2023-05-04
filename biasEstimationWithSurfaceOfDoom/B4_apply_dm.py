# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-05-03 12:19:36
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-03 13:30:07

### apply the dm from varShear and PSF modelling 

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# the surface of doom
inpath_doom = './SurfaceOfDoom/m_surface_ZB0p1_SNR20_R20_skills_v07D7p1_LF_321_kidsPhotometry_gold.csv'
mc_surface = pd.read_csv(inpath_doom)

# +++++++++++++++++++++++++++ correct the m 

ratio_blend = mc_surface['WeiBlend'].values / mc_surface['WeiTot'].values
for col_m in ['m1', 'm2']:

    mc_surface.loc[:, f'{col_m}_final'] = mc_surface[f'{col_m}_raw'].values \
                                        + mc_surface[f'var_d{col_m}'].values * ratio_blend\
                                        + mc_surface[f'PSFmodelling_d{col_m}'].values

    mc_surface.loc[:, f'{col_m}_final_err'] = np.sqrt(
                                                np.square(mc_surface[f'{col_m}_raw_err'].values) \
                                                + np.square(mc_surface[f'var_d{col_m}_err'].values * ratio_blend)\
                                                + np.square(mc_surface[f'PSFmodelling_d{col_m}_err'].values)\
                                                )
# save back
# print('>>> sorted raw/final (1)', np.sort(mc_surface['m1_raw'].values/mc_surface['m1_final'].values))
# print('>>> sorted raw/final (2)', np.sort(mc_surface['m2_raw'].values/mc_surface['m2_final'].values))
mc_surface.to_csv(inpath_doom, index=False, float_format='%.6f')
print(f'results saved back to {inpath_doom}')
