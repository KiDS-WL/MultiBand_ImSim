# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-05-13 10:17:47
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-10 14:25:40

### apply the shear catalogue selection

import os
import re

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ general info

# input 
inpath = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7ten_LF_321_kidsPhotometry_everything_col_flag_dr5.feather'
outpath = inpath.replace('_everything_col_flag_dr5.feather', '_shear_noSG_noWeiCut_newCut_dr5.feather')

# +++++++++++++++++++++++++++++ workhorse

# load the catalogue
cata = pd.read_feather(inpath)
print('number in ori', len(cata))

# apply selection
cata = cata[(cata['flag_gaap']==0)&(cata['flag_asteroid']==0)&(cata['flag_binary']==0)&(cata['flag_LF_noWeiCut']==0)]
cata.drop(columns=['flag_gaap', 'flag_asteroid', 'flag_binary', 'flag_LF_noWeiCut'], inplace=True)
cata.reset_index(drop=True, inplace=True)
print('number after selection', len(cata))

# save
cata.to_feather(outpath)
print('saved to', outpath)

# number in ori 89303769
# number after selection 47995929
# saved to /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7ten_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_dr5.feather
# Elapsed:7:23.78,User=811.429,System=3828.955,CPU=1045.6%.