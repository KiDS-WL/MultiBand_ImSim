# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-08-29 09:20:26
# @Last Modified by:   lshuns
# @Last Modified time: 2022-08-29 09:25:14

### delete copied images

import pandas as pd
import numpy as np 
import shutil
import pathlib

import re
import os
import glob

# where saved the selected tiles
test_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7_PSFmodelling/'
shear_labels = ['m283m283',  'm283p283',  'p283m283',  'p283p283']

# loop over shear tags to delete

for shear_label in shear_labels:

    print('for ', shear_label)

    subdir = os.path.join(test_dir, shear_label, 'images', 'original')

    # delete images
    dir_tmp_list = glob.glob(os.path.join(subdir, f'chips_*'))
    print('number of images found', len(dir_tmp_list))
    for dir_tmp in dir_tmp_list:
        shutil.rmtree(dir_tmp)