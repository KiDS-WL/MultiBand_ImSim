# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-09-13 14:39:40
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-17 14:39:54

### delete images

import pandas as pd
import numpy as np 
import shutil
import pathlib

import re
import os
import glob

# where saved the selected tiles
# test_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7_nU10/'
# test_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7_nU20/'
# test_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7_qU5/'
# test_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7_qU10/'
# test_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7_sizeU5/'
# test_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7_sizeU10/'

parts = ['part0']
shear_labels = ['m283m283',  'm283p283',  'p283m283',  'p283p283']

# loop over and delete
for part in parts:
    for shear_label in shear_labels:

        print('for ', part, shear_label)

        subdir = os.path.join(test_dir, part, shear_label, 'images')
        shutil.rmtree(subdir)