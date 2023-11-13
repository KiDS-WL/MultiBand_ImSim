# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-01-26 13:54:14
# @Last Modified by:   lshuns
# @Last Modified time: 2023-02-04 14:39:38

### combine the i2 input with other bands

import os 
import glob
import pathlib

import pandas as pd 
import numpy as np

### general info

main_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7'

# parts = ['part0', 'part1', 'part2', 'part3', 'part4', 'part5']
# parts = ['part0', 'part2', 'part3', 'part4', 'part5']
parts = ['part1']
shear_tags = ['m283p283', 'm283m283', 'p283p283', 'p283m283']

# where to find original catalogues
folder0 = 'input'

# where to find new catalogues
folder1 = 'input_dr5_iS'
## used columns
cols1 = ['i2_input']

# where to save
folder2 = 'input_dr5_ten'

# check column
col_common = 'index_input'

### loop over and copy
for part in parts:
    for shear_tag in shear_tags:

        # original 
        inpath0 = os.path.join(main_dir, part, shear_tag, 'catalogues', folder0)
        print('original catalogues from', inpath0)
        ## find all catalogues
        inpath0_files = glob.glob(os.path.join(inpath0, '*.feather'))
        print('    number of files found', len(inpath0_files))

        # new
        inpath1 = os.path.join(main_dir, part, shear_tag, 'catalogues', folder1)
        print('new catalogues from', inpath1)

        # for save
        inpath2 = os.path.join(main_dir, part, shear_tag, 'catalogues', folder2)
        print('save to', inpath2)
        pathlib.Path(inpath2).mkdir(parents=True, exist_ok=True)

        # loop over files
        for inpath0_file in inpath0_files:
            # the original cata
            cata0 = pd.read_feather(inpath0_file)
            # the new cata
            cata1 = pd.read_feather(os.path.join(inpath1, os.path.basename(inpath0_file)))

            # check if consistent
            if not np.array_equal(cata0[col_common].values, cata1[col_common].values):
                raise Exception(inpath0_file, 'inconsistent!')

            # combine
            cata0 = pd.concat([cata0, cata1[cols1]], axis=1)
            del cata1

            # save
            outpath = os.path.join(inpath2, os.path.basename(inpath0_file))
            cata0.to_feather(outpath)
            # print('saved to', outpath)
            del cata0

# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/m283p283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/m283p283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/m283p283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/m283m283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/m283m283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/m283m283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/p283p283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/p283p283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/p283p283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/p283m283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/p283m283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/p283m283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/m283p283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/m283p283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/m283p283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/m283m283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/m283m283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/m283m283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/p283p283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/p283p283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/p283p283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/p283m283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/p283m283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/p283m283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/m283p283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/m283p283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/m283p283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/m283m283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/m283m283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/m283m283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/p283p283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/p283p283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/p283p283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/p283m283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/p283m283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/p283m283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/m283p283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/m283p283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/m283p283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/m283m283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/m283m283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/m283m283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/p283p283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/p283p283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/p283p283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/p283m283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/p283m283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/p283m283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/m283p283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/m283p283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/m283p283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/m283m283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/m283m283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/m283m283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/p283p283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/p283p283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/p283p283/catalogues/input_dr5_ten
# original catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/p283m283/catalogues/input
#     number of files found 36
# new catalogues from /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/p283m283/catalogues/input_dr5_iS
# save to /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/p283m283/catalogues/input_dr5_ten
# Elapsed:3:31.36,User=291.831,System=199.835,CPU=232.6%.
