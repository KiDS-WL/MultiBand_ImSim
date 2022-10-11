# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-03-23 13:39:24
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-10 14:18:46

### a combined catalogue with everthing from ImSim pipeline

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ general info

# skills version
skills_version = 'skills_v07D7'

# lensfit version
LF_version = '321'

# combined catalogue info
# combined_suffix = 'combined'
combined_suffix = 'combined_kids_filters'

# output file type
## fits, feather or csv
out_type = 'feather'

# main directory contains all the outputs
main_dir = f'/disks/shear16/ssli/ImSim/output/{skills_version}'

# does main dir include sub-dir from different input
split_ran = True
split_ran_tags = ['part0', 'part1', 'part2', 'part3', 'part4', 'part5']

# shear and rotation info
unique_shear_tags = ['m283m283', 'm283p283', 'p283m283', 'p283p283']
unique_rots = [0., 90.]

# out info
# outpath = os.path.join(main_dir, f'{skills_version}_LF_{LF_version}_everything.feather')
outpath = os.path.join(main_dir, f'{skills_version}_LF_{LF_version}_kidsPhotometry_everything.feather')

# +++++++++++++++++++++++++++++ workhorse

## subdir is split run
if split_ran:
    subdirs = []
    for part_tag in split_ran_tags:
        for shear_tag in unique_shear_tags:
            subdirs.append(os.path.join(main_dir, part_tag, shear_tag))

## subdir is target
else:
    subdirs = [os.path.join(main_dir, x) for x in unique_shear_tags]

print('Total number of subdirs:', len(subdirs))

## check if directories exist
for subdir in subdirs:
    if not os.path.isdir(subdir):
        raise Exception(f'{subdir} does not exist!')

## loop over shear tags
cata_final = [] 
for run_tag in unique_shear_tags:
    print('Running for', run_tag)

    # select dir belong to this shear
    subdirs_selec = [subdir for subdir in subdirs if run_tag in subdir]
    print('number of subdirs:', len(subdirs_selec))

    # ++++++++++ 1. loop over all outputs
    for subdir in subdirs_selec:

        print('working on subdir', subdir)

        # basic info for input shear values
        with open(os.path.join(subdir, 'basic_info.txt'), 'r') as opened_file:
            all_lines = opened_file.readlines()
        useful_line = [line for line in all_lines if 'g_cosmic' in line][0]
        try:
            g1 = float(useful_line.split()[2])
        except ValueError:
            g1 = float(useful_line.split()[2][:-1])
        g2 = float(useful_line.split()[3])

        # main catalogues from simulation
        file_list = glob.glob(os.path.join(subdir, 'catalogues', f'*_{combined_suffix}.*'))
        if not file_list:
            raise Exception(f'Cannot find any combined catalogues in tag {run_tag}!\n\
            make sure taskID=7 is performed in the main pipeline!')
        print(f'Number of files: {len(file_list)}' )
        for file in file_list:
            file_name = os.path.basename(file)
            file_type = file[-3:]
            if file_type == 'csv':
                cata = pd.read_csv(file)
            elif file_type == 'her':
                cata = pd.read_feather(file)
            elif file_type == 'its':
                with fits.open(file) as hdul:
                    cata = Table(hdul[1].data).to_pandas()
            else:
                raise Exception(f'Not supported input file type! {file}')

            # some columns for tile info
            ## input cosmic shear
            cata.loc[:, 'g1_in'] = g1
            cata.loc[:, 'g2_in'] = g2
            ## run tag
            cata.loc[:, 'run_tag'] = run_tag        
            ## tile noise info
            cata.loc[:, 'tile_label'] = re.search(r'tile(.*)_rot', file_name).group(1)
            ## rotation info
            cata.loc[:, 'gal_rot'] = float(re.search(r'_rot(\d+)', file_name).group(1))

            # collect
            cata_final.append(cata)
            del cata

cata_final = pd.concat(cata_final, ignore_index=True)
print(f'Total number of source {len(cata_final)}')
cata_final.to_feather(outpath)
print('combined cata saved to', outpath)

########### original Shark photometry
# Total number of subdirs: 24
# Running for m283m283
# number of subdirs: 6
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/m283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part1/m283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/m283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/m283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/m283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/m283m283
# Number of files: 36
# Running for m283p283
# number of subdirs: 6
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/m283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part1/m283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/m283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/m283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/m283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/m283p283
# Number of files: 36
# Running for p283m283
# number of subdirs: 6
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/p283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part1/p283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/p283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/p283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/p283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/p283m283
# Number of files: 36
# Running for p283p283
# number of subdirs: 6
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/p283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part1/p283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/p283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/p283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/p283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/p283p283
# Number of files: 36
# Total number of source 89303975
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_everything.feather
# Elapsed:14:00.20,User=1444.176,System=4452.990,CPU=701.8%.

# ####### kids photometry
# Total number of subdirs: 24
# Running for m283m283
# number of subdirs: 6
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/m283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part1/m283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/m283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/m283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/m283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/m283m283
# Number of files: 36
# Running for m283p283
# number of subdirs: 6
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/m283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part1/m283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/m283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/m283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/m283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/m283p283
# Number of files: 36
# Running for p283m283
# number of subdirs: 6
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/p283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part1/p283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/p283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/p283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/p283m283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/p283m283
# Number of files: 36
# Running for p283p283
# number of subdirs: 6
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part0/p283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part1/p283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part2/p283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part3/p283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part4/p283p283
# Number of files: 36
# working on subdir /disks/shear16/ssli/ImSim/output/skills_v07D7/part5/p283p283
# Number of files: 36
# Total number of source 89303975
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_everything.feather
# Elapsed:13:05.10,User=1466.580,System=3862.607,CPU=678.7%.
