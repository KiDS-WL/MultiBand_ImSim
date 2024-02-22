# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-07-08 10:43:44
# @Last Modified by:   lshuns
# @Last Modified time: 2023-02-28 08:59:40

### copy selected images and catalogues from the fiducial runs

import pandas as pd
import numpy as np 
import shutil
import pathlib

import re
import os
import glob

# where to find the fiducial simulations
main_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7/'
## folder info
parts = ['part0', 'part1', 'part2', 'part3', 'part4', 'part5']
shear_labels = ['m283m283',  'm283p283',  'p283m283',  'p283p283']

# where to save the selected tiles
test_dir = '/disks/shear16/ssli/ImSim/output/skills_v07D7p1_PSFmodelling/'

# selected tiles
inpath_selected = './tile_forTEST.csv'
tile_selec = pd.read_csv(inpath_selected)['label'].values
print('selected tiles', tile_selec)

# loop over tiles and find info
for tile in tile_selec:

    print('>>>>', tile)

    # loop over parts to find tiles
    for part in parts:
        main_dir_part = os.path.join(main_dir, part)
        ## get all tile info in this part
        files_tmp = glob.glob(os.path.join(main_dir_part, shear_labels[0], 'catalogues', 'input', 'noise_info_tile*.csv'))
        tile_labels_tmp = [re.search(r'tile(.*).csv', file).group(1) for file in files_tmp]

        if (tile in tile_labels_tmp):
            print('finding in ', part)

            # copy all shear labels
            for shear_label in shear_labels:

                print('for ', shear_label)

                subdir = os.path.join(main_dir_part, shear_label)

                # the folders for test
                subdir_test = os.path.join(test_dir, shear_label)
                pathlib.Path(subdir_test).mkdir(parents=True, exist_ok=True)
                subdir_test_cata = os.path.join(subdir_test, 'catalogues')
                pathlib.Path(subdir_test_cata).mkdir(parents=True, exist_ok=True)
                subdir_test_ima = os.path.join(subdir_test, 'images')
                pathlib.Path(subdir_test_ima).mkdir(parents=True, exist_ok=True)

                # get the basic info
                if not os.path.exists(os.path.join(subdir_test, 'basic_info.txt')):
                    shutil.copy(os.path.join(subdir, 'basic_info.txt'), subdir_test)

                # get catalogues
                ## input 
                subdir_test_cata_tmp = os.path.join(subdir_test_cata, 'input')
                pathlib.Path(subdir_test_cata_tmp).mkdir(parents=True, exist_ok=True)
                files_tmp = glob.glob(os.path.join(subdir, 'catalogues', 'input', f'*tile{tile}*'))
                print('number of inputs found', len(files_tmp))
                for file in files_tmp:
                    if not os.path.exists(os.path.join(subdir_test_cata_tmp, os.path.basename(file))):
                        shutil.copy(file, subdir_test_cata_tmp)
                ## CrossMatch
                subdir_test_cata_tmp = os.path.join(subdir_test_cata, 'CrossMatch')
                pathlib.Path(subdir_test_cata_tmp).mkdir(parents=True, exist_ok=True)
                files_tmp = glob.glob(os.path.join(subdir, 'catalogues', 'CrossMatch', f'tile{tile}_*'))
                print('number of CrossMatch found', len(files_tmp))
                for file in files_tmp:
                    if not os.path.exists(os.path.join(subdir_test_cata_tmp, os.path.basename(file))):
                        shutil.copy(file, subdir_test_cata_tmp)
                ## SExtractor
                subdir_test_cata_tmp = os.path.join(subdir_test_cata, 'SExtractor')
                pathlib.Path(subdir_test_cata_tmp).mkdir(parents=True, exist_ok=True)
                files_tmp = glob.glob(os.path.join(subdir, 'catalogues', 'SExtractor', f'tile{tile}_*'))
                print('number of SExtractor found', len(files_tmp))
                for file in files_tmp:
                    if not os.path.exists(os.path.join(subdir_test_cata_tmp, os.path.basename(file))):
                        shutil.copy(file, subdir_test_cata_tmp)

                # get images
                subdir_test_ima_tmp = os.path.join(subdir_test_ima, 'original')
                pathlib.Path(subdir_test_ima_tmp).mkdir(parents=True, exist_ok=True)
                dir_tmp_list = glob.glob(os.path.join(subdir, 'images', 'original', f'chips_tile{tile}_*'))
                print('number of images found', len(dir_tmp_list))
                for dir_tmp in dir_tmp_list:
                    if not os.path.exists(os.path.join(subdir_test_ima_tmp, os.path.basename(dir_tmp))):
                        shutil.copytree(dir_tmp, os.path.join(subdir_test_ima_tmp, os.path.basename(dir_tmp)))

# Elapsed:4:08:07.70,User=192.058,System=1515.111,CPU=11.4%.