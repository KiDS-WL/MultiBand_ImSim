# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-11-02 09:53:06
# @Last Modified by:   lshuns
# @Last Modified time: 2023-01-19 16:08:17

### save PSF and noise info together

import os

import numpy as np
import pandas as pd

# ++++++++++++++++++ I/O

# noise file (from B_0_weight2rms.py)
noise_file = '../../test_noise_dr5_iS.csv'

# the measured moffat profiles (from B_1_*.py and B_2_*.py)
band_list = ['i1', 'i2']
inpath_list = ['/disks/shear15/ssli/KiDS_noise_seeing/psf2moffat_i1/moffat_i1.csv', 
                '/disks/shear15/ssli/KiDS_noise_seeing/psf2moffat_i2/moffat_i2.csv']

# ++++++++++++++++++ workhorse

# main noise file
data_final = pd.read_csv(noise_file)
print('number in main file', len(data_final))
# print(data_final.columns)

for i_band, band in enumerate(band_list):

    print('>>> for band', band)

    ## clean previous psf info
    if f'MeasSeeing_{band}' in data_final.columns:
        data_final.drop(columns=[f'MeasSeeing_{band}', f'MeasBeta_{band}', f'seeing_e1_{band}'], inplace=True)

    ## collect measured seeing & Moffat beta parameter
    data = pd.read_csv(inpath_list[i_band])
    print('number in psf file', len(data))
    data.drop(columns=['no'], inplace=True)

    ## pixel to arcsec
    data[f"MeasSeeing_{band}"] *= 0.2

    ## merge
    data_final = data_final.merge(data, on='label')
    print('number after merge', len(data_final))
    del data

# print(data_final.columns)
# save back to the original file
data_final.to_csv(noise_file, index=False)
print(f'PSF info saved back to {noise_file}')

# number in main file 1347
# >>> for band i1
# number in psf file 1347
# number after merge 1347
# >>> for band i2
# number in psf file 1347
# number after merge 1347
# PSF info saved back to ../../test_noise_dr5_iS.csv
