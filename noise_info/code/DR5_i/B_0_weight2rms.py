# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-11-01 14:01:02
# @Last Modified by:   lshuns
# @Last Modified time: 2023-01-19 15:47:44

### calculate rms from AW weight map
###### the resulted rms is in flux value with zeropoint=30
######      only central region is used (10,000 x 10,000)
######      mask is NOT applied
######  THELI = AW * (0.214/0.2)**2.
######  individual exposure = THELI * (5**0.5) for g r i

import os
import sys

import numpy as np
import pandas as pd

from astropy.io import fits

# +++++++++++++++++++++++ I/O

# zeropoint
dmag_info = pd.read_csv('../../kids_dr5_dmag_iS.csv')

# where to find the weight images
inDir_wei = '/disks/shear15/KiDS_DR5/Awe/weight_images'

# output
outpath = '../../test_noise_dr5_iS.csv'

# +++++++++++++++++++++++ workhorse

# extract AW rms 
df = pd.DataFrame({'label': ['00']*len(dmag_info),
                    'rmsAW_i1': np.zeros(len(dmag_info)),
                    'rmsAW_i2': np.zeros(len(dmag_info))})
for i_tile, row in dmag_info.iterrows():

    label = row['label']
    df.loc[i_tile, 'label'] = label
    print('>>>', label)

    for i_label in ('i1', 'i2'):
        zeropoint = row[f'dmag_{i_label}']

        # weight images
        infile = os.path.join(inDir_wei, f'KIDS_{label}_{i_label}.weight.fits.gz')
        with fits.open(infile) as hdul:
            data = hdul[0].data
        ## central region
        X_median = int(len(data[0, :])/2.)
        Y_median = int(len(data[:, 0])/2.)
        data = data[(Y_median-5000):(Y_median+5000), (X_median-5000):(X_median+5000)]

        # mask out bad pixels
        data = data[(data!=0)]

        # weight to rms
        rms = 1./np.sqrt(np.nanmedian(data))
        del data
        rms *= 10**(12-0.4*zeropoint) # to dmag=30

        # save
        df.loc[i_tile, f'rmsAW_{i_label}'] = rms
        print('    +++', i_label, rms)

# rms modification for different scenario
## the THELI rms
for i_label in ('i1', 'i2'): 
    df.loc[:, f'rmsTHELI_{i_label}'] = df[f'rmsAW_{i_label}'].values * (0.214/0.2)**2.
    df.loc[:, f'rmsExpo_{i_label}'] = df[f'rmsTHELI_{i_label}'].values * (5**0.5)

# save
df.to_csv(outpath, index=False)
print(f'combined results saved to {outpath}')

# combined results saved to ../../test_noise_dr5_iS.csv
# Elapsed:29:54:20.70,User=92953.297,System=14431.363,CPU=99.7%.