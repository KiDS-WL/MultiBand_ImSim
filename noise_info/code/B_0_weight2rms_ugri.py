# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-10-03 09:50:05
# @Last modified by:   lshuns
# @Last modified time: 2021-03-25, 20:28:54

### calculate rms from AW weight map with parallel
###### the resulted rms is in flux value with zeropoint=30
######      only central region is used (10,000 x 10,000)
######      mask is applied (==0)
######  THELI = AW * (0.214/0.2)**2.
######  individual exposure = THELI * (5**0.5) for g r i, * (4**0.5) for u

import os
import sys
import logging
import subprocess

import numpy as np
import pandas as pd
import multiprocessing as mp

from astropy.io import fits

logging.basicConfig(format='%(name)s : %(levelname)s - %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)

# +++++++++++++++++++++++ I/O

# zeropoint
dmag_ugri = pd.read_csv('../kids_dr4_dmag_ugri.csv')

# weight images
inDir_wei = '/disks/shear15/KiDS/KiDS-1000/Awe/weight_images/'

# mask images
inDir_msk = '/disks/shear15/KiDS/KiDS-1000/Awe/mask_images/'

# output
outDir = '../'

# +++++++++++++++++++++++ workhorse
# the function
def rmsFunc(indir_wei, indir_msk, out_dir, band, dmag_info):

    label_list = []
    rms_list = []
    for i_tile, row in dmag_info.iterrows():

        label = row['label']
        zeropoint = row[f'dmag_{band}']

        # weight images
        infile = indir_wei + f'KIDS_{label}_{band}_wei.fits'
        with fits.open(infile) as hdul:
            data = hdul[0].data
        ## central region
        X_median = int(len(data[0, :])/2.)
        Y_median = int(len(data[:, 0])/2.)
        data = data[(Y_median-5000):(Y_median+5000), (X_median-5000):(X_median+5000)]

        # mask images
        infile = indir_msk + f'KIDS_{label}_{band}_msk.fits'
        with fits.open(infile) as hdul:
            data_msk = hdul[0].data
        ## central region
        X_median = int(len(data_msk[0, :])/2.)
        Y_median = int(len(data_msk[:, 0])/2.)
        data_msk = data_msk[(Y_median-5000):(Y_median+5000), (X_median-5000):(X_median+5000)]

        # mask out bad pixels
        mask_tmp = (data!=0) & (data_msk==0)

        # weight to rms
        rms = 1./np.sqrt(np.nanmedian(data[mask_tmp]))
        rms *= 10**(12-0.4*zeropoint) # to dmag=30
        label_list.append(label)
        rms_list.append(rms)
        logger.debug('{:} {:} {:}'.format(band, label, rms))

    df = pd.DataFrame({'label': label_list,
                        'rms':rms_list})
    outfile = out_dir + f'test_AW_rms_{band}.csv'
    df.to_csv(outfile, index=False)
    logger.info('rms info saved as {:}'.format(outfile))

# extract rms for each band
## parallel implemented
bands = ['u', 'g', 'r', 'i']
p_list = []
for band in bands:
    dmag = dmag_ugri[['label', f'dmag_{band}']]

    p = mp.Process(target=rmsFunc, args=(inDir_wei, inDir_msk, outDir, band, dmag))
    p.start()
    p_list.append(p)
for p in p_list:
    p.join()

# combine to one file
for band in bands:
    infile = outDir + f'test_AW_rms_{band}.csv'
    data = pd.read_csv(infile)
    data.rename(columns={'rms': f'rmsAW_{band}'}, inplace=True)

    if band == 'u':
        data_final = data
    else:
        data_final = data_final.merge(data, on='label')

# rms modification for different scenario
for band in bands:

    ## the THELI rms
    data_final.loc[:, f'rmsTHELI_{band}'] = data_final[f'rmsAW_{band}'].values * (0.214/0.2)**2.

    # the exposure rms
    if band == 'u':
        data_final.loc[:, f'rmsExpo_{band}'] = data_final[f'rmsTHELI_{band}'].values * (4**0.5)
    else:
        data_final.loc[:, f'rmsExpo_{band}'] = data_final[f'rmsTHELI_{band}'].values * (5**0.5)

# save
outpath = outDir + 'test_noise_ugri.csv'
data_final.to_csv(outpath, index=False)
logger.info(f'combined results saved to {outpath}')

# delete intermediate files
for band in bands:
    infile = outDir + f'test_AW_rms_{band}.csv'
    os.remove(infile)
    logger.info(f'{infile} removed.')
