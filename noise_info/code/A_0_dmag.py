# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-10-03 09:53:04
# @Last modified by:   lshuns
# @Last modified time: 2021-04-26, 20:15:59

### get the dmag for u,g,r,i

import os
import sys
import logging
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from astropy.io import fits
import multiprocessing as mp

logging.basicConfig(format='%(name)s : %(levelname)s - %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)

# +++++++++++++++++++ I/O

inDir = '/disks/shear15/KiDS/KiDS-1000/ESO-DR4-photometry-catalogues/'

# +++++++++++++++++++ workhorse
bands = ['u', 'g', 'r', 'i']
label_list = []
zero_u = []
zero_g = []
zero_r = []
zero_i = []

# +++ main
pathlist = Path(inDir).rglob('KiDS_DR4.0_*_ugriZYJHKs_cat.fits')
for path in pathlist:
    # because path is object not string
    path_in_str = str(path)
    label = path_in_str[71:-20]
    label_list.append(label)

    with fits.open(path_in_str) as hdul:
        zero = hdul[0].header

    logger.info(path_in_str)

    for band in bands:
        zeropoint = zero['DMAG_{:}'.format(band.upper())]
        if band == 'u':
            zero_u.append(zeropoint)
        elif band == 'g':
            zero_g.append(zeropoint)
        elif band == 'r':
            zero_r.append(zeropoint)
        elif band == 'i':
            zero_i.append(zeropoint)

pathlist = Path(inDir).rglob('KiDS_DR4.1_*_ugriZYJHKs_cat.fits')
for path in pathlist:
    # because path is object not string
    path_in_str = str(path)
    label = path_in_str[71:-20]
    label_list.append(label)

    with fits.open(path_in_str) as hdul:
        zero = hdul[0].header

    logger.info(path_in_str)

    for band in bands:
        zeropoint = zero['DMAG_{:}'.format(band.upper())]
        if band == 'u':
            zero_u.append(zeropoint)
        elif band == 'g':
            zero_g.append(zeropoint)
        elif band == 'r':
            zero_r.append(zeropoint)
        elif band == 'i':
            zero_i.append(zeropoint)

final = pd.DataFrame({'label': label_list,
                        'dmag_u':zero_u,
                        'dmag_g':zero_g,
                        'dmag_r':zero_r,
                        'dmag_i':zero_i
                        })
final.to_csv('../kids_dr4_dmag_ugri.csv', index=False)
