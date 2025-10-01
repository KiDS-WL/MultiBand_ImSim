# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-11-19 16:55:45
# @Last Modified by:   lshuns
# @Last Modified time: 2021-11-19 17:20:04

### build noise files for sensitivity test

import numpy as np
import pandas as pd


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 1. PSF ellipticity

# I/O
## labels used
infile4labels = '../skills_fiducial/noise_selec_135_0_N18_part0.csv'
## psf info
infile4psf = '../kids_dr4_psf_moffat_fromcoeffs.csv'
outfile = '../skills_sensitivity/kids_dr4_psf_moffat_fromcoeffs_part0.csv'

# desired psf factors
f_psf_e_list = [0.4, 0.6, 0.8, 1.2, 1.4, 1.6]

# 1. select target labels
labels = pd.read_csv(infile4labels)['label'].values
cata = pd.read_csv(infile4psf)
cata = cata[np.isin(cata['label'].values, labels)]
cata.reset_index(drop=True, inplace=True)

# 2. get modified psf e
for f_psf_e in f_psf_e_list:
    cata.loc[:, 'seeing_e1_'+str(f_psf_e).replace('.', 'p')+'_r'] = f_psf_e * cata['seeing_e1_r'].values
    cata.loc[:, 'seeing_e2_'+str(f_psf_e).replace('.', 'p')+'_r'] = f_psf_e * cata['seeing_e2_r'].values

# 3. save
cata.to_csv(outfile, index=False)
print('saved to', outfile)