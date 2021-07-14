# @Author: lshuns
# @Date:   2021-03-25, 17:51:54
# @Last modified by:   ssli
# @Last modified time: 2021-03-26, 15:27:37

### collect PSF info for ugri
###     File seeing & e from esodr4_stats.dat
###     Measured seeing & Moffat beta parameter from AW images by modified GAaP code

import os

import numpy as np
import pandas as pd

# ++++++++++++++++++ I/O

# noise file
noise_file = '../test_noise_ugri.csv'

# directory to the original KiDS noise files
### esodr4_stats.dat, moffat.{band}
inDir = '/disks/shear10/ssli/KiDS_noise_seeing/'

# ++++++++++++++++++ workhorse
bands = ['u', 'g', 'r', 'i']

# main noise file
data_final = pd.read_csv(noise_file)

# DR4 state info
file = os.path.join(inDir, 'esodr4_stats.dat')
data_psf = pd.read_csv(file, names=['label', 'filter', 'seeing', 'e', 'mag1', 'mag2', 'mag3'])
print(f'File psf loaded from {file}')

# get the reported psf & ellipticity
for (index, row) in data_final.iterrows():

    tile_label = row['label']

    for band in bands:

        # seeing
        data_final.loc[data_final['label']==tile_label, f'FileSeeing_{band}'] = data_psf.loc[(data_psf['label']==f'KIDS_{tile_label}') & (data_psf['filter']==f'OCAM_{band}_SDSS'), 'seeing'].values[0]

        # ellipticity
        data_final.loc[data_final['label']==tile_label, f'seeing_e_{band}'] = data_psf.loc[(data_psf['label']==f'KIDS_{tile_label}') & (data_psf['filter']==f'OCAM_{band}_SDSS'), 'e'].values[0]

# collect measured seeing & Moffat beta parameter
for band in bands:

    file = os.path.join(inDir, f'moffat.{band}')
    data = pd.read_csv(file, sep=" ", header=None)
    print(f'Measured psf loaded from {file}')
    data.columns = ["label", 'no', f"MeasSeeing_{band}", f"MeasBeta_{band}"]
    data.drop(columns=['no'], inplace=True)

    data.loc[:, 'label'] = data['label'].str.extract(r'KIDS_(.*)_[ugri]').values

    data[f"MeasSeeing_{band}"] *= 0.2 # pixel to arcsec

    print(data)
    print(data_final)

    ## merge
    data_final = data_final.merge(data, on='label')

# save back to the original file
data_final.to_csv(noise_file, index=False)
print(f'PSF info saved back to {noise_file}')
