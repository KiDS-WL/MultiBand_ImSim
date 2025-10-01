# @Author: ssli
# @Date:   2021-05-11, 15:31:50
# @Last modified by:   ssli
# @Last modified time: 2021-05-17, 22:07:04

### combine noise files for ImSim

import numpy as np
import pandas as pd

# ++++++++++ I/O

# original noise files
noise_files = ['../test_noise_ugri.csv', '../test_noise_NIR.csv']

# psf e info from shear catalogue
psf_e_file = '../test_psf_e_fromShearCata.csv'

# output
outfile = '../kids_dr4_noise_ugriZYJHKs.csv'

# +++++++++++ workhorse

# 1. combine noise info
noise_info = None
for noise_file in noise_files:
    tmp = pd.read_csv(noise_file)
    if noise_info is None:
        noise_info = tmp
    else:
        noise_info = noise_info.merge(tmp, on='label')
del tmp

# 2. dummy columns for ImSim
for band in ['Z', 'Y', 'J', 'H', 'Ks']:
    noise_info.loc[:, f'rmsExpo_{band}'] = np.array(noise_info[f'rms_{band}'])

# 3. psf e for each exposure from shear catalogue
psf_e = pd.read_csv(psf_e_file)
## fill nan with mean value
for index, row in psf_e[psf_e.isnull().any(axis=1)].iterrows():
    for i in range(5):
        if pd.isna(row[f'PSF_e_mean_expo{i}']):
            print('nan in', index, f'PSF_e_mean_expo{i}')
            psf_e.loc[index, f'PSF_e_mean_expo{i}'] = psf_e.loc[index, f'PSF_e_mean']
## rename for better cross match
psf_e.loc[:, 'label'] = psf_e['label'].str.replace('p', '.').values
psf_e.loc[:, 'label'] = psf_e['label'].str.replace('m', '-').values
psf_e.loc[:, 'label'] = psf_e['label'].str.extract(r'KIDS_(.*)').values
## get desired columns
psf_e = psf_e[['label']+[f'PSF_e_mean_expo{i}' for i in range(5)]]
## change (1-q)/(1+q) to 1-q
for i in range(5):
    tmp = np.array(psf_e[f'PSF_e_mean_expo{i}'])
    psf_e.loc[:, f'PSF_e_mean_expo{i}'] = 1 - (1-tmp)/(1+tmp)
## cross match
noise_info = noise_info.merge(psf_e, on='label')

# 4. rename columns for ImSim
rename_col = {}
old_names = [f'PSF_e_mean_expo{i}' for i in range(5)]
new_names = [f'seeing_e1_r_expo{i}' for i in range(5)]
for i_name, old_name in enumerate(old_names):
    rename_col[old_name] = new_names[i_name]
noise_info.rename(columns=rename_col, inplace=True)

# 5. save
noise_info.to_csv(outfile, index=False)
print(f'combined noise info saved as {outfile}')
