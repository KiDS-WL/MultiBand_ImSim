# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-01-19 17:43:13
# @Last Modified by:   lshuns
# @Last Modified time: 2023-01-19 17:48:42

### remove tiles with zero weight map

import pandas as pd 

# load noise 
noise_file = '../../test_noise_dr5_iS.csv'
cata = pd.read_csv(noise_file)
print('number original', len(cata))

# drop zero weight
cata = cata[(cata['rmsAW_i1']<1e3)&(cata['rmsAW_i2']<1e3)]
print('number after selection', len(cata))

# save back to the original file
cata.to_csv(noise_file, index=False)
print(f'PSF info saved back to {noise_file}')

# number original 1347
# number after selection 1346
# PSF info saved back to ../../test_noise_dr5_iS.csv
