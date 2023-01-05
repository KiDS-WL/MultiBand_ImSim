# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-07-08 10:27:24
# @Last Modified by:   lshuns
# @Last Modified time: 2022-08-26 11:06:46

### select tiles for the test
####### 7 tiles covering PSF size distributions

import pandas as pd 
import numpy as np 

# ++++++++++++++++ I/O

# catalogue with all noise info
inpath = '../../noise_info/skills_fiducial/noise_selec_combined.csv'
cata0 = pd.read_csv(inpath)
## used columns
cata0 = cata0[['label', 'InputSeeing_r']]

# where to save selected
outpath = './tile_forTEST.csv'

# random seed
rng_seed = 502
Nsample = 23

# ++++++++++++++++ workhorse

# selected quantile
quantiles = [0, 0.05, 0.16, 0.5, 0.84, 0.95, 1.0]
cata_selected = cata0[np.isin(cata0['InputSeeing_r'], cata0['InputSeeing_r'].quantile(quantiles, interpolation='lower'))]

# add more randomly selected tiles
### do not repeat those already selected
cata0 = cata0[~np.isin(cata0['InputSeeing_r'].values, cata_selected['InputSeeing_r'].values)]
cata0.reset_index(drop=True, inplace=True)
cata0 = cata0.sample(n=Nsample, random_state=rng_seed)

# combine
cata_selected = pd.concat([cata_selected, cata0], ignore_index=True)
del cata0
print(cata_selected)

# save 
cata_selected.to_csv(outpath, index=False)
print('saved to ', outpath)

#           label  InputSeeing_r
# 0     136.0_1.5       0.738957
# 1    180.0_-1.5       0.840314
# 2     181.0_1.5       0.538631
# 3    229.0_-1.5       0.396071
# 4     8.2_-32.1       0.938986
# 5    43.5_-32.1       0.622086
# 6   333.8_-29.2       0.476176
# 7     231.0_1.5       0.649746
# 8    12.7_-30.2       0.464206
# 9     134.0_1.5       0.770465
# 10    9.5_-33.1       0.570216
# 11   45.9_-32.1       0.570968
# 12    180.5_2.5       0.554745
# 13    136.0_0.5       0.608580
# 14    229.0_0.5       0.722255
# 15   45.4_-31.2       0.574215
# 16   178.0_-1.5       0.509327
# 17   10.6_-32.1       0.493088
# 18    133.0_1.5       0.595648
# 19   135.0_-1.5       0.581178
# 20   45.6_-29.2       0.622978
# 21   137.0_-0.5       0.590776
# 22   232.0_-0.5       0.672112
# 23   178.0_-0.5       0.486061
# 24   46.3_-33.1       0.701185
# 25  332.9_-32.1       0.650636
# 26  337.0_-30.2       0.911850
# 27   179.0_-1.5       0.638005
# 28    179.0_1.5       0.561644
# 29   137.0_-1.5       0.587224
# saved to  ./tile_forTEST.csv