# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-10-27 10:39:10
# @Last Modified by:   lshuns
# @Last Modified time: 2023-01-17 15:55:52

### get the dmag for i2

import pandas as pd 
import numpy as np 

# the RD5 zero points
inpath = '/net/grecht/data2/ssli_files/chrome_downloads/dmags.txt'
cata0 = pd.read_csv(inpath,
                 delim_whitespace=True, #separator is whitespace
                 header=None, #no header
                 usecols=[0,1,2],
                 names=['label','band','dmag']) #set columns names

# get the i1 band values
cata_i1 = cata0[cata0['band']=='i1']
cata_i1 = cata_i1.drop(columns=['band']).rename(columns={'dmag': 'dmag_i1'})
print('number of tiles (i1)', len(cata_i1))

# get the i2 band values
cata_i2 = cata0[cata0['band']=='i2']
del cata0
cata_i2 = cata_i2.drop(columns=['band']).rename(columns={'dmag': 'dmag_i2'})
print('number of tiles (i2)', len(cata_i2))

# merge i1 and i2
cata = cata_i1.merge(cata_i2, on='label')
print('number in merge', len(cata))

# rename the labels
cata.loc[:, 'label'] = [label[5:].replace('p', '.').replace('m', '-') for label in cata['label'].values]

# save
cata.to_csv('../../kids_dr5_dmag_iS.csv', index=False)

### out info
# number of tiles (i1) 1347
# number of tiles (i2) 1347
# number in merge 1347
