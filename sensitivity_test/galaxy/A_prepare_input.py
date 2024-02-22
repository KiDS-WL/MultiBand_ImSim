# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-03-12 17:06:30
# @Last Modified by:   lshuns
# @Last Modified time: 2023-03-13 14:57:38

### prepare input catalogue by shifting galaxy properties
###### 1. learn err/val vs r_mag_auto relation from input catalogue (linear interpolation)
###### 2. shift the input para based on r_SDSS_apparent_corr

import pandas as pd 
import numpy as np 

from scipy import interpolate

########## COSMOS data catalogue
inpath = '/disks/shear10/ssli/ImSim/input/COSMOS_cata/cosmos_shape_z_uBVriZYJHKs.feather'
cosmos_cata = pd.read_feather(inpath)
print('COSMOS ori', len(cosmos_cata))
### select
###### 0. discard too small or too big galaxies
mask_re = (cosmos_cata['RE_GALFIT_HI']>=1e-2) & (cosmos_cata['RE_GALFIT_HI']<=10.)
###### 1. good shape
mask_galfit = (cosmos_cata['FLAG_GALFIT_HI']==0)
###### 2. has magnitude
mask_mag = (cosmos_cata['r_mag_auto']>0)
### apply
cosmos_cata = cosmos_cata[mask_galfit & mask_re & mask_mag]
del mask_galfit, mask_re
cosmos_cata.reset_index(drop=True, inplace=True)
print('COSMOS selected', len(cosmos_cata))

########## original SKiLLS input catalogue
inpath = '/disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0.feather'
cata0 = pd.read_feather(inpath)

# saved name
outpath = inpath.replace('.feather', '_shifted.feather')

########## morphology parameters
save_names = ['size', 'n', 'q']
limits = [[1e-2, 10.], [0.5, 6.0], [0.05, 1.0]]
skills_cols = ['Re_arcsec', 'shape/sersic_n', 'BA']
cosmos_cols = ['RE_GALFIT_HI', 'N_GALFIT_HI', 'BA_GALFIT_HI']
err_cols = ['REERR_GALFIT_HI', 'NERR_GALFIT_HI', 'BAERR_GALFIT_HI']

########## the xval
cosmos_xcol = 'r_mag_auto'
skills_xcol = 'r_SDSS_apparent_corr'

########## selected used columns
cosmos_cata = cosmos_cata[cosmos_cols+err_cols+[cosmos_xcol]]

########## start learning and changing
for i_col, skills_col in enumerate(skills_cols):

    save_name = save_names[i_col]
    limit = limits[i_col]
    cosmos_col = cosmos_cols[i_col]
    err_col = err_cols[i_col]

    ## >>>> 1. learn err/val vs r_mag_auto relation from input catalogue
    test_df = pd.DataFrame({'xval': cosmos_cata[cosmos_xcol].values,
                            'yval': cosmos_cata[err_col].values/cosmos_cata[cosmos_col].values})
    # group 
    test_df.loc[:, 'bin'] = pd.qcut(test_df['xval'].values, 10, labels=False)
    # get the median
    test_df_median = test_df.groupby('bin').median()
    del test_df

    ## >>>> 2. interpolate with linear
    ### extrapolate with constant
    ftmp = interpolate.interp1d(test_df_median['xval'].values, test_df_median['yval'].values,
                kind='linear', bounds_error=False,
                fill_value=(test_df_median['yval'].values[0], test_df_median['yval'].values[-1]))
    del test_df_median

    ## >>>> 3. change values in SKiLLS
    # get the shift amount
    y_learned = ftmp(cata0[skills_xcol])
    del ftmp
    ## check legitimate
    if np.min(y_learned) < 0:
        raise Exception('y_learned negative, something is wrong!')
    if np.max(y_learned) > 1:
        raise Exception('y_learned larger than 1, something is wrong!')
    cata0.loc[:, f'{skills_col}_err'] = y_learned
    # shift up
    cata0.loc[:, f'{skills_col}_U'] = (1+y_learned) * cata0[skills_col].values
    ### limits
    cata0.loc[:, f'{skills_col}_U'] = np.where(cata0[f'{skills_col}_U'].values<limit[0], 
                                                limit[0], cata0[f'{skills_col}_U'].values)
    cata0.loc[:, f'{skills_col}_U'] = np.where(cata0[f'{skills_col}_U'].values>limit[1], 
                                                limit[1], cata0[f'{skills_col}_U'].values)
    # shift down
    cata0.loc[:, f'{skills_col}_D'] = (1-y_learned) * cata0[skills_col].values
    ### limits
    cata0.loc[:, f'{skills_col}_D'] = np.where(cata0[f'{skills_col}_D'].values<limit[0], 
                                                limit[0], cata0[f'{skills_col}_D'].values)
    cata0.loc[:, f'{skills_col}_D'] = np.where(cata0[f'{skills_col}_D'].values>limit[1], 
                                                limit[1], cata0[f'{skills_col}_D'].values)
    del y_learned

## save results
cata0.to_feather(outpath)
print('saved', outpath)

# COSMOS ori 304688
# COSMOS selected 62764
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_shifted.feather
# Elapsed:0:08.98,User=23.362,System=7.623,CPU=344.9%.
