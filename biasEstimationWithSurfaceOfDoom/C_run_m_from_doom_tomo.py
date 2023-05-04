# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-05-01 14:15:02
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-03 17:10:03

### example of how to use m_from_doom_func.py 
###### this example is for estimating m for tomographic bins

import numpy as np
import pandas as pd
import statsmodels.api as sm 

from scipy import stats
from astropy.io import fits
from astropy.table import Table

# self-defined module
from m_from_doom_func import mCalFunc_from_doom

# >>>>>>>>>>>>>> I/O

## the surface of doom
inpath_doom = './SurfaceOfDoom/m_surface_ZB0p1_SNR20_R20_skills_v07D7p1_LF_321_kidsPhotometry_gold.csv'
col_m1, col_m2 = 'm1_final', 'm2_final'

## the catalogue needs m
inpath_cata = '/disks/shear10/ssli/K1000CS/LF321_Outputs/K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_goldclasses.cat.A2.feather'
col_ZB, col_SNR, col_R, col_weight = 'Z_B', 'model_SNratio', 'R', 'AlphaRecalC_weight'
col_goldFlag = 'Flag_SOM_Fid_NONE'

## the defined tomo bins
tomo_edges = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]

## where to save the results
outpath = './test_results/m_tomo_K1000_LF_glab_321_v2_A12_gold.csv'

# ### running info
# Number of sources in the data 23401764
# number after gold selection 23401764
# selected objects (weight>0) 23401764
# ZB (0.1, 1.2]: -0.003827599088373934, 0.0038923554387873035
# ZB (0.1, 0.3]: -0.019833358532438103, 0.019752404643884218
# ZB (0.3, 0.5]: -0.024093381810425335, 0.008003787297216272
# ZB (0.5, 0.7]: -0.01558064120531115, 0.007104566559357633
# ZB (0.7, 0.9]: 0.015195799816212304, 0.00558385046438898
# ZB (0.9, 1.2]: 0.031190475958772577, 0.005999995779093551
# results saved to ./test_results/m_tomo_K1000_LF_glab_321_v2_A12_gold.csv
# Elapsed:5:40.48,User=233.948,System=286.546,CPU=152.8%.

# >>>>>>>>>>>>>>>> workhorse

# load data catalogue
file_type = inpath_cata[-3:]
if file_type == 'csv':
    cata_data = pd.read_csv(inpath_cata)
elif file_type == 'her':
    cata_data = pd.read_feather(inpath_cata)
elif file_type == 'its':
    with fits.open(inpath_cata) as hdul:
        cata_data = hdul[1].data
else:
    raise Exception(f'Not supported input file type! {inpath_cata}')
print('Number of sources in the data', len(cata_data))
### select gold class
if col_goldFlag is not None:
    cata_data = cata_data[(cata_data[col_goldFlag].values>0)]
    cata_data.reset_index(drop=True, inplace=True)
    print('number after gold selection', len(cata_data))
### select weight
cata_data = cata_data[cata_data[col_weight]>0]
cata_data.reset_index(drop=True, inplace=True)
print('selected objects (weight>0)', len(cata_data))

# load surface of doom
cata_doom = pd.read_csv(inpath_doom)

# calculate m
m_res = pd.DataFrame(-999., 
                    index = np.arange(len(tomo_edges)), 
                    columns = ['m1', 'm2', 'm1_err', 'm2_err',
                                'Z_B_min', 'Z_B_max', 'Nwei'])
## for the whole sample
#### only for those within the bin
zbin_min = tomo_edges[0]
zbin_max = tomo_edges[-1]
cata_selec = cata_data[(cata_data[col_ZB]>zbin_min)&(cata_data[col_ZB]<=zbin_max)]
m1, m2, m1_err, m2_err = mCalFunc_from_doom(cata_selec, cata_doom, 
                                col_ZB, col_SNR, col_R, col_weight,
                                col_m1, col_m2)
print(f'ZB ({zbin_min}, {zbin_max}]: {(m1+m2)/2.}, {(m1_err+m2_err)/2.}')
m_res.loc[0, 'm1'] = m1
m_res.loc[0, 'm2'] = m2
m_res.loc[0, 'm1_err'] = m1_err
m_res.loc[0, 'm2_err'] = m2_err
m_res.loc[0, 'Z_B_min'] = zbin_min
m_res.loc[0, 'Z_B_max'] = zbin_max
m_res.loc[0, 'Nwei'] = np.sum(cata_selec[col_weight].values)
del cata_selec
## for tomo bins
for i_zbin in range(len(tomo_edges)-1):
    zbin_min = tomo_edges[i_zbin]
    zbin_max = tomo_edges[i_zbin+1]
    ## select sample
    cata_selec = cata_data[(cata_data[col_ZB]>zbin_min)&(cata_data[col_ZB]<=zbin_max)]
    ## calculate m from surface of doom
    m1, m2, m1_err, m2_err = mCalFunc_from_doom(cata_selec, cata_doom, 
                col_ZB, col_SNR, col_R, col_weight,
                col_m1, col_m2)
    print(f'ZB ({zbin_min}, {zbin_max}]: {(m1+m2)/2.}, {(m1_err+m2_err)/2.}')
    ## save the results
    m_res.loc[(i_zbin+1), 'm1'] = m1
    m_res.loc[(i_zbin+1), 'm2'] = m2
    m_res.loc[(i_zbin+1), 'm1_err'] = m1_err
    m_res.loc[(i_zbin+1), 'm2_err'] = m2_err
    m_res.loc[(i_zbin+1), 'Z_B_min'] = zbin_min
    m_res.loc[(i_zbin+1), 'Z_B_max'] = zbin_max
    m_res.loc[(i_zbin+1), 'Nwei'] = np.sum(cata_selec[col_weight].values)
    del cata_selec

m_res.to_csv(outpath, index=False, float_format='%.6f')
print(f'results saved to {outpath}')
