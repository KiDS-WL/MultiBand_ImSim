# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-07-23 11:25:09
# @Last Modified by:   lshuns
# @Last Modified time: 2022-08-22 14:34:10

### prepare input catalogue by shifting galaxy properties

import pandas as pd 
import numpy as np 

import os 

########## original catalogue
inpath = '/disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0.feather'
cata0 = pd.read_feather(inpath)
# print(cata0.columns)

########## shift
shift_cols = ['Re_arcsec', 'shape/sersic_n', 'BA']
limits = [[1e-2, 10.], [0.5, 6.0], [0.05, 1.0]]
save_names = ['size', 'n', 'q']

# # 1 sigma
# shift_amounts = [0.05, 0.1, 0.05]

# # 2 sigma
# shift_amounts = [0.1, 0.2, 0.1]

# 3 sigma
shift_amounts = [0.15, 0.3, 0.15]

for i_col, shift_col in enumerate(shift_cols):

    save_name = save_names[i_col]
    shift_amount = shift_amounts[i_col]
    limit = limits[i_col]

    save_amount = shift_amount*100

    print('>>>>>>>>>>> for', save_name)

    ## increase
    outpath = inpath.replace('.feather', f'_{save_name}U{save_amount:.0f}.feather')
    cata_tmp = cata0.copy()
    cata_tmp.loc[:, shift_col] *= (1+shift_amount)
    ### limits
    cata_tmp.loc[:, shift_col] = np.where(cata_tmp[shift_col].values<limit[0], 
                                                limit[0], cata_tmp[shift_col].values)
    cata_tmp.loc[:, shift_col] = np.where(cata_tmp[shift_col].values>limit[1], 
                                                limit[1], cata_tmp[shift_col].values)
    ### save
    print('>>> ', cata0[shift_col].values)
    print('>>> ', cata_tmp[shift_col].values)
    print('+++ min, max', np.min(cata_tmp[shift_col].values), np.max(cata_tmp[shift_col].values))
    cata_tmp.to_feather(outpath)
    del cata_tmp
    print('saved', outpath)
    ## decrease
    outpath = inpath.replace('.feather', f'_{save_name}D{save_amount:.0f}.feather')
    cata_tmp = cata0.copy()
    cata_tmp.loc[:, shift_col] *= (1 - shift_amount)
    ### limits
    cata_tmp.loc[:, shift_col] = np.where(cata_tmp[shift_col].values<limit[0], 
                                                limit[0], cata_tmp[shift_col].values)
    cata_tmp.loc[:, shift_col] = np.where(cata_tmp[shift_col].values>limit[1], 
                                                limit[1], cata_tmp[shift_col].values)
    ### save
    print('>>> ', cata_tmp[shift_col].values)
    print('+++ min, max', np.min(cata_tmp[shift_col].values), np.max(cata_tmp[shift_col].values))
    cata_tmp.to_feather(outpath)
    print('saved', outpath)

############################ 1 sigma
# >>>>>>>>>>> for size
# >>>  [0.5990288  0.2574687  0.32050198 ... 1.06220223 0.34539227 6.09180839]
# >>>  [0.62898024 0.27034214 0.33652708 ... 1.11531234 0.36266188 6.39639881]
# +++ min, max 0.02872157734254494 10.0
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_sizeU5.feather
# >>>  [0.56907736 0.24459527 0.30447688 ... 1.00909212 0.32812265 5.78721797]
# +++ min, max 0.02598618902420732 9.439208648742284
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_sizeD5.feather
# >>>>>>>>>>> for n
# >>>  [0.5        1.08866797 1.06168242 ... 1.29778034 1.17381028 0.91325004]
# >>>  [0.55       1.19753477 1.16785066 ... 1.42755837 1.2911913  1.00457505]
# +++ min, max 0.55 6.0
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_nU10.feather
# >>>  [0.5        0.97980117 0.95551417 ... 1.1680023  1.05642925 0.82192504]
# +++ min, max 0.5 5.4
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_nD10.feather
# >>>>>>>>>>> for q
# >>>  [0.46000001 0.55000001 0.28       ... 0.66185985 0.34202672 0.3855206 ]
# >>>  [0.48300001 0.57750001 0.294      ... 0.69495284 0.35912806 0.40479663]
# +++ min, max 0.05 1.0
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_qU5.feather
# >>>  [0.43700001 0.52250001 0.266      ... 0.62876686 0.32492539 0.36624457]
# +++ min, max 0.05 1.0
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_qD5.feather
# Elapsed:0:37.65,User=78.735,System=34.857,CPU=301.6%.

# ############################# 2 sigma
# >>>>>>>>>>> for size
# >>>  [0.5990288  0.2574687  0.32050198 ... 1.06220223 0.34539227 6.09180839]
# >>>  [0.65893168 0.28321557 0.35255217 ... 1.16842245 0.37993149 6.70098923]
# +++ min, max 0.030089271501713744 10.0
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_sizeU10.feather
# >>>  [0.53912592 0.23172183 0.28845178 ... 0.95598201 0.31085304 5.48262755]
# +++ min, max 0.024618494865038518 8.942408193545322
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_sizeD10.feather
# >>>>>>>>>>> for n
# >>>  [0.5        1.08866797 1.06168242 ... 1.29778034 1.17381028 0.91325004]
# >>>  [0.6        1.30640156 1.2740189  ... 1.5573364  1.40857233 1.09590005]
# +++ min, max 0.6 6.0
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_nU20.feather
# >>>  [0.5        0.87093437 0.84934593 ... 1.03822427 0.93904822 0.73060003]
# +++ min, max 0.5 4.800000000000001
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_nD20.feather
# >>>>>>>>>>> for q
# >>>  [0.46000001 0.55000001 0.28       ... 0.66185985 0.34202672 0.3855206 ]
# >>>  [0.50600001 0.60500001 0.308      ... 0.72804583 0.3762294  0.42407266]
# +++ min, max 0.05 1.0
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_qU10.feather
# >>>  [0.41400001 0.49500001 0.252      ... 0.59567386 0.30782405 0.34696854]
# +++ min, max 0.05 1.0
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_qD10.feather
# Elapsed:1:18.83,User=48.925,System=44.170,CPU=118.0%.

# ############################## 3 sigma
# >>>>>>>>>>> for size
# >>>  [0.5990288  0.2574687  0.32050198 ... 1.06220223 0.34539227 6.09180839]
# >>>  [0.68888312 0.29608901 0.36857727 ... 1.22153257 0.39720111 7.00557965]
# +++ min, max 0.03145696566088255 10.0
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_sizeU15.feather
# >>>  [0.50917448 0.2188484  0.27242668 ... 0.9028719  0.29358343 5.17803713]
# +++ min, max 0.02325080070586971 8.44560773834836
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_sizeD15.feather
# >>>>>>>>>>> for n
# >>>  [0.5        1.08866797 1.06168242 ... 1.29778034 1.17381028 0.91325004]
# >>>  [0.65       1.41526836 1.38018714 ... 1.68711444 1.52595336 1.18722505]
# +++ min, max 0.65 6.0
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_nU30.feather
# >>>  [0.5        0.76206758 0.74317769 ... 0.90844624 0.82166719 0.63927503]
# +++ min, max 0.5 4.199999999999999
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_nD30.feather
# >>>>>>>>>>> for q
# >>>  [0.46000001 0.55000001 0.28       ... 0.66185985 0.34202672 0.3855206 ]
# >>>  [0.52900001 0.63250001 0.322      ... 0.76113883 0.39333073 0.44334869]
# +++ min, max 0.05 1.0
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_qU15.feather
# >>>  [0.39100001 0.46750001 0.238      ... 0.56258087 0.29072271 0.32769251]
# +++ min, max 0.05 1.0
# saved /disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_qD15.feather
# Elapsed:1:13.56,User=50.312,System=43.322,CPU=127.2%.